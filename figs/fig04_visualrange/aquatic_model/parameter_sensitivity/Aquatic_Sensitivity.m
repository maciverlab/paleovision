function Aquatic_Sensitivity
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visual range calculations based on alternative parameter values
%%
%% Title                : A massive increase in visual range preceded the origin of terrestrial vertebrates
%% Authors              : Ugurcan Mugan, Malcolm A. MacIver
%% Authors' Affiliation : Northwestern University
%% This work is licensed under the Creative Commons Attribution 4.0 International License. 
%% To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
%% January 2017
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global BIGEYEROOT
%% Intialize variables and sensitivity variables
    run Parameters.m
    load('Parameters.mat')
    DrangeAquatic=linspace(minpupil,maxpupil,30); 
    depth=8; epsilon=5e-4;

    sensitivityParams.XAquatic=0.08;
    sensitivityParams.qAquatic=0.5;
    sensitivityParams.DtAquatic1=0.0116;
    sensitivityParams.DtAquatic2=1.6;
    sensitivityParams.d=10e-6;
    sensitivityParams.M=3;
    aquaticequiv={'XAquatic','qAquatic','DtAquatic','DtAquatic','d','M'};
    parameters={'XAquatic','qAquatic','DtAquatic1', 'DtAquatic2','d','M'};
    
     conditions={'Daylight'};
     viewing={'Upward','Horizontal'};
%% Given condition, viewing solve for range for all pupil diameter 
    visualRangeAquatic=zeros(length(DrangeAquatic),2,length(parameters),length(conditions));
    actRangeAquatic=visualRangeAquatic;
    rInit=[5,3];
    for s=1:length(parameters)
           eval(sprintf('%s=%d',aquaticequiv{s},sensitivityParams.(parameters{s})));
        for c=1:length(conditions);
            a=aAquatic.(conditions{c}); 
            b=bAquatic.(conditions{c});
            C0=C0Aquatic.(conditions{c});
            pAbsorb=pAbsorbAquatic.(conditions{c});
            for j=1:length(viewing);
                K=KAquatic.(conditions{c}).(viewing{j}); 
                L=LAquatic.(conditions{c}).(viewing{j}); 
                r=rInit(c);
                for i=1:length(DrangeAquatic)
                    D=DrangeAquatic(i);
                    delta=10^(floor(log10(r))-4);

                    P=10;
                    while abs(P-1)>epsilon
                        P=firingThresh(depth,lambda,...
                            pAbsorb,a,b,K,L,...
                            r,D,XAquatic,DtAquatic,qAquatic,d,k,len,T,M,R,C0);
                        if P>1
                            r=r-delta;
                        else
                            r=r+delta;
                        end
                        clc;
                        fprintf('parameter: %s\npupil diameter: %f\nrange: %f\nerror: %f\n',...
                            parameters{s},D,r,abs(P-1));
                    end
                visualRangeAquatic(i,j,s,c)=r;
                end
            end
        end
        load('Parameters')
    end

save([BIGEYEROOT 'fig04_visualrange/aquatic_model/parameter_sensitivity/meteoAquaticSensitivityParameters.mat'],'visualRangeAquatic','DrangeAquatic');
%% Given condition and viewing calculate contrast threshold for a luminance and pupil diameter    
    load('meteoAquaticSensitivityParameters.mat');
    actRangeAquatic=zeros(size(visualRangeAquatic,1),size(visualRangeAquatic,2),size(visualRangeAquatic,3),length(conditions));
    epsilon=5e-2;
    for c=1:length(conditions)     
        a=aAquatic.(conditions{c}); 
        b=bAquatic.(conditions{c}); 
        C0=C0Aquatic.(conditions{c});
        for j=1:length(viewing);
            K=KAquatic.(conditions{c}).(viewing{j})(:,depth);
            L=LAquatic.(conditions{c}).(viewing{j})(:,depth);
            [~,index]=max(L);
            lambdaMax=lambda(index);
            B=BAquatic.(conditions{c}).(viewing{j})(depth); 
            for s=1:length(parameters)
                afunc=@(l) interp1(lambda,a,l,'pchip'); 
                bfunc=@(l) interp1(lambda,b,l,'pchip');
                Kfunc=@(l) interp1(lambda,K,l,'pchip');
                for i=1:length(DrangeAquatic)
                    D=DrangeAquatic(i);
                    tempRange=visualRangeAquatic(i,j,s,c);

                    %Apparent contrast for horizontal and upward viewing
                    Cr(i)=C0*exp((Kfunc(lambdaMax)-afunc(lambdaMax)-bfunc(lambdaMax)).*tempRange);                
                    %angular size of object horizontal and upward in mrad
                    angularSize(i)=atan(T/(2*tempRange))*2*10^3;

                    Kt(i)=liminalContrast(A,B,angularSize(i));

                    if (10^(Kt(i))<abs(Cr(i)))
                        actRangeAquatic(i,j,s,c)=tempRange;
                    else    
                        tempVisualRange=linspace(tempRange,tempRange/10,tempRange*1000);
                        ind=1;
                        while(10^(Kt(i))-abs(Cr(i)))>epsilon
                            tempRange=tempVisualRange(ind);
                            angularSize(i)=atan(T/(2*tempRange))*2*10^3;
                            Cr(i)=C0*exp((Kfunc(lambdaMax)-afunc(lambdaMax)-bfunc(lambdaMax)).*tempRange);
                            Kt(i)=liminalContrast(A,B,angularSize(i));
                            ind=ind+1;
                            clc;
                            fprintf('parameter: %s\npupil diameter: %f\nrange: %f\nerror: %f\n',...
                                parameters{s},D,tempRange,10^(Kt(i))-abs(Cr(i)));
                        end
                        actRangeAquatic(i,j,s,c)=tempRange;
                    end
                end
            end
        end
    end
    visualRangeAquatic=actRangeAquatic;   
save([BIGEYEROOT 'fig04_visualrange/aquatic_model/parameter_sensitivity/visibilityAquaticSensitivityParameters.mat'],'visualRangeAquatic','DrangeAquatic')

function Kt = liminalContrast(A,L,angularSize)
%Equations should be the same as in the aquatic case, quick annotation
%Define functions
    logLa=@(L) log10(L);
    iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();
    deltat=@(L) iif( -2<logLa(L) && logLa(L)<=4, @() (1.22*((555e-9)/A)*10^3),...
        logLa(L)<=-2, @() (1.22*((507e-9)/A))*10^3); %Supplementary Appendix, pg. 5

    Psi=@(L) (logLa(L)+3.5)/2.75;%Supplementary Appendix, Eq. 10 part 2 pg. 5
    %Equations from Aksyutov
    logdeltarStar= @(L) iif(logLa(L)>1.535, @() 0.928,...
        -0.75<=logLa(L) && logLa(L)<=1.535, @() 1-(0.5*(Psi(L)^(-3.2))),...
        -2.089<=logLa(L) && logLa(L)<-0.75, @()-0.28+0.1*Psi(L)^(3.2),...
        logLa(L)<-2.089, @() 0.14*logLa(L)+0.442); 
    X=@(deltar,L) -(-log10(deltar)+logdeltarStar(L));
    Xm=@(L) iif(logLa(L)>1.5, @() 1.48-0.48*logLa(10^1.5),...
        logLa(L)>0.11 && logLa(L)<=1.5, @() 1.48-0.48*logLa(L),...
        logLa(L)>=-1.46 && logLa(L)<=0.11, @() 1.40-0.78*logLa(L),...
        logLa(L)<-1.46, @() 1.40-0.78*logLa(10^-1.46)); 
    Ym=@(L) 0.75*Xm(L)-0.32; 
    Xr=@(deltar,L) X(deltar,L)./Xm(L); 
    Yr=@(deltar,L) iif(Xr(deltar,L)<=0.56, @() ((1.6+0.23*logLa(L)).*(Xr(deltar,L)-0.05)).^(0.5),...
        Xr(deltar,L)>0.56, @() (1-((0.42-0.26*logLa(L)).*(1-Xr(deltar,L))))^(0.5)); 
    Y=@(deltar,L) Yr(deltar,L).*Ym(L);

 % Calculate contrast threshold
    angularResolution=deltat(L);
    deltarparam=angularSize/angularResolution;
    logStardeltar=logdeltarStar(L);
    if log10(deltarparam)<=logStardeltar
        Kt=-2*log10(deltarparam); %Supplementary Appendix, pg. 5
    else %Conditions from Aksyutov
        Xparam=X(deltarparam,L);
        Xmparam=Xm(L);
        if Xparam>0 && Xparam<=0.20 
            Yval=1.45*Xparam;
            Kt =-2*logStardeltar-Yval; 
        elseif Xparam>0.20
            Yval=Y(deltarparam,L); 
            Kt=-2*logStardeltar-Yval; 
        end
    end        


function  solution=firingThresh(depth,lambda,photoreceptorAbsorption,a,b,KAll,LAll,r,A,X,Dt,q,d,k,len,T,M,R,C0)
    lambda1=lambda(1); lambda2=lambda(end);
    K=KAll(:,depth);
    L=LAll(:,depth);

    alphaInterp=@(l) interp1(lambda,photoreceptorAbsorption,l,'pchip');
    aInterp=@(l) interp1(lambda,a,l,'pchip');
    bInterp=@(l) interp1(lambda,b,l,'pchip');
    LInterp=@(l) interp1(lambda,L,l,'pchip');
    KInterp=@(l) interp1(lambda,K,l,'pchip');

    Nfalse=((T*M*A)/(2*r*d))^2*X*Dt;  %Supplementary Appendix, pg 3 for explanation.
    %Supplementary Appendix, Eq 7,8 pg. 5
    %Integrand
    RhFunc=@(l) LInterp(l).*l.*(1-exp(-k*alphaInterp(l)*len));
    RoFunc=@(l) LInterp(l).*l.*(1-exp(-k*alphaInterp(l)*len)).*(1+(C0*exp((KInterp(l)-aInterp(l)-bInterp(l))*r)));    
    %Integral
    Rh=integral(RhFunc,lambda1,lambda2);
    Ro=integral(RoFunc,lambda1,lambda2);
    %pg 7, Eq 7 for h, and O
    Nh=((pi/4)^2)*(A^2)*((T/r)^2)*q*Dt*Rh;
    No=((pi/4)^2)*(A^2)*((T/r)^2)*q*Dt*Ro;
    %Supplementary Appendix, Eq. 9 pg. 5
    solution=(R*sqrt(No+Nh+2*Nfalse))/(abs(No-Nh));    