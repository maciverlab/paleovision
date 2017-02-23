function Aquatic_daylightContrastRange
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation of visual range for finned and digited pupil diameters based on optical properties and contrast threshold with varying contrast
%% CONDITION: DAYLIGHT
%%
%% Title                : A massive increase in visual range preceded the origin of terrestrial vertebrates
%% Authors              : Ugurcan Mugan, Malcolm A. MacIver
%% Authors' Affiliation : Northwestern University
%% This work is licensed under the Creative Commons Attribution 4.0 International License. 
%% To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
%% January 2017
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global BIGEYEROOT
%% Initialize Variables
    run Parameters.m
    load('Parameters.mat')
    load OM_TF_ST.mat
    load FinnedDigitedOrbitLength.mat
    
    pupil_TF = [mean(noElpistoOrb)-std(noElpistoOrb) mean(noElpistoOrb)+std(noElpistoOrb)].*0.449;
    pupil_ST = [mean(noSecAqOrb)-std(noSecAqOrb) mean(noSecAqOrb)+std(noSecAqOrb)].*0.449;
    finnedpupil=mean(noElpistoOrb)*.449;
    digitedpupil=mean(noSecAqOrb)*.449;
    DrangeAquatic=sort([finnedpupil,digitedpupil]*1e-3);
    C0Range=linspace(-1,4,20);
    viewing={'Upward','Horizontal'};
    
    tol=5e-5; depth=8;
    rInit=[5.7,3.7];
    visualRangeAquatic=zeros(length(C0Range),length(DrangeAquatic),length(viewing));    
    
    aValue=aAquatic.Daylight; 
    bValue=bAquatic.Daylight;
    aFunc=@(l) interp1(lambda,aValue,l,'pchip'); 
    bFunc=@(l) interp1(lambda,bValue,l,'pchip');
    pAbsorb=pAbsorbAquatic.Daylight;

%% Calculate visual range based on optical properties
    for v=length(viewing)
        KValue=KAquatic.Daylight.(viewing{v})(:,depth);
        LValue=LAquatic.Daylight.(viewing{v})(:,depth);
        BValue=BAquatic.Daylight.(viewing{v})(depth);
        [~,indexMax]=max(LValue);
        lambdaMax=lambda(indexMax);               
        KFunc=@(l) interp1(lambda,KValue,l,'pchip');
        r=rInit(v);
        for j=1:length(DrangeAquatic);
            D=DrangeAquatic(j);
            for i=1:length(C0Range);
                C0=C0Range(i);
                delta=10^(floor(log10(r))-5);
             
                P=10;
                while abs(P-1)>tol
                    P=firingThresh(lambda,...
                         pAbsorb,aValue,bValue,KValue,LValue,...
                         r,D,XAquatic,DtAquatic,qAquatic,d,k,len,T,M,R,C0);
                     if P>1
                         r=r-delta;
                     else
                         r=r+delta;
                     end
                     clc;
                     fprintf('pupil diameter: %f\ncontrast: %f\nrange: %f\nerror: %f\n',...
                         D, C0, r,abs(P-1) );
                end            
                tempRange=r;                
%% Calculate range based on contrast threshold    
                Cr(i)=C0*exp((KFunc(lambdaMax)-aFunc(lambdaMax)-bFunc(lambdaMax)).*tempRange);
                angularSize(i)=atan(T/(2*tempRange))*2*10^3;                 
                Kt(i)=liminalContrast(D,BValue,angularSize(i));

                if (10^(Kt(i))<abs(Cr(i)))
                    visualRangeAquatic(i,j,v)=tempRange;
                else    
                    tempVisualRange=linspace(tempRange,floor(tempRange/10),ceil(tempRange)*1000);
                    ind=1;
                    while(10^(Kt(i)) - abs(Cr(i)))>5e-2
                        tempRange=tempVisualRange(ind);
                        angularSize(i)=atan(T/(2*tempRange))*2*10^3;
                        Cr(i)=C0*exp((KFunc(lambda_down)-aFunc(lambda_down)-bFunc(lambda_down)).*tempRange);
                        Kt(i)=liminalContrast(D,BValue,angularSize(i));
                        ind=ind+1;      
                        clc;
                        fprintf('pupil diameter: %f\ncontrast: %f\nrange: %f\nerror: %f\n',...
                            D, C0, tempRange, (10^(Kt(i))-abs(Cr(i))));
                    end
                    visualRangeAquatic(i,j,v)=tempRange;
                end
            end
        end
    end
    
save([BIGEYEROOT, 'figExt07_contrast/aquatic_model/daylightVisibilityAquaticContrast.mat'],'C0Range','DrangeAquatic','visualRangeAquatic');
            
function  solution=firingThresh(lambda,photoreceptorAbsorption,a,b,K,L,r,D,X,Dt,q,d,k,len,T,M,R,C0)
    lambda1=lambda(1); lambda2=lambda(end);
    alphaInterp=@(l) interp1(lambda,photoreceptorAbsorption,l,'pchip');
    aInterp=@(l) interp1(lambda,a,l,'pchip');
    bInterp=@(l) interp1(lambda,b,l,'pchip');
    LInterp=@(l) interp1(lambda,L,l,'pchip');
    KInterp=@(l) interp1(lambda,K,l,'pchip');
    
    Nfalse=((T*M*D)/(2*r*d))^2*X*Dt;%Supplementary Appendix, pg 3 for explanation.
    %Supplementary Appendix, Eq 7,8 pg. 5
    %Integrand
    RhFunc=@(l) LInterp(l).*l.*(1-exp(-k*alphaInterp(l)*len));
    RoFunc=@(l) LInterp(l).*l.*(1-exp(-k*alphaInterp(l)*len)).*(1+(C0*exp((KInterp(l)-aInterp(l)-bInterp(l))*r)));
    %Integral
    Rh=integral(RhFunc,lambda1,lambda2);
    Ro=integral(RoFunc,lambda1,lambda2);    
    Nh=((pi/4)^2)*(D^2)*((T/r)^2)*q*Dt*Rh;
    No=((pi/4)^2)*(D^2)*((T/r)^2)*q*Dt*Ro;
    %Supplementary Appendix, Eq. 9 pg. 5
    solution=(R*sqrt(No+Nh+2*Nfalse))/(abs(No-Nh));  
    
function Kt = liminalContrast(A,L,angularSize)
%Define functions
    logLa=@(L) log10(L);
    iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();
    deltat=@(L) iif( -2<logLa(L) && logLa(L)<=4, @() (1.22*((555e-9)/A)*10^3),...
        logLa(L)<=-2, @() (1.22*((507e-9)/A))*10^3); %Supplementary Appendix, pg. 5

    Psi=@(L) (logLa(L)+3.5)/2.75; %Supplementary Appendix, Eq. 10 part 1 pg. 5
    logdeltarStar= @(L) iif(logLa(L)>1.535, @() 0.928,...
        -0.75<=logLa(L) && logLa(L)<=1.535, @() 1-(0.5*(Psi(L)^(-3.2))),...
        -2.089<=logLa(L) && logLa(L)<-0.75, @()-0.28+0.1*Psi(L)^(3.2),...
        logLa(L)<-2.089, @() 0.14*logLa(L)+0.442);%Supplementary Appendix, Eq. 10 part 2 pg. 5
    %Equations from Aksyutov
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
