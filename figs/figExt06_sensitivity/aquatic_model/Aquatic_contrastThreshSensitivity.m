function Aquatic_contrastThreshSensitivity
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inclusion of contrast threshold for aquatic vision with different water models
%%
%% Title                : A massive increase in visual range preceded the origin of terrestrial vertebrates
%% Authors              : Ugurcan Mugan, Malcolm A. MacIver
%% Authors' Affiliation : Northwestern University
%% This work is licensed under the Creative Commons Attribution 4.0 International License. 
%% To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
%% January 2017
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global BIGEYEROOT
    run ParametersSensitivity.m
    load('ParametersSensitivity.mat');
    if ~exist('meteoAquaticSensitivity.mat','file')==2
        Aquatic_firingThreshSensitivity;
    end
    
    load('meteoAquaticSensitivity.mat');
    
    conditions={'HighTurbidity','Clear','AbsDom','ScatDom'};
    viewing={'Upward','Horizontal'};
    waterDepth=8;
    DrangeAquatic=linspace(minpupil,maxpupil,30);   
    actRangeAquaticSensitivity=0*visualRangeAquaticSensitivity;
    lambda1=lambda(1); lambda2=lambda(length(lambda));
    for c=1:length(conditions)
        aValue=a.(conditions{c}); 
        bValue=b.(conditions{c});
        for v=1:length(viewing)                         
            KValue=K.(conditions{c}).(viewing{v})(:,waterDepth); 
            BValue=BLum.(conditions{c}).(viewing{v})(waterDepth);          
            LValue=L.(conditions{c}).(viewing{v})(:,waterDepth);
            [~,index]=max(LValue); 
            lmax=lambda(index);

            aFunc=@(l) interp1(lambda,aValue,l,'pchip'); 
            bFunc=@(l) interp1(lambda,bValue,l,'pchip');
            KFunc=@(l) interp1(lambda,KValue,l,'pchip');
            for i=1:length(DrangeAquatic)
               D=DrangeAquatic(i);
               tempRange=visualRangeAquaticSensitivity(i,c,v);
               epsilon=5e-2;
                if strcmp(conditions{c},'AbsDom')
                    epsilon=5e-5;
                end
                
                Cr(i)=C0*exp((KFunc(lmax)-aFunc(lmax)-bFunc(lmax)).*tempRange);
                angularSize(i)=atan(T/(2*tempRange))*2*10^3;
                Kt(i)=liminalContrast(D,BValue,angularSize(i));
                if (10^(Kt(i))<=abs(Cr(i)))
                    actRangeAquaticSensitivity(i,c,v)=tempRange;
                else
                    tempVisualRange=linspace(tempRange,tempRange/10,tempRange*100);
                    ind=1;
                    while 10^(Kt(i))-abs(Cr(i))>epsilon
                        tempRange=tempVisualRange(ind);
                        angularSize(i)=atan(T/(2*tempRange))*2*10^3;
                        Cr(i)=C0*exp((KFunc(lmax)-aFunc(lmax)-bFunc(lmax)).*tempRange);
                        Kt(i)=liminalContrast(D,BValue,angularSize(i));
                 
                        ind=ind+1;
                        clc
                        fprintf('condition: %s\nviewing: %s\npupil diameter: %f\nrange: %f\nerror: %f\n',...
                            conditions{c}, viewing{v},D,tempRange,10^(Kt(i))-abs(Cr(i)));
                    end
                    actRangeAquaticSensitivity(i,c,v)=tempRange;
                end
            end
        end
    end
    
    visualRangeAquaticSensitivity=actRangeAquaticSensitivity;   
save([BIGEYEROOT, 'figExt06_sensitivity/aquatic_model/visibilityAquaticSensitivity.mat'],'visualRangeAquaticSensitivity','DrangeAquatic');
                
                
function Kt = liminalContrast(A,L,angularSize)
%Define functions
    logLa=@(L) log10(L);
    iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();
    deltat=@(L) iif( -2<logLa(L) && logLa(L)<=4, @() (1.22*((555e-9)/A)*10^3),...
    logLa(L)<=-2, @() (1.22*((507e-9)/A))*10^3);   %Supplementary Appendix, pg. 5

    Psi=@(L) (logLa(L)+3.5)/2.75;%Supplementary Appendix, Eq. 10 part 1 pg. 5
    logdeltarStar= @(L) iif(logLa(L)>1.535, @() 0.928,...
        -0.75<=logLa(L) && logLa(L)<=1.535, @() 1-(0.5*(Psi(L)^(-3.2))),...
        -2.089<=logLa(L) && logLa(L)<-0.75, @()-0.28+0.1*Psi(L)^(3.2),...
        logLa(L)<-2.089, @() 0.19*logLa(L)+0.442);%Supplementary Appendix, Eq. 10 part 2 pg. 5
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
        Kt=-2*log10(deltarparam);%Supplementary Appendix, pg. 5
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
