function Aerial_contrastRangeRelation
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
%% Intialize variables
    run Parameters.m
    load('Parameters.mat')
    load OM_TF_ST.mat
    load FinnedDigitedOrbitLength.mat
    
    pupil_TF = [mean(noElpistoOrb)-std(noElpistoOrb) mean(noElpistoOrb)+std(noElpistoOrb)].*0.449;
    pupil_ST = [mean(noSecAqOrb)-std(noSecAqOrb) mean(noSecAqOrb)+std(noSecAqOrb)].*0.449;
    finnedpupil=mean(noElpistoOrb)*.449;
    digitedpupil=mean(noSecAqOrb)*.449;
    DrangeAerial=sort([finnedpupil,digitedpupil]*1e-3);
    C0Range=linspace(-1,4,20);   
    epsilon=5e-5; r=1.4e4;
    
    Wlambdaylambda=csvread('Wlambda.csv'); %Spectral radiance emittence, Table from Middleton
    WlambdaylambdaInterp= @(lambda) interp1(Wlambdaylambda(:,1),Wlambdaylambda(:,2),lambda,'pchip'); %value checked with mathematica
    
    lambda1=0.4; lambda2=0.7;%visible wavelength
    %extinction coeff from Middleton (originally in 1/km coverted to 1/m)
    sigma=@(lambda) ((1.1e-3*lambda.^(-4))+(0.008*lambda.^(-2.09)))/(1e3); 
    
    Bh=BAerial.Daylight*integral(WlambdaylambdaInterp,lambda1,lambda2); %Supplementary Appendix Materials and Methods, pg.2
    %Intensitiy Parameter coefficient, supplement pg 5 line 550
    %Full equation to get Rh is on line 549 Rh=beta*B_h,o(lambda)
    %horizon space-light photons/m^2 s sr, value checked with mathematica
    Rh=((1.31e3)/0.89)*Bh*(1e6)^2; %value checked with mathematica
    
%% Calculate range based on optical properties
    Dt=DtAerial.Daylight; F=FAerial.Daylight;
    X=XAerial; q=qAerial.Daylight;

    visualRangeAerial=zeros(length(C0Range),length(DrangeAerial));        
    for i=1:length(DrangeAerial)
        D=DrangeAerial(i);
        for c=1:length(C0Range)
            delta=10^(floor(log10(r))-5);
            C0=C0Range(c);                       
            P=10;
            while abs(P-1)>=epsilon
                P=firingThreshRange(...
                     WlambdaylambdaInterp,D,r,C0,BAerial.Daylight,...
                     T,F,Rh,Dt,q,k,len,X,d,R);
                if P>1
                    r=r-delta;
                else
                    r=r+delta;
                end
                clc;
                fprintf('pupil diameter: %f\ncontrast: %f\nrange: %f\nerror: %f\n',...
                    D,C0,r,abs(P-1));
            end
            visualRangeAerial(c,i)=r;
            tempRange=visualRangeAerial(c,i);

%% Calculate range based on contrast threshold
            CrFunc=@(lambda) exp(-sigma(lambda).*tempRange);
            Cr= C0*integral(CrFunc,lambda1,lambda2);

            angularSize=atan(T/(2*tempRange))*2*10^3;
            Kt=liminalContrast(D,BAerial.Daylight,angularSize);

            if 10^(Kt) <= abs(Cr)
                 visualRangeAerial(c,i)=tempRange;
            else
                delta=delta*1e1;
                while(10^(Kt)>abs(Cr))
                    tempRange=tempRange-delta;
                    angularSize=atan(T/(2*tempRange))*2*10^3;
                    Cr=C0*integral(CrFunc,lambda1,lambda2);
                    Kt=liminalContrast(D,BAerial.Daylight,angularSize);
                    clc
                    fprintf('pupil diameter: %f\ncontrast: %f\nrange: %f\nerror: %f\n',...
                    D,C0,tempRange,10^(Kt)>abs(Cr));
                end
                visualRangeAerial(c,i)=tempRange;
            end
        end
    end
save([BIGEYEROOT,'figExt07_contrast/aerial_model/daylightVisibilityAerialContrast.mat'],'visualRangeAerial','C0Range','DrangeAerial');                

    
function  solution=firingThreshRange(WHandle,D,r,C0,B,T,F,Rh,Dt,q,k,len,X,d,R)
    lambda1=0.4; lambda2=0.7; %visible wavelength
    %extinction coeff from Middleton (originally in 1/km coverted to 1/m)
    sigma=@(lambda) ((1.1e-3*lambda.^(-4))+(0.008*lambda.^(-2.09)))/(1e3);
    
    Nh=(pi/4)^2*(T/r)^2*D^2*Rh*Dt*q*((k*len)/(2.3+(k*len))); %Supplementary Appendix pg. 3
    Nfalse=((T*F*D)/(r*d))^2*X*Dt;

    Bofunc=@(lambda) WHandle(lambda).*(1+(C0.*exp(-sigma(lambda).*r)));
    Bo= B*integral(Bofunc,lambda1,lambda2);
    Ro=((1.31e3)/0.89)*Bo*(1e6)^2; %Supplementary Appendix, Eq. 3 pg. 2
    
    No=(pi/4)^2*(T/r)^2*D^2*Ro*Dt*q*((k*len)/(2.3+(k*len))); %Supplementary Appendix, combination of Eq. 5 and Rh, first part of Eq. 6

    solution=(R*sqrt(No+Nh+2*Nfalse))/(abs(No-Nh));

function Kt = liminalContrast(A,L,angularSize)
%% Function definitions, section: contrast threshold as a function of pupil diameter
    %delta0 is labeled deltar
    %initial definitions -------
    logLa=@(L) log10(L);
    iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();
    deltat=@(L) iif( -2<logLa(L) && logLa(L)<=4, @() (1.22*((555e-9)/A)*10^3),...
    logLa(L)<=-2, @() (1.22*((507e-9)/A))*10^3); %Supplementary Appendix, pg. 5
    %----------
    %Equations adapted from Aksyutov 
    Psi=@(L) (logLa(L)+3.5)/2.75;  %Supplementary Appendix, pg.5 Eq. 10 part 1
    logdeltarStar= @(L) iif(logLa(L)>1.535, @() 0.928,...
        -0.75<=logLa(L) && logLa(L)<=1.535, @() 1-(0.5*(Psi(L)^(-3.2))),...
        -2.089<=logLa(L) && logLa(L)<-0.75, @()-0.28+0.1*Psi(L)^(3.2),...
        logLa(L)<-2.089, @() 0.14*logLa(L)+0.442);%Supplementary Appendix, pg. 5 Eq. 10 part 2
    %Not provided in supplement as seen in Aksyutov 2002
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
    
    %% Calculate contrast threshold
    angularResolution=deltat(L);
    deltarparam=angularSize/angularResolution;
    logStardeltar=logdeltarStar(L);
    if log10(deltarparam)<=logStardeltar
        Kt=-2*log10(deltarparam); %condition: Supplementary Appendix, pg. 5
    else %Adapted from Aksyutov
        Xparam=X(deltarparam,L);
        if Xparam>0 && Xparam<=0.20
            Yval=1.45*Xparam;
            Kt =-2*logStardeltar-Yval;
        elseif Xparam>0.20
            Yval=Y(deltarparam,L);
            Kt=-2*logStardeltar-Yval;
        end
    end         