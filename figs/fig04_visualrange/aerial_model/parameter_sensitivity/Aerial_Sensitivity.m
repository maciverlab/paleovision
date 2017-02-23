function Aerial_Sensitivity
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visual range based on alternative parameters, except contrast threshold
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
    Wlambdaylambda=csvread('Wlambda.csv'); %Spectral radiance emittence, Table from Middleton
    WlambdaylambdaInterp= @(lambda) interp1(Wlambdaylambda(:,1),Wlambdaylambda(:,2),lambda,'pchip');
    
    %extinction coeff from Middleton (originally in 1/km coverted to 1/m)
    sigma=@(lambda) ((1.1e-3*lambda.^(-4))+(0.008*lambda.^(-2.09)))/(1e3);
    lambda1=0.4; lambda2=0.7; %visible range wavelength
    conditions={'Daylight'};
    
%% Alternative Parameter Values
    sensitivityParams.XAerial=0.011;
    sensitivityParams.qAerial=0.3;
    sensitivityParams.DtAerial1=0.0116;
    sensitivityParams.DtAerial2=1.6;
    sensitivityParams.d=10e-6;
    sensitivityParams.FAerial=3/2;
    
    parameters={'XAerial','qAerial','DtAerial1', 'DtAerial2','d','FAerial'};    
    %% Calcualte visual range based on firing threshold
    DrangeAerial=linspace(minpupil,maxpupil,25);
    epsilon=1e-4;
    
    visualRangeAerial=zeros(length(DrangeAerial),length(parameters),length(conditions));
    for c=1:length(conditions)
        %following Middleton:
        Bh=BAerial.(conditions{c})*integral(WlambdaylambdaInterp,lambda1,lambda2); %horizon luminance, %Supplementary Appendix Materials and Methods, pg.2
        %Intensitiy Parameter coefficient, supplement pg 5 line 550
        %Full equation to get Rh is on line 549 Rh=beta*B_h,o(lambda)
        %horizon space-light photons/m^2 s sr, value checked with mathematica
        Rh=((1.31e3)/0.89)*Bh*(1e6)^2; %Supplementary Appendix, pg.3 and Table S1
        C0=C0Aerial.(conditions{c});
        for s=1:length(parameters)
            if strcmp(parameters{s},'XAerial') || strcmp(parameters{s},'d')
                aerialEquiv=parameters{s};
            else
                aerialEquivTemp=strcat(regexp(parameters{s},'[a-z]','ignorecase','match'));
                aerialEquiv='';
                for reg=1:length(aerialEquivTemp); aerialEquiv=strcat(aerialEquiv,aerialEquivTemp{reg}); end
                aerialEquiv=strcat(aerialEquiv,'.',conditions{c});
            end            
            eval(sprintf('%s=%d',aerialEquiv,sensitivityParams.(parameters{s})));
            DtAerialVal=DtAerial.(conditions{c});
            FAerialVal=FAerial.(conditions{c});
            qAerialVal=qAerial.(conditions{c});
            
            r=8000;
            for i=1:length(DrangeAerial)
                D=DrangeAerial(i);
                delta=10^(floor(log10(r))-4);
                
                P=0;
                while(abs(P-1)>epsilon)                           
                    Nh=(pi/4)^2*(T/r)^2*D^2*Rh*DtAerialVal*qAerialVal*...
                        ((k*len)/(2.3+(k*len))); %Supplementary Appendix, combination of Eq. 5 and Rh, first part of Eq. 6
                    if ~strcmp(conditions{c},'Daylight')
                        Nh=(pi/4)^2*(T/r)^2*D^2*Rh*Dt*q*(1-exp(-k*len));
                    end
                    Nfalse=((T*FAerialVal*D)/(r*d))^2*XAerial*DtAerialVal; %Supplementary Appendix pg. 3
                   
                    Bofunc=@(lambda) WlambdaylambdaInterp(lambda).*...
                        (1+(C0.*exp(-sigma(lambda).*r)));
                    Bo= BAerial.(conditions{c})*integral(Bofunc,lambda1,lambda2);                  
                    Ro=((1.31e3)/0.89)*Bo*(1e6)^2; %Supplementary Appendix, Eq. 3 pg. 2

                    No=(pi/4)^2*(T/r)^2*D^2*Ro*DtAerialVal*qAerialVal*...
                        ((k*len)/(2.3+(k*len)));%Supplementary Appendix, combination of Eq. 5 and Rh, first part of Eq. 6
                    if~strcmp(conditions{c},'Daylight')
                         No=(pi/4)^2*(T/r)^2*D^2*Ro*Dt*q*(1-exp(-k*len));
                    end
                    P=(R*sqrt(No+Nh+2*Nfalse))/(abs(No-Nh)); %Supplementary Appendix, Eq. 4 (pg. 2) and Eq. 6 (pg. 3)
                    
                    if P>1
                        r=r-delta;
                    else
                        r=r+delta;
                    end
                    
                    clc;
                    fprintf('condition: %s\nparameter: %s\npupil diameter: %f\nrange: %f\nerror: %f\n',...
                        conditions{c}, parameters{s}, DrangeAerial(i), r, abs(P-1));
                end
                visualRangeAerial(i,s,c)=r;
            end
            load('Parameters.mat')
        end
    end
    save([BIGEYEROOT 'fig04_visualrange/aerial_model/parameter_sensitivity/meteoAerialParameterSensitivity.mat'],'visualRangeAerial','DrangeAerial')
    %% Calculate visual range based on contrast threshold
    load('meteoAerialParameterSensitivity.mat')

    actRangeAerial=zeros(size(visualRangeAerial,1),size(visualRangeAerial,2), size(visualRangeAerial,3));
    tempRangePrevAll=[10 0.01 0.001];
    
    for c=1:length(conditions);
        C0=C0Aerial.(conditions{c});
        Bh=BAerial.(conditions{c});
        for s=1:length(parameters);
            tempRangePrev=tempRangePrevAll(c);
            for i=1:length(DrangeAerial);
                D=DrangeAerial(i);
                tempRange=visualRangeAerial(i,s,c);

                CrFunc=@(lambda) exp(-sigma(lambda).*tempRange);
                Cr(i)= C0*integral(CrFunc,lambda1,lambda2); %Apparent contrast of object

                angularSize(i)=atan(T/(2*tempRange))*2*10^3; %angular size of object in mrad
                Kt(i)=liminalContrast(D,Bh,angularSize(i)); %contrast threshold given angular size and luminance
                if (10^(Kt(i))) <= abs(Cr(i))
                    actRangeAerial(i)=tempRange;
                else
                    tempVisualRange=linspace(tempRange,tempRangePrev,max(visualRangeAerial(:,s,c))*3);
                    j=1;
                    %while apparent contrast is below liminal contrast the object
                    %is invisible, decrease range until apparent contrast is
                    %greater than the liminal contrast
%                     while (10^(Kt(i))-abs(Cr(i))>1e-4)
                        tempRange=tempVisualRange(j); %decrease visual range
                        %Recalculate angular size, contrast threshold and apparent
                        %contrast based on the decreased range
                        angularSize(i)=atan(T/(2*tempRange))*2*10^3;
                        Cr(i)=C0*integral(CrFunc,lambda1,lambda2);
                        Kt(i)=liminalContrast(D,Bh,angularSize(i));
                        actRangeAerial(i,s,c)=tempRange;
                        j=j+1;
                        clc;
                        fprintf('condition: %s\nparameter: %s\npupil diameter: %f\nrange %f\nerror: %f\n',...
                            conditions{c}, parameters{s}, DrangeAerial(i), tempRange,10^(Kt(i))-abs(Cr(i)));
                    end
                    tempRangePrev=tempRange;
                end
            end
        end
    end
        
    visualRangeAerial=actRangeAerial;
save([BIGEYEROOT 'fig04_visualrange/aerial_model/parameter_sensitivity/visibilityAerialParameterSensitivity.mat'], 'visualRangeAerial','DrangeAerial');
       
    function Kt = liminalContrast(A,L,angularSize)

    %% Function definitions, section: contrast threshold as a function of pupil diameter
    %delta0 is labeled deltar
    %initial definitions -------
    logLa=@(L) log10(L);
    iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();
    deltat=@(L) iif( -2<logLa(L) && logLa(L)<=4, @() (1.22*((555e-9)/A)*10^3),...
    logLa(L)<=-2, @() (1.22*((507e-9)/A))*10^3); %supp pg 14 line 1730
    %----------
    %Equations adapted from Aksyutov
    Psi=@(L) (logLa(L)+3.5)/2.75; %Supplementary Appendix, pg.5 Eq. 10 part 1
    
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
    if log10(deltarparam)<=logStardeltar %condition: Supplementary Appendix, pg. 5
        Kt=-2*log10(deltarparam);
    else %Adapted from Aksyutov
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
