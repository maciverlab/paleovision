function Aerial_contrastThreshSensitivity
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visual range based on the y-axis shift in contrast threshold
%%
%% Title                : A massive increase in visual range preceded the origin of terrestrial vertebrates
%% Authors              : Ugurcan Mugan, Malcolm A. MacIver
%% Authors' Affiliation : Northwestern University
%% This work is licensed under the Creative Commons Attribution 4.0 International License. 
%% To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
%% January 2017
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global BIGEYEROOT
    %% Initialize variables
    run Parameters.m
    load('Parameters.mat')
    if ~exist('meteoAerial.mat','file')==2
        Aerial_firingThresh;
    end
    load('meteoAerial.mat')
    
    %Extinction coefficient from Middeleton, converted from 1/km to 1/m
    sigma=@(lambda) ((1.1e-3*lambda.^(-4))+(0.008*lambda.^(-2.09)))/(1e3);
    lambda1=0.4; lambda2=0.7;%visible range wavelength
    conditions={'Daylight'};

    %% Calculate visual range from contrast threshold
    actRangeAerial=zeros(size(visualRangeAerial,1),length(conditions));
    tempRangePrevAll=[10,0.01,0.001];
    for c=1:length(conditions);
        C0=C0Aerial.(conditions{c});   
        tempRangePrev=tempRangePrevAll(c);
        for i=1:length(DrangeAerial);
            D=DrangeAerial(i);
            tempRange=visualRangeAerial(i,c);

            CrFunc=@(lambda) exp(-sigma(lambda).*tempRange); %contrast attenuation
            Cr(i)= C0*integral(CrFunc,lambda1,lambda2); %Apparent contrast of object
            %Apparent contrast value accounting for all wavelengths,
            %C0Aerial.Daylight is the contrast of the object

            angularSize(i)=atan(T/(2*tempRange))*2*10^3; %angular size of object in mrad
            Kt(i)=liminalContrast(D,BAerial.(conditions{c}),angularSize(i)); %contrast threshold given angular size and luminance
            if (10^(Kt(i))*1.15) <= abs(Cr(i))
                actRangeAerial(i)=tempRange;
            else
                tempVisualRange=linspace(tempRange,tempRangePrev,max(visualRangeAerial(:,c))*3);
                j=1;
                %while apparent contrast is below liminal contrast the object
                %is invisible, decrease range until apparent contrast is
                %greater than the liminal contrast
                while((10^(Kt(i)))*1.15 > abs(Cr(i)))
                    tempRange=tempVisualRange(j); %decrease visual range
                    %Recalculate angular size, contrast threshold and apparent
                    %contrast based on the decreased range
                    angularSize(i)=atan(T/(2*tempRange))*2*10^3;
                    Cr(i)=C0*integral(CrFunc,lambda1,lambda2);
                    Kt(i)=liminalContrast(D,BAerial.(conditions{c}),angularSize(i));
                    actRangeAerial(i)=tempRange;
                    j=j+1;
                    clc;
                    fprintf('condition: %s\npupil diameter: %f\nrange: %f\nerror: %f\n',...
                        conditions{c}, D, tempRange, (10^(Kt(i)))*1.15-abs(Cr(i)));
                end
                tempRangePrev=tempRange;
            end
        end
    end

    visualRangeAerial=actRangeAerial;
save([BIGEYEROOT 'fig04_visualrange/aerial_model/parameter_sensitivity/visibilityAerialContrastSensitivity.mat'], 'visualRangeAerial','DrangeAerial');
   
function Kt = liminalContrast(A,L,angularSize)

    %% Function definitions, section: contrast threshold as a function of pupil diameter
    %delta0 is labeled deltar
    %initial definitions -------
    logLa=@(L) log10(L);
    iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();
    deltat=@(L) iif( -2<logLa(L) && logLa(L)<=4, @() (1.22*((555e-9)/A)*10^3),...
    logLa(L)<=-2, @() (1.22*((507e-9)/A))*10^3);%Supplementary Appendix, pg. 5
    %----------
    %Equations adapted from Aksyutov
    Psi=@(L) (logLa(L)+3.5)/2.75; %Supplementary Appendix, pg.5 Eq. 10 part 1
    logdeltarStar= @(L) iif(logLa(L)>1.535, @() 0.928,...
        -0.75<=logLa(L) && logLa(L)<=1.535, @() 1-(0.5*(Psi(L)^(-3.2))),...
        -2.089<=logLa(L) && logLa(L)<-0.75, @()-0.28+0.1*Psi(L)^(3.2),...
        logLa(L)<-2.089, @() 0.14*logLa(L)+0.442);%Supplementary Appendix, pg. 5 Eq. 10 part 2
    
    %Not provided in supplement as seen in Aksyutov 200
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
        Xmparam=Xm(L);
        if Xparam>0 && Xparam<=0.20
            Yval=1.45*Xparam; 
            Kt =-2*logStardeltar-Yval;
        elseif Xparam>0.20
            Yval=Y(deltarparam,L);
            Kt=-2*logStardeltar-Yval;
        end
    end        