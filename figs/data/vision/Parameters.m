function Parameters
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Statistical summary of raw socket length data based on grouping
%%
%% Title                : A massive increase in visual range preceded the origin of terrestrial vertebrates
%% Authors              : Ugurcan Mugan, Malcolm A. MacIver
%% Authors' Affiliation : Northwestern University
%% This work is licensed under the Creative Commons Attribution 4.0 International License. 
%% To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
%% January 2017
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% function for initializing the parameters for the visual model
% a comprehensive list with values and explanations can be found in
% Supplementary Appendix Table 1
    close all;
    global BIGEYEROOT
    parametersPath=[BIGEYEROOT 'data/vision'];
    %% Common Parameters:
    k=0.035; % photoreceptor absorbtion, units 1/micrometers
    len=57;  % length of photoreceptor, in micrometers
    T=.1; % units: m, prey width
    f = @(rho,phi,theta) rho.^2.*sin(phi); %volume equation in spherical coordinates
    R=1.96; %units: n/a, reliability coefficient for 95% confidence, used typical value
    d=3e-6; %units: m, photoreceptor diameter, used typical value
    minpupil=1e-3; % largest diameter of pupil, meters
    maxpupil=25e-3; % smallest diameter of pupil, meters

    % SENSORY VOLUME PARAMS
    % aerial and water half elevation angle of sensory volume
    elevationCoastal=pi/6; %30 deg 

    % azimuth is here assuming 35 degree overlap typical of fish;
    % definitely underestimate for terrestrial case
    azimuthCoastal=(305)*(pi/180);
    % azimuth for crocodile vision (2*monoocular - binocular vision)
    azimuthAir = ((2*156)-25)*(pi/180);
    
    % sphere volume integral bounds for aerial and aquatic vision
    elevationMin=pi/2-elevationCoastal;
    elevationMax=pi/2+elevationCoastal;
    elevationMinAir=pi/2;
    elevationMaxAir=pi/2+elevationCoastal;

    azimuthMin=0;
    azimuthMax=azimuthCoastal;
    azimuthMaxAir=azimuthAir;

    %% Aerial model Parameters:
    % wavelength span (symbol: ?) (units: ?m)
    lambda=(355:10:745)';
    % Photoeceptor detection efficiency (symbol: ?)
    qAerial.Daylight=0.5; 
    qAerial.Moonlight=0.36; 
    qAerial.Starlight=0.36; 

    % Luminance and (All values from Middleton) (symbol: L, unit: cd/m^2)
    Wlambdaylambda=csvread('Wlambda.csv'); %spectral radiance data from Middleton
    BAerial.Daylight=1e3; % Daylight luminance 
    BAerial.Moonlight=1e-2; % Fairly bright clear sky moonlight
    BAerial.Starlight=1e-4; % Moonless clear night sky

    % Integration Time calculated based on Donner et al (symbol: ?t, unit: s)
    DtAerial.Daylight=(BAerial.Daylight)^-0.19; 
    DtAerial.Moonlight=(BAerial.Moonlight)^-0.19;
    DtAerial.Starlight=(BAerial.Starlight)^-0.19;

    % Object contrast (black: -1, snow: +4, bounds: [-1,inf)) (symbol: C_0)
    C0Aerial.Daylight=-1; 
    C0Aerial.Moonlight=-1; 
    C0Aerial.Starlight=-1;

    % Dark Noise (symbol: ?) (units: Rh/s)
    XAerial=0.08; % T=@23.5C

    % F-number (symbol: F)
    FAerial.Daylight=8.8; 
    FAerial.Moonlight=2.1; 
    FAerial.Starlight=2.1;

    % Intensitiy Parameter Coefficient (symbol: \bar{?})
    BbarAerial=(1.31e15)/0.89;

    %% Aquatic modelParameters:
    % Photoeceptor detection efficiency (symbol: ?)
    qAquatic=0.36;
    % Integration Time from Nilsson et al (symbol: ?t, unit: s)
    DtAquatic=1.16; 

     % Object contrast (black: -1, snow: +4, bounds: [-1,inf)) (symbol: C_0)
    C0Aquatic.Daylight=-1;
    C0Aquatic.Moonlight=-1; 
    C0Aquatic.Starlight=-1;

    % Dark Noise (symbol: ?) (units: Rh/s)
    XAquatic=0.011; %T=@16.5C

    %Mattheisien's Ratio (symbol: M) -> rearranged to give F=M/2 (F-number))
    M=2.55;

    %Aquatic Luminance Calculations for contrast threshold:
    lambdabar=380:10:770; %luminosity function wavelength domain
    % Photopic luminosity function (value at each sampled wavelength) from Mobley
    ybarAquaticLambda=[0.0000 0.0001 0.0004 0.0012 0.0040 0.0116 0.02300 .0380 0.0600 0.0910,...
        0.1390 0.2080 0.3230 0.5030 0.7100 0.8620 0.9540 0.9950 0.9950 0.9520,...
        0.8700 0.7570 0.6310 0.5030 0.3810 0.2650 0.1750 0.1070 0.0610 0.0320 0.0170 0.0082,...
        0.0041 0.0021 0.0011 0.0005 0.0003 0.0001 0.0001 0.0000];
    ybarAquaticInterp=@(l) interp1(lambdabar,ybarAquaticLambda,l,'pchip');

    % Aquatic daylight hydrolight data input (absorption, scattering,
    % spectral diffuse attenuation, radiance)
    
    % Absorption Coefficient
    aAquatic.Daylight=xlsread([parametersPath,'/hydrolight/base_sun/Hydrolight_BrownWater.xlsx'],'a_Model');
    aAquatic.Daylight=aAquatic.Daylight(:,1);
    
    % Scattering Coefficient
    bAquatic.Daylight=xlsread([parametersPath,'/hydrolight/base_sun/Hydrolight_BrownWater.xlsx'],'b_Model');
    bAquatic.Daylight=bAquatic.Daylight(:,1);
    
    % Diffuse Spectral Attenuation Coeff
    KAquatic.Daylight.Upward=xlsread([parametersPath,'/hydrolight/base_sun/Hydrolight_BrownWater.xlsx'],'Kd');
    KAquatic.Daylight.Upward=KAquatic.Daylight.Upward(:,:);

    KAquatic.Daylight.Horizontal=zeros(size(KAquatic.Daylight.Upward,1),size(KAquatic.Daylight.Upward,2));
    
    % Spectral Radiance
    % Upwelling, horizontal radiance
    LAquatic.Daylight.Upward=xlsread([parametersPath,'/hydrolight/base_sun/Hydrolight_BrownWater.xlsx'],'Ld');
    LAquatic.Daylight.Upward=LAquatic.Daylight.Upward(:,:)*5.03e15;

    LAquatic.Daylight.Horizontal=xlsread([parametersPath,'/hydrolight/base_sun/Hydrolight_BrownWater.xlsx'],'Lh_2');
    LAquatic.Daylight.Horizontal=LAquatic.Daylight.Horizontal(:,:)*5.03e15;
    
    % Luminance Calculations from ybar (line 91-98)
    ybarAquatic.Daylight=ybarAquaticInterp(lambda);
   
    BAquatic.Daylight.Upward=zeros(size(LAquatic.Daylight.Upward,2),1);
    BAquatic.Daylight.Horizontal=zeros(size(LAquatic.Daylight.Horizontal,2),1);
    for i=1:size(LAquatic.Daylight.Upward,2)
        tempD=(LAquatic.Daylight.Upward(:,i)/5.03e15).*ybarAquatic.Daylight.*lambda;
        BAquatic.Daylight.Upward(i)=trapz(tempD); 

        tempH=(LAquatic.Daylight.Horizontal(:,i)/5.03e15).*ybarAquatic.Daylight.*lambda;
        BAquatic.Daylight.Horizontal(i)=trapz(tempH);
    end

    % Aquatic moonlight hydrolight data input (absorption, scattering,
    % spectral diffuse attenuation, radiance)
    
    % Absorption Coefficient
    aAquatic.Moonlight=xlsread([parametersPath,'/hydrolight/base_moon/Mbase_moon.xls'],'a');
    aAquatic.Moonlight=aAquatic.Moonlight(:,1);
    
    % Scattering Coefficient
    bAquatic.Moonlight=xlsread([parametersPath,'/hydrolight/base_moon/Mbase_moon.xls'],'b');
    bAquatic.Moonlight=bAquatic.Moonlight(:,1);
    
    % Diffuse Spectral Attenuation Coeff
    KAquatic.Moonlight.Upward=xlsread([parametersPath,'/hydrolight/base_moon/Mbase_moon.xls'],'Kd');
    KAquatic.Moonlight.Upward=KAquatic.Moonlight.Upward(:,:);

    KAquatic.Moonlight.Horizontal=zeros(size(KAquatic.Moonlight.Upward,1),size(KAquatic.Moonlight.Upward,2));
    
    %Spectral Radiance
    %Upwelling, horizontal,downwelling radiance
    LAquatic.Moonlight.Upward=xlsread([parametersPath,'/hydrolight/base_moon/Mbase_moon.xlsx'],'Ld');
    LAquatic.Moonlight.Upward=LAquatic.Moonlight.Upward(:,:)*5.03e15;

    LAquatic.Moonlight.Horizontal=xlsread([parametersPath,'/hydrolight/base_moon/Mbase_moon.xlsx'],'Lh_2');
    LAquatic.Moonlight.Horizontal=LAquatic.Moonlight.Horizontal(:,:)*5.03e15;
    
    % Luminance Calculations from ybar (line 91-98)
    ybarAquatic.Moonlight=ybarAquaticInterp(lambda);

    BAquatic.Moonlight.Upward=zeros(size(LAquatic.Moonlight.Upward,2),1);
    BAquatic.Moonlight.Horizontal=zeros(size(LAquatic.Moonlight.Horizontal,2),1);
    for i=1:size(LAquatic.Moonlight.Upward,2)
        tempD=(LAquatic.Moonlight.Upward(:,i)/5.03e15).*ybarAquatic.Moonlight.*lambda;
        BAquatic.Moonlight.Upward(i)=trapz(lambda,tempD);

        tempH=(LAquatic.Moonlight.Horizontal(:,i)/5.03e15).*ybarAquatic.Moonlight.*lambda;
        BAquatic.Moonlight.Horizontal(i)=trapz(lambda,tempH);
    end

    % Aquatic starlight hydrolight data input (absorption, scattering,
    % spectral diffuse attenuation, radiance)
    % Absorption Coefficient
    aAquatic.Starlight=xlsread([parametersPath,'/hydrolight/base_stars/Mbase_stars.xls'],'a');
    aAquatic.Starlight=aAquatic.Starlight(:,1);
    % Scattering Coefficient
    bAquatic.Starlight=xlsread([parametersPath,'/hydrolight/base_stars/Mbase_stars.xls'],'b');
    bAquatic.Starlight=bAquatic.Starlight(:,1);
    % Diffuse Spectral Attenuation Coeff
    KAquatic.Starlight.Upward=xlsread([parametersPath,'/hydrolight/base_stars/Mbase_stars.xls'],'Kd');
    KAquatic.Starlight.Upward=KAquatic.Starlight.Upward(:,:);

    KAquatic.Starlight.Horizontal=zeros(size(KAquatic.Starlight.Upward,1),size(KAquatic.Starlight.Upward,2));
    
    % Spectral Radiance
    % Upwelling,horizontal,downwellling radiance
    LAquatic.Starlight.Upward=xlsread([parametersPath,'/hydrolight/base_stars/Mbase_stars.xls'],'Ld');
    LAquatic.Starlight.Upward=LAquatic.Starlight.Upward(:,:)*5.03e15;

    LAquatic.Starlight.Horizontal=xlsread([parametersPath,'/hydrolight/base_stars/Mbase_stars.xls'],'Lh_2');
    LAquatic.Starlight.Horizontal=LAquatic.Starlight.Horizontal(:,:)*5.03e15;
    % Luminance Calculations from ybar (line 91-98)
    ybarAquatic.Starlight=ybarAquaticInterp(lambda);

    BAquatic.Starlight.Upward=zeros(size(LAquatic.Starlight.Upward,2),1);
    BAquatic.Starlight.Horizontal=zeros(size(LAquatic.Starlight.Horizontal,2),1);
    for i=1:size(LAquatic.Starlight.Upward,2)
        tempD=(LAquatic.Starlight.Upward(:,i)/5.03e15).*ybarAquatic.Starlight.*lambda;
        BAquatic.Starlight.Upward(i)=trapz(lambda,tempD);

        tempH=(LAquatic.Starlight.Horizontal(:,i)/5.03e15).*ybarAquatic.Starlight.*lambda;
        BAquatic.Starlight.Horizontal(i)=trapz(lambda,tempH);
    end

    % Calcuation of photoreceptor absoprtion curve based on the SSH
    % template from Warrant et al
    A=1; a0A=800; a1A=3.1;
    B=0.5; a0B=176; a1B=1.52;

    [~, ind_Su]=max(LAquatic.Daylight.Upward);
    lambdaMax_Su=lambda(ind_Su(1));
    [~,ind_M]=max(LAquatic.Moonlight.Upward);
    lambdaMax_M=lambda(ind_M(1));
    [~,ind_St]=max(LAquatic.Starlight.Upward);
    lambdaMax_St=lambda(ind_St(1));

    pAbsorbAquatic.Daylight=A*exp(-a0A*(log10(lambda./lambdaMax_Su)).^2.*...
        (1+a1A*log10(lambda./lambdaMax_Su)+(3*a1A^2/8).*log10(lambda./lambdaMax_Su).^2)+...
        B*exp(-a0B*(log10(lambda./368)).^2.*...
        (1+a1B*log10(lambda./368)+(3*a1B^2/8)*log10(lambda./368))));

    pAbsorbAquatic.Moonlight=A*exp(-a0A*(log10(lambda./lambdaMax_M)).^2.*...
        (1+a1A*log10(lambda./lambdaMax_M)+(3*a1A^2/8).*log10(lambda./lambdaMax_M).^2)+...
         B*exp(-a0B*(log10(lambda./368)).^2.*...
        (1+a1B*log10(lambda./368)+(3*a1B^2/8)*log10(lambda./368))));

    pAbsorbAquatic.Starlight=A*exp(-a0A*(log10(lambda./lambdaMax_St)).^2.*...
        (1+a1A*log10(lambda./lambdaMax_St)+(3*a1A^2/8).*log10(lambda./lambdaMax_St).^2)+...
         B*exp(-a0B*(log10(lambda./368)).^2.*...
        (1+a1B*log10(lambda./368)+(3*a1B^2/8)*log10(lambda./368))));
    
%% save parameters for later use    
    save([parametersPath,'/Parameters.mat'])