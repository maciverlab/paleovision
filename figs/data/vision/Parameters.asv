function Parameters
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
    minpupil=0.001; % largest diameter of pupil, meters
    maxpupil=0.025; % smallest diameter of pupil, meters

    % SENSORY VOLUME PARAMS
    % aerial and water half elevation angle of sensory volume
    elevationCoastal=pi/6; %30 deg 

    % azimuth is here assuming 35 degree overlap typical of fish;
    % definitely underestimate for terrestrial case
    azimuthCoastal=(305)*(pi/180);
    azimuthAir = azimuthCoastal;

    elevationMin=pi/2-elevationCoastal;
    elevationMax=pi/2+elevationCoastal;
    elevationMinAir=pi/2;
    elevationMaxAir=pi/2+elevationCoastal;

    azimuthMin=0;
    azimuthMax=azimuthCoastal;
    azimuthMaxAir=azimuthAir;

    %% AEIRAL MODEL
    lambda=(355:10:745)';
    %photoeceptor detection efficiency
    qAerial.Daylight=0.5; qAerial.Moonlight=0.36; qAerial.Starlight=0.36; 
    %qAerialVals=[qAerial_Daylight,qAerial_Moonlight,qAerial_Starlight];

    %Luminance (All values from Middleton)
    Wlambdaylambda=csvread('Wlambda.csv');
    BAerial.Daylight=1e3; % daylight luminance in cd/m^2 
    BAerial.Moonlight=1e-2; %Fairly brigh moonlight in cd/m^2
    BAerial.Starlight=1e-4; %moonless clear night sky in cd/m^2
    %BAerial=[BAerial_Daylight,BAerial_Moonlight,BAerial_Starlight];

    %Integration Time
    DtAerial.Daylight=(BAerial.Daylight)^-0.19; %Donner etal 1994
    DtAerial.Moonlight=(BAerial.Moonlight)^-0.19;
    DtAerial.Starlight=(BAerial.Starlight)^-0.19;
    %DtAerial=[DtAerial_Daylight, DtAerial_Moonlight, DtAerial_Starlight];

    % Contrast parameters. Miller uses +/- 0.5, +/-1, and +/-2 as span 
    C0Aerial.Daylight=-1; %aerial daylight contrast value, target black
    C0Aerial.Moonlight=-1; %aerial moonlight contrast value, target black 
    C0Aerial.Starlight=-1; %aerial starlight contrast value, taget black
    %C0Aerial=[C0Aerial_Daylight, C0Aerial_Moonlight, C0Aerial_Starlight];

    %Dark Noise
    XAerial=0.08;% units: photons/s, dark noise rate Rh/photoreceptor @23.5C

    %F-number
    FAerial.Daylight=8.8; %F-number: focal length/pupil diamter(D) for bright light
    FAerial.Moonlight=2.1; FAerial.Starlight=2.1; %F-number: focal lenght/D for light starved
    %FAerial=[FAerial_Daylight,FAerial_Moonlight,FAerial_Starlight];
    
    %Intensitiy Parameter Coefficient
    BbarAerial=(1.31e15)/0.89;

    %% AQUATIC MODEL
    qAquatic=0.36; %Efficiency
    DtAquatic=1.16; %Integration time
    Dt=1.16;

    %Contrast parameters
    C0Aquatic.Daylight=-1; %aquatic contrast value, black target
    C0Aquatic.Moonlight=-1; C0Aquatic.Starlight=-1;

    %Dark-noise
    XAquatic=0.011; %Rh/photoreceptor dark-noise rate @16.5C

    %Mattheisien's Ratio
    M=2.55; %units: n/a, ratio of focal length and pupil radius (2f/A), set to Matthiessen's ratio

    ybarAquaticLambda=[0.0000 0.0001 0.0004 0.0012 0.0040 0.0116 0.02300 .0380 0.0600 0.0910,...
        0.1390 0.2080 0.3230 0.5030 0.7100 0.8620 0.9540 0.9950 0.9950 0.9520,...
        0.8700 0.7570 0.6310 0.5030 0.3810 0.2650 0.1750 0.1070 0.0610 0.0320 0.0170 0.0082,...
        0.0041 0.0021 0.0011 0.0005 0.0003 0.0001 0.0001 0.0000]; %photopic luminosity function, Mobley Light and Water book
    lambdabar=380:10:770; %luminosity function wavelength domain
    ybarAquaticInterp=@(l) interp1(lambdabar,ybarAquaticLambda,l,'pchip');

    % DAYLIGHT
    %Absorption Coefficient
    aAquatic.Daylight=xlsread([parametersPath,'/hydrolight/base_sun/Hydrolight_BrownWater.xlsx'],'a_Model');
    aAquatic.Daylight=aAquatic.Daylight(:,1);
    %Scattering Coefficient
    bAquatic.Daylight=xlsread([parametersPath,'/hydrolight/base_sun/Hydrolight_BrownWater.xlsx'],'b_Model');
    bAquatic.Daylight=bAquatic.Daylight(:,1);
    %Diffuse Spectral Attenuation Coeff
    KuAquatic.Daylight=xlsread([parametersPath,'/hydrolight/base_sun/Hydrolight_BrownWater.xlsx'],'Ku');
    KuAquatic.Daylight=KuAquatic.Daylight(:,:);

    KdAquatic.Daylight=xlsread([parametersPath,'/hydrolight/base_sun/Hydrolight_BrownWater.xlsx'],'Kd');
    KdAquatic.Daylight=KdAquatic.Daylight(:,:);

    KhAquatic.Daylight=zeros(size(KdAquatic.Daylight,1),size(KdAquatic.Daylight,2));
    %Spectral Radiance
    %Upwelling, horizontal, downwelling radiance
    LuAquatic.Daylight=xlsread([parametersPath,'/hydrolight/base_sun/Hydrolight_BrownWater.xlsx'],'Lu');
    LuAquatic.Daylight=LuAquatic.Daylight(:,:)*5.03e15;

    LhAquatic.Daylight=xlsread([parametersPath,'/hydrolight/base_sun/Hydrolight_BrownWater.xlsx'],'Lh_2');
    LhAquatic.Daylight=LhAquatic.Daylight(:,:)*5.03e15;

    LdAquatic.Daylight=xlsread([parametersPath,'/hydrolight/base_sun/Hydrolight_BrownWater.xlsx'],'Ld');
    LdAquatic.Daylight=LdAquatic.Daylight(:,:)*5.03e15;
    %Luminance
    ybarAquatic.Daylight=ybarAquaticInterp(lambda);

    BuAquatic.Daylight=zeros(size(LuAquatic.Daylight,2),1);
    BhAquatic.Daylight=zeros(size(LhAquatic.Daylight,2),1);
    BdAquatic.Daylight=zeros(size(LdAquatic.Daylight,2),1);
    for i=1:size(LuAquatic.Daylight,2)
        tempU=(LuAquatic.Daylight(:,i)/5.03e15).*ybarAquatic.Daylight.*lambda;
        BuAquatic.Daylight(i)=trapz(lambda,tempU); 

        tempH=(LhAquatic.Daylight(:,i)/5.03e15).*ybarAquatic.Daylight.*lambda;
        BhAquatic.Daylight(i)=trapz(tempH);

        tempD=(LdAquatic.Daylight(:,i)/5.03e15).*ybarAquatic.Daylight.*lambda;
        BdAquatic.Daylight(i)=trapz(tempD);
    end

    % MOONLIGHT
    %Absorption Coefficient
    aAquatic.Moonlight=xlsread([parametersPath,'/hydrolight/base_moon/Mbase_moon.xls'],'a');
    aAquatic.Moonlight=aAquatic.Moonlight(:,1);
    %Scattering Coefficient
    bAquatic.Moonlight=xlsread([parametersPath,'/hydrolight/base_moon/Mbase_moon.xls'],'b');
    bAquatic.Moonlight=bAquatic.Moonlight(:,1);
    %Diffuse Spectral Attenuation Coeff
    KuAquatic.Moonlight=xlsread([parametersPath,'/hydrolight/base_moon/Mbase_moon.xls'],'Ku');
    KuAquatic.Moonlight=KuAquatic.Moonlight(:,:);

    KdAquatic.Moonlight=xlsread([parametersPath,'/hydrolight/base_moon/Mbase_moon.xls'],'Kd');
    KdAquatic.Moonlight=KdAquatic.Moonlight(:,:);

    KhAquatic.Moonlight=zeros(size(KdAquatic.Moonlight,1),size(KdAquatic.Moonlight,2));
    %Spectral Radiance
    %Upwelling, horizontal,downwelling radiance
    LuAquatic.Moonlight=xlsread([parametersPath,'/hydrolight/base_moon/Mbase_moon.xlsx'],'Lu');
    LuAquatic.Moonlight=LuAquatic.Moonlight(:,:)*5.03e15;

    LdAquatic.Moonlight=xlsread([parametersPath,'/hydrolight/base_moon/Mbase_moon.xlsx'],'Ld');
    LdAquatic.Moonlight=LdAquatic.Moonlight(:,:)*5.03e15;

    LhAquatic.Moonlight=xlsread([parametersPath,'/hydrolight/base_moon/Mbase_moon.xlsx'],'Lh_2');
    LhAquatic.Moonlight=LhAquatic.Moonlight(:,:)*5.03e15;
    %Luminance
    ybarAquatic.Moonlight=ybarAquaticInterp(lambda);

    BuAquatic.Moonlight=zeros(size(LuAquatic.Moonlight,2),1);
    BhAquatic.Moonlight=zeros(size(LhAquatic.Moonlight,2),1);
    BdAquatic.Moonlight=zeros(size(LdAquatic.Moonlight,2),1);
    for i=1:size(LuAquatic.Moonlight,2)
        tempU=(LuAquatic.Moonlight(:,i)/5.03e15).*ybarAquatic.Moonlight.*lambda;
        BuAquatic.Moonlight(i)=trapz(lambda,tempU);

        tempH=(LhAquatic.Moonlight(:,i)/5.03e15).*ybarAquatic.Moonlight.*lambda;
        BhAquatic.Moonlight(i)=trapz(lambda,tempH);

        tempD=(LdAquatic.Moonlight(:,i)/5.03e15).*ybarAquatic.Moonlight.*lambda;
        BdAquatic.Moonlight(i)=trapz(lambda,tempD);

    end

    % STARLIGHT
    %Absorption Coefficient
    aAquatic.Starlight=xlsread([parametersPath,'/hydrolight/base_stars/Mbase_stars.xls'],'a');
    aAquatic.Starlight=aAquatic.Starlight(:,1);
    %Scattering Coefficient
    bAquatic.Starlight=xlsread([parametersPath,'/hydrolight/base_stars/Mbase_stars.xls'],'b');
    bAquatic.Starlight=bAquatic.Starlight(:,1);
    %Diffuse Spectral Attenuation Coeff
    KuAquatic.Starlight=xlsread([parametersPath,'/hydrolight/base_stars/Mbase_stars.xls'],'Ku');
    KuAquatic.Starlight=KuAquatic.Starlight(:,:);

    KdAquatic.Starlight=xlsread([parametersPath,'/hydrolight/base_stars/Mbase_stars.xls'],'Kd');
    KdAquatic.Starlight=KdAquatic.Starlight(:,:);

    KhAquatic.Starlight=zeros(size(KdAquatic.Starlight,1),size(KdAquatic.Starlight,2));
    %Spectral Radiance
    %Upwelling,horizontal,downwellling radiance
    LuAquatic.Starlight=xlsread([parametersPath,'/hydrolight/base_stars/Mbase_stars.xls'],'Lu');
    LuAquatic.Starlight=LuAquatic.Starlight(:,:)*5.03e15;

    LdAquatic.Starlight=xlsread([parametersPath,'/hydrolight/base_stars/Mbase_stars.xls'],'Ld');
    LdAquatic.Starlight=LdAquatic.Starlight(:,:)*5.03e15;

    LhAquatic.Starlight=xlsread([parametersPath,'/hydrolight/base_stars/Mbase_stars.xls'],'Lh_2');
    LhAquatic.Starlight=LhAquatic.Starlight(:,:)*5.03e15;
    %Luminance
    ybarAquatic.Starlight=ybarAquaticInterp(lambda);

    BuAquatic.Starlight=zeros(size(LuAquatic.Starlight,2),1);
    BhAquatic.Starlight=zeros(size(LhAquatic.Starlight,2),1);
    BdAquatic.Starlight=zeros(size(LdAquatic.Starlight,2),1);
    for i=1:size(LuAquatic.Starlight,2)
        tempU=(LuAquatic.Starlight(:,i)/5.03e15).*ybarAquatic.Starlight.*lambda;
        BuAquatic.Starlight(i)=trapz(lambda,tempU);

        tempH=(LhAquatic.Starlight(:,i)/5.03e15).*ybarAquatic.Starlight.*lambda;
        BhAquatic.Starlight(i)=trapz(lambda,tempH);

        tempD=(LdAquatic.Starlight(:,i)/5.03e15).*ybarAquatic.Starlight.*lambda;
        BdAquatic.Starlight(i)=trapz(lambda,tempD);
    end

    % PHOTORECEPTOR ABSORPTION
    A=1; a0A=800; a1A=3.1;
    B=0.5; a0B=176; a1B=1.52;

    [~, ind_Su]=max(LdAquatic.Daylight);
    lambdaMax_Su=lambda(ind_Su(1));
    [~,ind_M]=max(LdAquatic.Moonlight);
    lambdaMax_M=lambda(ind_M(1));
    [~,ind_St]=max(LdAquatic.Starlight);
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
    save([parametersPath,'/Parameters.mat'])