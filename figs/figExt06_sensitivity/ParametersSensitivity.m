function ParametersSensitivity
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter initializations from Hydrolight data
%%
%% Title                : A massive increase in visual range preceded the origin of terrestrial vertebrates
%% Authors              : Ugurcan Mugan, Malcolm A. MacIver
%% Authors' Affiliation : Northwestern University
%% This work is licensed under the Creative Commons Attribution 4.0 International License. 
%% To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
%% January 2017
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%global BIGEYEROOT
parametersPath=[BIGEYEROOT 'data/vision/'];
close all;

    q=0.36; %units: n/a, detection efficien
    eta=0.5; %Walzl et al 2007
    Dt=1.16; %units: s, integration time, used typical value
    k=0.035; % photoreceptor absorbtion, units 1/micrometers
    len=57;  % length of photoreceptor, in micrometers
    X=0.011; %units: photons/s, dark-noise rate/photoreceptor @16.5degrees Celsius
    R=1.96; %units: n/a, reliability coefficient for 95% confidence, used typical value
    d=3e-6; %units: m, photoreceptor diameter, used typical value
    M=2.55; %units: n/a, ratio of focal length and pupil radius (2f/A), set to Matthiessen's ratio
    T=.1; % units: m, prey width

    minpupil=0.001; % largest diameter of pupil, meters
    maxpupil=0.03; % smallest diameter of pupil, meters
    C0=-1;
    ybarAquaticLambda=[0.0000 0.0001 0.0004 0.0012 0.0040 0.0116 0.02300 .0380 0.0600 0.0910,...
        0.1390 0.2080 0.3230 0.5030 0.7100 0.8620 0.9540 0.9950 0.9950 0.9520,...
        0.8700 0.7570 0.6310 0.5030 0.3810 0.2650 0.1750 0.1070 0.0610 0.0320 0.0170 0.0082,...
        0.0041 0.0021 0.0011 0.0005 0.0003 0.0001 0.0001 0.0000]; %photopic luminosity function, Mobley Light and Water book
    lambdabar=380:10:770; %luminosity function wavelength domain
    ybarAquaticInterp=@(l) interp1(lambdabar,ybarAquaticLambda,l,'pchip');
    %% Sensory sector volume parameters

    elevationCoastal=pi/3; %60 deg 

    elevationMin=pi/2-elevationCoastal;
    elevationMax=pi/2+elevationCoastal;
    azimuthMin=0;
    azimuthCoastal=(2*120-35)*(pi/180); %viewing azimuth
    azimuthMax=azimuthCoastal;

    f = @(rho,phi,theta) rho.^2.*sin(phi); %volume equation in spherical coordinates

    %% Beam attenuation coefficients, absorbtion
    %HighTurbidity
    a.HighTurbidity=xlsread('hydrolight/HighTurbidity/MHighTurbidity.xls','a');
    a.HighTurbidity=a.HighTurbidity(:,1);
    %CLEAR
    a.Clear=xlsread('hydrolight/Clear/MClear.xls','a');
    a.Clear=a.Clear(:,1);
    %HIGH ABSORPTION
    a.AbsDom=xlsread('hydrolight/AbsDom/MAbsDom.xls','a');
    a.AbsDom=a.AbsDom(:,1);
    %HIGH SCATTERING
    a.ScatDom=xlsread('hydrolight/ScatDom/MScatDom.xls','a');
    a.ScatDom=a.ScatDom(:,1);
    %% Scattering 
    %HIGH TURBIDITY
    b.HighTurbidity=xlsread('hydrolight/HighTurbidity/MHighTurbidity.xls','b');
    b.HighTurbidity=b.HighTurbidity(:,1);
    %CLEAR
    b.Clear=xlsread('hydrolight/Clear/MClear.xls','b');
    b.Clear=b.Clear(:,1);
    %HIGH ABSORPTION
    b.AbsDom=xlsread('hydrolight/AbsDom/MAbsDom.xls','b');
    b.AbsDom=b.AbsDom(:,1);
    %HIGH SCATTERING
    b.ScatDom=xlsread('hydrolight/ScatDom/MScatDom.xls','b');
    b.ScatDom=b.ScatDom(:,1);
    %% Diffuse attenuation coeff, K-func
    %HIGH TURBIDITY
    K.HighTurbidity.Upward=xlsread('hydrolight/HighTurbidity/MHighTurbidity.xls','Kd');
    K.HighTurbidity.Upward=K.HighTurbidity.Upward(:,:);

    K.HighTurbidity.Horizontal=zeros(size(K.HighTurbidity.Upward,1),size(K.HighTurbidity.Upward,2));
    %CLEAR
    K.Clear.Upward=xlsread('hydrolight/Clear/MClear.xls','Kd');
    K.Clear.Upward=K.Clear.Upward(:,:);

    K.Clear.Horizontal=zeros(size(K.Clear.Upward,1),size(K.Clear.Upward,2));
    %HIGH ABSORPTION
    K.AbsDom.Upward=xlsread('hydrolight/AbsDom/MAbsDom.xls','KLu');
    K.AbsDom.Upward=K.AbsDom.Upward(:,:);

    K.AbsDom.Horizontal=zeros(size(K.AbsDom.Upward,1),size(K.AbsDom.Upward,2));
    %HIGH SCATTERING
    K.ScatDom.Upward=xlsread('hydrolight/ScatDom/MScatDom.xls','Kd');
    K.ScatDom.Upward=K.ScatDom.Upward(:,:);

    K.ScatDom.Horizontal=zeros(size(K.ScatDom.Upward,1),size(K.ScatDom.Upward,2));

    %% Spectral radiance values
    lambda=(355:10:745)';
    %High TURBIDITY
    L.HighTurbidity.Upward=xlsread('hydrolight/HighTurbidity/MHighTurbidity.xls','Ld');
    L.HighTurbidity.Upward=L.HighTurbidity.Upward(:,:)*5.03e15;

    L.HighTurbidity.Horizontal=xlsread('hydrolight/HighTurbidity/MHighTurbidity.xls','Lh_2');
    L.HighTurbidity.Horizontal=L.HighTurbidity.Horizontal(:,:)*5.03e15;
    %CLEAR
    L.Clear.Upward=xlsread('hydrolight/Clear/MClear.xls','Ld');
    L.Clear.Upward=L.Clear.Upward(:,:)*5.03e15;

    L.Clear.Horizontal=xlsread('hydrolight/Clear/MClear.xls','Lh_2');
    L.Clear.Horizontal=L.Clear.Horizontal(:,:)*5.03e15;
    %HIGH ABSORPTION
    L.AbsDom.Upward=xlsread('hydrolight/AbsDom/MAbsDom.xls','Ld');
    L.AbsDom.Upward=L.AbsDom.Upward(:,:)*5.03e15;

    L.AbsDom.Horizontal=xlsread('hydrolight/AbsDom/MAbsDom.xls','Lh_2');
    L.AbsDom.Horizontal=L.AbsDom.Horizontal(:,:)*5.03e15;
    %HIGH SCATTERING
    L.ScatDom.Upward=xlsread('hydrolight/ScatDom/MScatDom.xls','Ld');
    L.ScatDom.Upward=L.ScatDom.Upward(:,:)*5.03e15;

    L.ScatDom.Horizontal=xlsread('hydrolight/ScatDom/MScatDom.xls','Lh_2');
    L.ScatDom.Horizontal=L.ScatDom.Horizontal(:,:)*5.03e15;

    %% Luminance calculations
    cond={'HighTurbidity','Clear','AbsDom','ScatDom'};
    ybar=ybarAquaticInterp(lambda);
    for j=1:length(cond)
        for i=1:size(L.ScatDom.Upward,2)
            tempH=(L.(cond{j}).Horizontal(:,i)/5.03e15).*ybar.*lambda;
            BLum.(cond{j}).Horizontal(i)=trapz(tempH);

            tempD=(L.(cond{j}).Upward(:,i)/5.03e15).*ybar.*lambda;
            BLum.(cond{j}).Upward(i)=trapz(tempD);
        end
    end
    %% Photoreceptor absorbtion based on Warrant et al
    A=1; a0A=800; a1A=3.1;
    B=0.5; a0B=176; a1B=1.52;

    [~,ind_B]=max(L.HighTurbidity.Upward);
    lambdaMax_B=lambda(ind_B(1));
    [~,ind_C]=max(L.Clear.Upward);
    lambdaMax_C=lambda(ind_C(1));
    [~,ind_HA]=max(L.AbsDom.Upward);
    lambdaMax_HA=lambda(ind_HA(1));
    [~,ind_HS]=max(L.ScatDom.Upward);
    lambdaMax_HS=lambda(ind_HS(1));

    pAbsorb.HighTurbidity=A*exp(-a0A*(log10(lambda./lambdaMax_B)).^2.*...
        (1+a1A*log10(lambda./lambdaMax_B)+(3*a1A^2/8).*log10(lambda./lambdaMax_B).^2)+...
        B*exp(-a0B*(log10(lambda./368)).^2.*...
        (1+a1B*log10(lambda./368)+(3*a1B^2/8)*log10(lambda./368))));

    pAbsorb.Clear=A*exp(-a0A*(log10(lambda./lambdaMax_C)).^2.*...
        (1+a1A*log10(lambda./lambdaMax_C)+(3*a1A^2/8).*log10(lambda./lambdaMax_C).^2)+...
        B*exp(-a0B*(log10(lambda./368)).^2.*...
        (1+a1B*log10(lambda./368)+(3*a1B^2/8)*log10(lambda./368))));

    pAbsorb.AbsDom=A*exp(-a0A*(log10(lambda./lambdaMax_HA)).^2.*...
        (1+a1A*log10(lambda./lambdaMax_HA)+(3*a1A^2/8).*log10(lambda./lambdaMax_HA).^2)+...
        B*exp(-a0B*(log10(lambda./368)).^2.*...
        (1+a1B*log10(lambda./368)+(3*a1B^2/8)*log10(lambda./368))));

    pAbsorb.ScatDom=A*exp(-a0A*(log10(lambda./lambdaMax_HS)).^2.*...
        (1+a1A*log10(lambda./lambdaMax_HS)+(3*a1A^2/8).*log10(lambda./lambdaMax_HS).^2)+...
        B*exp(-a0B*(log10(lambda./368)).^2.*...
        (1+a1B*log10(lambda./368)+(3*a1B^2/8)*log10(lambda./368))));
    
save([BIGEYEROOT 'figEXT06_sensitivity/ParametersSensitivity.mat'])