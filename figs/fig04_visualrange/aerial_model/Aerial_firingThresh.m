function Aerial_firingThresh
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation of visual range based on optical properties
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
    conditions={'Daylight','Moonlight','Starlight'};
    DrangeAerial=linspace(minpupil,maxpupil,25);
    
    Wlambdaylambda=csvread('Wlambda.csv'); %Spectral radiance emittence, Table from Middleton
    WlambdaylambdaInterp= @(lambda) interp1(Wlambdaylambda(:,1),Wlambdaylambda(:,2),lambda,'pchip'); 
    %extinction coeff from Middleton (originally in 1/km coverted to 1/m)
    sigma=@(lambda) ((1.1e-3*lambda.^(-4))+(0.008*lambda.^(-2.09)))/(1e3);
    lambda1=0.4; lambda2=0.7; %visible range wavelength
    %% Calculate visual range given viewing condition and pupil diameter
    visualRangeAerial=zeros(length(DrangeAerial),length(conditions));
    for c=1:length(conditions)
        %following Middleton:
        Bh=BAerial.(conditions{c})*integral(WlambdaylambdaInterp,lambda1,lambda2); %horizon luminance
        %Supplementary Appendix Materials and Methods, pg.2

        %Intensitiy Parameter coefficient, supplement pg 5 line 550
        %Full equation to get Rh is on line 549 Rh=beta*B_h,o(lambda)
        %horizon space-light photons/m^2 s sr, value checked with mathematica
        Rh=((1.31e3)/0.89)*Bh*(1e6)^2; %Supplementary Appendix, pg.3 and Table S1
        
        if strcmp(conditions{c},'Daylight'); %Arrange gridding based on target diameter and viewing condition
            if T<0.05; 
                minvisualrange=10; maxvisualrange=2000;
            else
                minvisualrange=100; maxvisualrange=50000;
            end
        elseif strcmp(conditions{c},'Moonlight');
            if T<0.05
                minvisualrange=1; maxvisualrange=500;
            else
                minvisualrange=50; maxvisualrange=3300;
            end
        elseif strcmp(conditions{c},'Starlight');
            if T<0.05
                minvisualrange=1; maxvisualrange=200;
            else
                minvisualrange=1; maxvisualrange=750;
            end
        end
        
        visualRangeTemp=linspace(minvisualrange,maxvisualrange,1000);        
        for j=1:length(DrangeAerial)
            D=DrangeAerial(j);
            P=zeros(length(visualRangeTemp),1);
            for i=1:length(visualRangeTemp)
                r=visualRangeTemp(i);
                P(i)=firingThreshold(BAerial.(conditions{c}),Rh,D,T,r,sigma,WlambdaylambdaInterp,...
                    XAerial,qAerial.(conditions{c}),DtAerial.(conditions{c}),k,len,d,...
                    FAerial.(conditions{c}),C0Aerial.(conditions{c}),...
                    lambda1, lambda2,BbarAerial,R,conditions{c});           
            end
            IDX=knnsearch(P,1,'distance','euclidean'); %find closest range value to soln 1
            visualRangeAerial(j,c)=visualRangeTemp(IDX); %store range value
            clc;
            fprintf('viewing condition: %s\niteration: %d\npupil diameter: %f\nvisual range: %f\n',...
                conditions{c},j,D,visualRangeAerial(j,c));
        end
    end

save([BIGEYEROOT 'fig04_visualrange/aerial_model/meteoAerial.mat'],'visualRangeAerial','DrangeAerial')
    
function soln=firingThreshold(B,Rh,D,T,r,sigma,WlambdaylambdaInterp,X,q,Dt,k,len,d,F,C0,lambda1,lambda2,BbarAerial,R,cond)
    
    Nh=(pi/4)^2*(T/r)^2*D^2*Rh*Dt*q*((k*len)/(2.3+(k*len))); %Supplementary Appendix, combination of Eq. 5 and Rh, first part of Eq. 6
    if ~strcmp(cond,'Daylight')
         Nh=(pi/4)^2*(T/r)^2*D^2*Rh*Dt*q*(1-exp(-k*len));
    end

    Nfalse=((T*F*D)/(r*d))^2*X*Dt ; %Supplementary Appendix pg. 3

    Bofunc=@(lambda) WlambdaylambdaInterp(lambda).*...
        (1+(C0.*exp(-sigma(lambda).*r)));
    Bo= B*integral(Bofunc,lambda1,lambda2);
    Ro=((1.31e3)/0.89)*Bo*(1e6)^2; %Supplementary Appendix, Eq. 3 pg. 2

    No=(pi/4)^2*(T/r)^2*D^2*Ro*Dt*q*((k*len)/(2.3+(k*len)));  %Supplementary Appendix, combination of Eq. 5 and Rh, first part of Eq. 6
    if~strcmp(cond,'Daylight')
        No=(pi/4)^2*(T/r)^2*D^2*Ro*Dt*q*(1-exp(-k*len));
    end
    
    soln=(R*sqrt(No+Nh+2*Nfalse))/(abs(No-Nh)); %Supplementary Appendix, Eq. 4 (pg. 2) and Eq. 6 (pg. 3)