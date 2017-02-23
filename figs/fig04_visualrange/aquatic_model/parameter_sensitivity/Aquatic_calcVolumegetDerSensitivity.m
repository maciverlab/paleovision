function [visualRangeAquatic, visualVolumeAquatic, drdAAquatic, dVdAAquatic,DrangeAquatic] = Aquatic_calcVolumegetDerSensitivity(CONTRASTTHRESH)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Volume and derivative calculations for visual range obtained from parameter perturbations
%%
%% Title                : A massive increase in visual range preceded the origin of terrestrial vertebrates
%% Authors              : Ugurcan Mugan, Malcolm A. MacIver
%% Authors' Affiliation : Northwestern University
%% This work is licensed under the Creative Commons Attribution 4.0 International License. 
%% To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
%% January 2017
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Intialize Variables
    run Parameters.m
    load('Parameters.mat')
    if CONTRASTTHRESH
        load('visibilityAquaticSensitivityParameters.mat')
        temp=visualRangeAquatic;
        load('visibilityAquaticSensitivityContrast.mat');
        visualRangeAquatic=reshape(visualRangeAquatic,[30,2,1]);
    else
        load('meteoAquaticSensitivityParameters.mat')
    end
    temp(:,:,size(temp,3)+1)=visualRangeAquatic; visualRangeAquatic=temp;               
    drdAAquatic=visualRangeAquatic;
    dVdAAquatic=visualRangeAquatic;
    visualVolumeAquatic=visualRangeAquatic;

%% Calculate sector volume
    for i=1:size(visualRangeAquatic,3);
        for j=1:size(visualRangeAquatic,2);
            for k=1:size(visualRangeAquatic,1);
                visualVolumeAquatic(k,j,i)=integral3(f,0,visualRangeAquatic(k,j,i),...
                elevationMin,elevationMax,...
                azimuthMin,azimuthMax);
            end
        end
    end

%% Calculate derivatives
    for i=1:size(visualRangeAquatic,3);
        for j=1:size(visualRangeAquatic,2);
            drdAAquatic(:,j,i)=derivative(DrangeAquatic*10^3,visualRangeAquatic(:,j,i));
            dVdAAquatic(:,j,i)=derivative(DrangeAquatic*10^3,visualVolumeAquatic(:,j,i));
        end
    end
    
function [ dydx ] = derivative( x,y )
    dydx(1)=(y(2)-y(1))/(x(2)-x(1));
    for n=2:length(y)-1
        dydx(n)=(y(n+1)-y(n-1))/(x(n+1)-x(n-1));
    end
    dydx(length(y))=(y(end)-y(end-1))/(x(end)-x(end-1));