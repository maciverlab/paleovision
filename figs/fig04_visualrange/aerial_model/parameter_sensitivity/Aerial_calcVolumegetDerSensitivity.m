function [visualRangeAerial,visualVolumeAerial,drdAAerial,dVdAAerial,DrangeAerial]=Aerial_calcVolumegetDerSensitivity(CONTRASTTHRESH)

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

%% Initialize variables    
    run Parameters.m
    load('Parameters.mat')

    if CONTRASTTHRESH
        fileNames={'visibilityAerialParameterSensitivity.mat'};
        load(fileNames{1});
        visualRangeParameter=visualRangeAerial;
        load('visibilityAerialContrastSensitivity.mat');                  
    else
        fileNames={'meteoAerialParameterSensitivity.mat'};        
        load(fileNames{1});                
    end
    
    visualRangeAerial=[visualRangeAerial visualRangeParameter];

    %% CALCULATE VOLUME
    for k=1:size(visualRangeAerial,3)
        for j=1:size(visualRangeAerial,2)
            for i=1:size(visualRangeAerial,1)
                visualVolumeAerial(i,k,j)=integral3(f,0,(visualRangeAerial(i,j,k)),...
                    elevationMinAir,elevationMaxAir,...
                    azimuthMin,azimuthMaxAir);
            end
        end
    end

    %% Calculate visual range and visual volume derivatives
    for k=1:size(visualRangeAerial,3)
        for j=1:size(visualRangeAerial,2)
            drdAAerial(:,k,j)=derivative(DrangeAerial*10^3,visualRangeAerial(:,j,k));
            dVdAAerial(:,k,j)=derivative(DrangeAerial*10^3,visualVolumeAerial(:,k,j));
        end
    end

visualRangeAerial=reshape(visualRangeAerial,[size(visualRangeAerial,1),size(visualRangeAerial,3),size(visualRangeAerial,2)]);

function [ dydx ] = derivative( x,y )
    dydx(1)=(y(2)-y(1))/(x(2)-x(1));
    for n=2:length(y)-1
        dydx(n)=(y(n+1)-y(n-1))/(x(n+1)-x(n-1));
    end
    dydx(length(y))=(y(end)-y(end-1))/(x(end)-x(end-1));
