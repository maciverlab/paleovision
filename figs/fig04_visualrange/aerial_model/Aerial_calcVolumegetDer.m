function [visualRangeAerial,visualVolumeAerial,drdAAerial,dVdAAerial,DrangeAerial]=Aerial_calcVolumegetDer(CONTRASTTHRESH)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation of visual volume, derivatives of visual range and visual volume
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
       load('visibilityAerial.mat')
    else
        load('meteoAerial.mat');       
    end

    %% Calculate spherical sector for visual volume
    
    for loop1=1:size(visualRangeAerial,2)
        for loop2=1:size(visualRangeAerial,1)
            visualVolumeAerial(loop2,loop1)=integral3(f,0,(visualRangeAerial(loop2,loop1)),...
                elevationMinAir,elevationMaxAir,...
                azimuthMin,azimuthMaxAir);
        end
    end

    %% Calculate visual range and visual volume derivatives
    
    for loop1=1:size(visualRangeAerial,2);
        drdAAerial(:,loop1)=derivative(DrangeAerial*10^3,visualRangeAerial(:,loop1));
        dVdAAerial(:,loop1)=derivative(DrangeAerial*10^3,visualVolumeAerial(:,loop1));
    end
    
function [ dydx ] = derivative( x,y )
    dydx(1)=(y(2)-y(1))/(x(2)-x(1));
    for n=2:length(y)-1
        dydx(n)=(y(n+1)-y(n-1))/(x(n+1)-x(n-1));
    end
    dydx(length(y))=(y(end)-y(end-1))/(x(end)-x(end-1));
