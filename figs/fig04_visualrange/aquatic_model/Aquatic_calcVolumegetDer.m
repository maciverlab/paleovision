function [visualRangeAquatic, visualVolumeAquatic, drdAAquatic, dVdAAquatic,DrangeAquatic] = Aquatic_calcVolumegetDer(CONTRASTTHRESH)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation of volume and derivatives for aquatic conditions
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
        load('visibilityAquatic.mat')
    else
        load('meteoAquatic.mat')
    end
                
    drdAAquatic=visualRangeAquatic;
    dVdAAquatic=visualRangeAquatic;
    visualVolumeAquatic=visualRangeAquatic;
%% Calculate visual volume based on parameters.m
    for i=1:size(visualRangeAquatic,3);
        for j=1:size(visualRangeAquatic,2);
            for k=1:size(visualRangeAquatic,1);
                visualVolumeAquatic(k,j,i)=integral3(f,0,visualRangeAquatic(k,j,i),...
                elevationMin,elevationMax,...
                azimuthMin,azimuthMax);
            end
        end
    end
%% Calculate visual range and volume derivatives    
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