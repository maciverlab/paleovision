function calculatePercentChange(CONTRASTTHRESH)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Percent change between original parameter values and alternate values
%%
%% Title                : A massive increase in visual range preceded the origin of terrestrial vertebrates
%% Authors              : Ugurcan Mugan, Malcolm A. MacIver
%% Authors' Affiliation : Northwestern University
%% This work is licensed under the Creative Commons Attribution 4.0 International License. 
%% To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
%% January 2017
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global BIGEYEROOT
    if nargin<1
        CONTRASTTHRESH=1;
    end
    [visualRangeParam.Aquatic, visualVolumeParam.Aquatic, drdDParam.Aquatic, dVdDParam.Aquatic] = ...
        Aquatic_calcVolumegetDerSensitivity(CONTRASTTHRESH);
    [rawVisualRange.Aquatic, rawVisualVolume.Aquatic, rawdrdD.Aquatic, rawdVdD.Aquatic] =...
        Aquatic_calcVolumegetDer(CONTRASTTHRESH);
    [visualRangeParam.Aerial, visualVolumeParam.Aerial,drdDParam.Aerial, dVdDParam.Aerial]=...
        Aerial_calcVolumegetDerSensitivity(CONTRASTTHRESH);
    [rawVisualRange.Aerial, rawVisualVolume.Aerial, rawdrdD.Aerial,rawdVdD.Aerial]=...
        Aerial_calcVolumegetDer(CONTRASTTHRESH);

    allVisualRange.Aquatic=cat(3,[rawVisualRange.Aquatic(:,1,1),rawVisualRange.Aquatic(:,1,2)],visualRangeParam.Aquatic);
    allVisualRange.Aerial=cat(3,rawVisualRange.Aerial(:,1),visualRangeParam.Aerial);
    alldrdD.Aquatic=cat(3,[rawdrdD.Aquatic(:,1,1),rawdrdD.Aquatic(:,1,2)],drdDParam.Aquatic);
    alldrdD.Aerial=cat(3,rawdrdD.Aerial(:,1),drdDParam.Aerial);
    allVisualVolume.Aquatic=cat(3,[rawVisualVolume.Aquatic(:,1,1),rawVisualVolume.Aquatic(:,1,2)],visualVolumeParam.Aquatic);
    allVisualVolume.Aerial=cat(3,rawVisualVolume.Aerial(:,1),visualVolumeParam.Aerial);
    alldVdD.Aquatic=cat(3,[rawdVdD.Aquatic(:,1,1),rawdVdD.Aquatic(:,1,2)],dVdDParam.Aquatic);
    alldVdD.Aerial=cat(3,rawdVdD.Aerial(:,1),dVdDParam.Aerial);

    visualRange.Aquatic.NoChange=allVisualRange.Aquatic(:,:,1);
    visualRange.Aerial.NoChange=allVisualRange.Aerial(:,:,1);
    drdD.Aquatic.NoChange=alldrdD.Aquatic(:,:,1);
    drdD.Aerial.NoChange=alldrdD.Aerial(:,:,1);
    visualVolume.Aquatic.NoChange=allVisualVolume.Aquatic(:,:,1);
    visualVolume.Aerial.NoChange=allVisualVolume.Aerial(:,:,1);
    dVdD.Aquatic.NoChange=alldVdD.Aquatic(:,:,1);
    dVdD.Aerial.NoChange=alldVdD.Aerial(:,:,1);
    parameters={'X','q','Dt1','Dt2','d','M','contrast'};
    conditions={'Aquatic','Aerial'};

    for c=1:length(conditions)
        for i=1:length(parameters)
           visualRange.(conditions{c}).(parameters{i})=allVisualRange.(conditions{c})(:,:,i+1);
           percChangeRange.(conditions{c}).(parameters{i})=(visualRange.(conditions{c}).(parameters{i})-visualRange.(conditions{c}).NoChange)./visualRange.(conditions{c}).NoChange;
           meanPercChangeRange.(conditions{c}).(parameters{i})=mean(percChangeRange.(conditions{c}).(parameters{i}));

           drdD.(conditions{c}).(parameters{i})=alldrdD.(conditions{c})(:,:,i+1);
           percChangedrdD.(conditions{c}).(parameters{i})=(drdD.(conditions{c}).(parameters{i})-drdD.(conditions{c}).NoChange)./drdD.(conditions{c}).NoChange;
           meanPercChangedrdD.(conditions{c}).(parameters{i})=mean(percChangedrdD.(conditions{c}).(parameters{i}));

           visualVolume.(conditions{c}).(parameters{i})=allVisualVolume.(conditions{c})(:,:,i+1);
           percChangeVolume.(conditions{c}).(parameters{i})=(visualVolume.(conditions{c}).(parameters{i})-visualVolume.(conditions{c}).NoChange)./visualVolume.(conditions{c}).NoChange;
           meanPercChangeVolume.(conditions{c}).(parameters{i})=mean(percChangeVolume.(conditions{c}).(parameters{i}));

           dVdD.(conditions{c}).(parameters{i})=alldVdD.(conditions{c})(:,:,i+1);
           percChangedVdD.(conditions{c}).(parameters{i})=(dVdD.(conditions{c}).(parameters{i})-dVdD.(conditions{c}).NoChange)./dVdD.(conditions{c}).NoChange;
           meanPercChangedVdD.(conditions{c}).(parameters{i})=mean(percChangedVdD.(conditions{c}).(parameters{i}));
        end
    end

save([BIGEYEROOT 'fig04_visualrange/figure_sensitivity/percChange.mat'],'visualRange','percChangeRange','meanPercChangeRange',...
    'drdD','percChangedrdD','meanPercChangedrdD',...
    'visualVolume','percChangeVolume','meanPercChangeVolume',...
    'dVdD','percChangedVdD','meanPercChangedVdD');
   
   
   

 