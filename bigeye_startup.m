%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Startup file for directories to be added to the path
%%
%% Title                : A massive increase in visual range preceded the origin of terrestrial vertebrates
%% Authors              : Ugurcan Mugan, Malcolm A. MacIver
%% Authors' Affiliation : Northwestern University
%% This work is licensed under the Creative Commons Attribution 4.0 International License. 
%% To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
%% January 2017
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up globals
global BIGEYEROOT

BIGEYEROOT=pwd;
BIGEYEROOT=[BIGEYEROOT '/'];
BIGEYEROOT=strrep(BIGEYEROOT, '\','/');

alldirs = dir([BIGEYEROOT  'figs/']);

for j=1:length(alldirs)    
    if isdir([BIGEYEROOT  'figs/' alldirs(j).name]) && ~strcmp(alldirs(j).name,'archive')...
            && ~strcmp(alldirs(j).name,'.')&&~strcmp(alldirs(j).name,'cover_submission') &&~strcmp(alldirs(j).name,'..');
        fprintf('adding folder and all subfolders... %s\n',alldirs(j).name);
        addpath(genpath([BIGEYEROOT 'figs/' alldirs(j).name]));            
    end
end
%rmpath(genpath([BIGEYEROOT 'figs/archive']))
BIGEYEROOT=[BIGEYEROOT ,'figs/'];

clearvars -except BIGEYEROOT
