function [MYr,varargout] = era2MYr(stage)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Convert era to the start and end MYr
%%
%% Title                : A massive increase in visual range preceded the origin of terrestrial vertebrates
%% Authors              : Chen Chen, Malcolm A. MacIver
%% Authors' Affiliation : Northwestern University
%% This work is licensed under the Creative Commons Attribution 4.0 International License. 
%% To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
%% January 2017
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT:
%%    stage - Input Stage, string, e.g. 'FAM(u)'
%% OUTPUT:
%%    MYr - [Year Start, Year End] (MYr ago)
%%    dur - (optional) Unit duration of the era (MYr), for each substage
%% USAGE:
%%    MYr = era2MYr('FAM(u)')
%% Handle Input
%(Data from ICS - International Chronostratigraphic Chart v 2015/01)

stgName = stage(1:3); % Name of the Stage

if length(stage) >= 5
    subStgLbl = stage(5);
else
    subStgLbl = '-';
end

switch stgName
    % From Late Silurian to Middle Devinian
    % (*** Substages and Acronyms are inferred ***)
    % (********** IT COULD BE WRONG!!!! **********)
    case 'GOR' % Gorstian - 1
        dur = 1.8;   % Duration (MYr)
        numOfSubStgs = 1;
        str = 427.4; % Start time (MYr ago)
        
    case 'LUD' % Ludfordian - 1
        dur = 2.6;   % Duration (MYr)
        numOfSubStgs = 1;
        str = 425.6; % Start time (MYr ago)
        
    case 'PRI' % Pridoli - 1 (This name is Epoch not Stage!)
        dur = 3.8;   % Duration (MYr)
        numOfSubStgs = 1;
        str = 423.0; % Start time (MYr ago)
        
    case 'LOC' % Lochkovian - 2
        dur = 8.4;   % Duration (MYr)
        numOfSubStgs = 2;
        str = 419.2; % Start time (MYr ago)

    case 'PRA' % Pragian - 2
        dur = 3.2;   % Duration (MYr)
        numOfSubStgs = 2;
        str = 410.8; % Start time (MYr ago)

    case 'EMS' % Emsian - 3
        dur = 14.3;   % Duration (MYr)
        numOfSubStgs = 3;
        str = 407.6; % Start time (MYr ago)

    case 'EIF' % Eifelian - 2
        dur = 5.6;   % Duration (MYr)
        numOfSubStgs = 2;
        str = 393.3; % Start time (MYr ago)

    % From Middle Devonian - Givetian (Start 387.7 MYr ago)
    %  To  Lower Jurassic  - Toarcian (182.7 MYr ago)
    case 'GIV' % Givetian - 2
        dur = 5;   % Duration (MYr)
        numOfSubStgs = 2;
        str = 387.7; % Start time (MYr ago)

    case 'FRS' % Frasnian - 3
        dur = 10.5;   % Duration (MYr)
        numOfSubStgs = 3;
        str = 382.7; % Start time (MYr ago)

    case 'FAM' % Famennian - 3
        dur = 13.3;   % Duration (MYr)
        numOfSubStgs = 3;
        str = 372.2; % Start time (MYr ago)

    case 'TOU' % Tournaisian - 3
        dur = 12.2;   % Duration (MYr)
        numOfSubStgs = 3;
        str = 358.9; % Start time (MYr ago)

    case 'VIS' % Visean - 3
        dur = 15.8;   % Duration (MYr)
        numOfSubStgs = 3;
        str = 346.7; % Start time (MYr ago)

    case 'SPK' % Serpukhovian - 2
        dur = 7.7;   % Duration (MYr)
        numOfSubStgs = 2;
        str = 330.9; % Start time (MYr ago)

    case 'BSH' % Bashkirian - 2
        dur = 8.0;   % Duration (MYr)
        numOfSubStgs = 2;
        str = 323.2; % Start time (MYr ago)

    case 'MOS' % Moscovian - 2
        dur = 8.2;   % Duration (MYr)
        numOfSubStgs = 2;
        str = 315.2; % Start time (MYr ago)

    case 'KAS' % Kasimovian - 1
        dur = 3.3;   % Duration (MYr)
        numOfSubStgs = 1;
        str = 307.0; % Start time (MYr ago)

    case 'GZE' % Gzelian - 1
        dur = 4.8;   % Duration (MYr)
        numOfSubStgs = 1;
        str = 303.7; % Start time (MYr ago)

    case 'ASS' % Asselian - 2
        dur = 3.9;   % Duration (MYr)
        numOfSubStgs = 2;
        str = 298.9; % Start time (MYr ago)

    case 'SAK' % Sakmarian - 2
        dur = 4.9;   % Duration (MYr)
        numOfSubStgs = 2;
        str = 295.0; % Start time (MYr ago)

    case 'ART' % Artinskian - 2
        dur = 6.6;   % Duration (MYr)
        numOfSubStgs = 2;
        str = 290.1; % Start time (MYr ago)

    case 'KUN' % Kungurian - 2
        dur = 11.2;   % Duration (MYr)
        numOfSubStgs = 2;
        str = 283.5; % Start time (MYr ago)

    case 'ROA' % Roadian - 1
        dur = 3.5;   % Duration (MYr)
        numOfSubStgs = 1;
        str = 272.3; % Start time (MYr ago)

    case 'WOR' % Wordian - 1
        dur = 3.7;   % Duration (MYr)
        numOfSubStgs = 1;
        str = 268.8; % Start time (MYr ago)

    case 'CAP' % Capitanian - 2
        dur = 5.3;   % Duration (MYr)
        numOfSubStgs = 2;
        str = 265.1; % Start time (MYr ago)

    case 'WUC' % Wuchiapingian - 2
        dur = 5.66;   % Duration (MYr)
        numOfSubStgs = 2;
        str = 259.8; % Start time (MYr ago)

    case 'CHX' % Changhsingian - 1
        dur = 1.97;   % Duration (MYr)
        numOfSubStgs = 1;
        str = 254.14; % Start time (MYr ago)

    case 'IND' % Induan - 1
        dur = 0.97;   % Duration (MYr)
        numOfSubStgs = 1;
        str = 252.17; % Start time (MYr ago)

    case 'OLE' % Olenekian - 2
        dur = 4.0;   % Duration (MYr)
        numOfSubStgs = 2;
        str = 251.2; % Start time (MYr ago)

    case 'ANS' % Anisian - 2
        dur = 5.2;   % Duration (MYr)
        numOfSubStgs = 2;
        str = 247.2; % Start time (MYr ago)

    case 'LAD' % Ladinian - 2
        dur = 5.0;   % Duration (MYr)
        numOfSubStgs = 2;
        str = 242.0; % Start time (MYr ago)

    case 'CRN' % Carnian - 2
        dur = 10.0;   % Duration (MYr)
        numOfSubStgs = 2;
        str = 237.0; % Start time (MYr ago)

    case 'NOR' % Norian - 3
        dur = 18.5;   % Duration (MYr)
        numOfSubStgs = 3;
        str = 227.0; % Start time (MYr ago)

    case 'RHT' % Rhaetian - 2
        dur = 7.2;   % Duration (MYr)
        numOfSubStgs = 2;
        str = 208.5; % Start time (MYr ago)

    case 'HET' % Hettangian - 1
        dur = 2.0;   % Duration (MYr)
        numOfSubStgs = 1;
        str = 201.3; % Start time (MYr ago)

    case 'SIN' % Sinemurian - 2
        dur = 8.5;   % Duration (MYr)
        numOfSubStgs = 2;
        str = 199.3; % Start time (MYr ago)

    case 'PLB' % Pliensbachian - 2
        dur = 8.1;   % Duration (MYr)
        numOfSubStgs = 2;
        str = 190.8; % Start time (MYr ago)

    case 'TOA' % Toarcian - 2
        dur = 8.6;   % Duration (MYr)
        numOfSubStgs = 2;
        str = 182.7; % Start time (MYr ago)

end

unitDur = double(dur / numOfSubStgs);
switch subStgLbl
    case 'l'
        MYr = [str,str-unitDur];
    case 'm'
        if numOfSubStgs == 2
            error('Stage %s DOES NOT have a middle substage!!',...
                stgName);
        end
        MYr = [str-unitDur,str-unitDur*2];
    case 'u'
        MYr = [str-unitDur*(numOfSubStgs-1),str-unitDur*numOfSubStgs];
    otherwise
        MYr = [str,str-dur];
end

MYr = MYr';
%% Output Duration if Requested
if nargout == 2
    varargout{1} = unitDur;
end