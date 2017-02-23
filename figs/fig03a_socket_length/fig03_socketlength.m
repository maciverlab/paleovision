function fig03_socket_length()
global BIGEYEROOT
close all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Statistical summary of raw socket length data based on grouping
%%
%% Title                : A massive increase in visual range preceded the origin of terrestrial vertebrates
%% Authors              : Malcolm A. MacIver, Ugurcan Mugan, Chen Chen
%% Authors' Affiliation : Northwestern University
%% This work is licensed under the Creative Commons Attribution 4.0 International License. 
%% To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
%% January 2017
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read or Load XLS Data
[~,~,XlsData.GS] = xlsread('bigEye.xlsx',1,'B3:C64');
[~,~,XlsData.Eye] = xlsread('bigEye.xlsx',1,'T3:T64');
[~,~,XlsData.Length] = xlsread('bigEye.xlsx',1,'P3:P64');
[~,~,XlsData.RefKey] = xlsread('bigEye.xlsx',1,'O3:O64');
[~,~,XlsData.startEnd] = xlsread('bigEye.xlsx',1,'D3:E64');

% Specify Needed Species Here
% Needed Species' Index in the Spreedsheet File
TFind = [1:21]' ;    % Tetrapodomorph Fish
STind = [25:62]' ; % Stem Tetrapods

global numTF
numTF=numel(TFind);

% Key column addresses in spreadsheet
stageBegin=4;    % start of stratigraphic period
stageEnd=5;      % end of stratigraphic period
genusCol = 2;    % genus
speciesCol = 3;  % species
refKeyCol = 15;  % bibtex ref key
lengthCol = 16;  % PPL
eyeCol = 20;     % OM

circleMax=60;    % max size of normalized eye, points

%% Convert to LaTeX Table
[TF_orbit_length_span, gsTF]=loc_ConvData(XlsData,TFind,lengthCol,refKeyCol,eyeCol,genusCol,speciesCol,stageBegin,stageEnd,'TFout.tex');
[ST_orbit_length_span, gsST]=loc_ConvData(XlsData,STind,lengthCol,refKeyCol,eyeCol,genusCol,speciesCol,stageBegin,stageEnd,'STout.tex');

% Hack to make table that includes both; still need above two lines since
% we need TF_orbit_length_span and ST* ... but really should rewrite

[all_orbit_length_span, allgsTF]=loc_ConvData(XlsData,[TFind; STind],lengthCol,refKeyCol,eyeCol,genusCol,speciesCol,stageBegin,stageEnd,'allout.tex');

% Data analysis on table information

% ratio of OM to PPL
datTF = (TF_orbit_length_span(:,1)./TF_orbit_length_span(:,2));
disp('Mean TF OM')
mean(TF_orbit_length_span(:,1))

disp('Mean ST OM')
mean(ST_orbit_length_span(:,1))

datST = (ST_orbit_length_span(:,1)./ST_orbit_length_span(:,2));

x=TF_orbit_length_span(:,1);
disp(['The median orbit length for TF was ' num2str(median(x)) '; the 1st quartile is ' ...
    num2str(prctile(x,25))])
disp(['; the 3rd quartile is ' num2str(prctile(x,75))])
disp(' ')

x=ST_orbit_length_span(:,1);
disp(['The median orbit length for ST was ' num2str(median(x)) '; the 1st quartile is ' ...
    num2str(prctile(x,25))])
disp(['; the 3rd quartile is ' num2str(prctile(x,75))])

%
%      ORBITDAT figure
%
%

fig_props.noYsubplots = 1;
fig_props.noXsubplots = 1;

fig_props.figW = 9;   % cm
fig_props.figH = 10;  % cm

fig_props.ml = 0.8;
fig_props.mt = 0.8;
fig_props.bottom_margin=2.5;

circleSize=5;
create_BE_figure
fig_props.sub_pW = fig_props.sub_pW-.5;
time_subsamp = 1;
time_limit = 0.4;
text_pos = [-5,2*time_limit/10,50];
text_color = [0 0 0];
text_size = 12;
pn = {'Color','FontSize','FontWeight',};
pv = {text_color,text_size,'bold'};

plotnoX = 1;
plotnoY = 1;
horiz_line_len=0.07;
ha3= create_BE_axes(plotnoX,plotnoY,fig_props);

% separate elpistostegalians

elpisto = {'Panderichthys','Tiktaalik','Elpistostege'};
noElpistoOrb=TF_orbit_length_span(:,1);
noElpistoSkl=TF_orbit_length_span(:,2);
elpistoIdx=zeros(length(TF_orbit_length_span),1);

for i=1:length(elpisto)
    idx=find(strcmp(gsTF(:,1), elpisto(i)) );
    elpistoIdx(idx)=1;
end

% break TF into elpistostegalians and non-epistos
noElpistoOrb(find(elpistoIdx))=[];
noElpistoSkl(find(elpistoIdx))=[];


onlyElpistoOrb=TF_orbit_length_span(:,1);
onlyElpistoSkl=TF_orbit_length_span(:,2);
onlyElpistoOrb(find(~elpistoIdx))=[];
onlyElpistoSkl(find(~elpistoIdx))=[];


% finned tetrapod minus Elpisto absolute size points
line(repmat(1/4,length(noElpistoOrb)),noElpistoOrb, ...
    'linestyle','none', ...
    'marker','o','markersize',circleSize,'markeredgecolor','black','markerfacecolor',[1 0 0])

meanOrb = mean(noElpistoOrb);
line([(1/4)-horiz_line_len (1/4)+horiz_line_len],[meanOrb meanOrb],'color',[1 0 0],'linewidth',1.5)
text((1/4)+horiz_line_len+0.01, meanOrb, [num2str(meanOrb,2) ' mm'])

hold on

% now plot only the Epistos
line(repmat(1/2,length(onlyElpistoOrb)),onlyElpistoOrb, ...
    'linestyle','none', ...
    'marker','o','markersize',circleSize,'markeredgecolor','black','markerfacecolor',[1 0 0])

meanOrb = mean(onlyElpistoOrb);
line([(1/2)-horiz_line_len (1/2)+horiz_line_len],[meanOrb meanOrb],'color',[1 0 0],'linewidth',1.5)
text((1/2)+horiz_line_len+0.01, meanOrb, [num2str(meanOrb,2) ' mm'])

% now only digited 

% first separate the adelospondyl-colosteid clade
secAqColor=[252 130 0]./255;
% Derive sets without seconarily aquatic animals:
secAquatic = {'Adelogyrinus','Adelospondylus','Acherontiscus','Colosteus','Greererpeton','Deltaherpeton'};
noSecAqOrb=ST_orbit_length_span(:,1);
noSecAqSkl=ST_orbit_length_span(:,2);
secAqIdx=zeros(length(ST_orbit_length_span),1);

for i=1:length(secAquatic)
    idx=find(strcmp(gsST(:,1), secAquatic(i)) );
    secAqIdx(idx)=1;
end

% break ST into non secondarily aquatic and aquatic for plotting and means
% all ST but no sec aquatic
noSecAqSkl(find(secAqIdx))=[];
noSecAqOrb(find(secAqIdx))=[];

% sec aquatic only
SecAqOrb=ST_orbit_length_span(:,1);
SecAqSkl=ST_orbit_length_span(:,2);
SecAqOrb(find(~secAqIdx))=[];
SecAqSkl(find(~secAqIdx))=[];

% all ST
Orb=ST_orbit_length_span(:,1);
Skl=ST_orbit_length_span(:,2);

%plot all ST without sec aq
line(repmat(3/4,length(noSecAqOrb),1),noSecAqOrb, ...
    'linestyle','none', ...
    'marker','o','markersize',circleSize,'markeredgecolor','black','markerfacecolor',[0 0 1])

% plot all ST without sec aq mean
meanOrb = mean(noSecAqOrb);
line([(3/4)-horiz_line_len (3/4)+horiz_line_len],[meanOrb meanOrb],'color',[0 0 1],'linewidth',1.5)
text((3/4)+horiz_line_len+0.01, meanOrb+.4, [num2str(meanOrb,2) ' mm'])

% plot all ST including sec aq mean
%plot all sec aq only
line(repmat(1,length(SecAqOrb),1),SecAqOrb, ...
    'linestyle','none', ...
    'marker','o','markersize',circleSize,'markeredgecolor','black','markerfacecolor',[0 0 1])

% plot only sec aq mean
meanOrb = mean(SecAqOrb);
line([(1)-horiz_line_len (1)+horiz_line_len],[meanOrb meanOrb],'color',[0 0 1],'linewidth',1.5)
text((1)+horiz_line_len+0.01, meanOrb, [num2str(meanOrb,2) ' mm'])

set(gca,'ylim',[0 70])
set(gca,'xlim',[0 1.25])
set(gca,'xtick', [1/4 1/2 3/4 1],'xticklabel',{'Finned','Finned-Transitional','Digited','Digited-Aquatic'})
set(gca,'XTickLabelRotation',-45');
ylabel('eye socket length (mm)')

[~,~,rawResGroups]=xlsread('BMresgroups.xlsx');
rawResGroups=rawResGroups(2:end,:);

indFinned=find(cellfun(@(x) strcmp(x,'finned'), rawResGroups(:,3),'UniformOutput',1));
finnedResVals=cell2mat(rawResGroups(indFinned,2));
indTransitional=find(cellfun(@(x) strcmp(x,'transitional'), rawResGroups(:,3),'UniformOutput',1));
transitionalResVals=cell2mat(rawResGroups(indTransitional,2));
indDigited=find(cellfun(@(x) strcmp(x,'digited'), rawResGroups(:,3),'UniformOutput',1));
digitedResVals=cell2mat(rawResGroups(indDigited,2));
indReturn=find(cellfun(@(x) strcmp(x,'return'), rawResGroups(:,3),'UniformOutput',1));
returnResVals=cell2mat(rawResGroups(indReturn,2));

% RELATIVE SIZES NOW
% finned tetrapod minus Elpisto absolute size points
    
relsizeTF= 100.*(noElpistoOrb./noElpistoSkl);

%only elpisto
rootFinned=0.02058466;

%only transitional
rootTransitional=0.1725637;

% digited tetrapods minus adelospondyl-colosteid clade
rootDigited=0.2038828;

% only sec aq
rootReturn=0.04520924;

filename=[BIGEYEROOT 'fig03_socket_length/orbitData.pdf'];
print(filename,'-painters','-dpdf','-r600');

x=TF_orbit_length_span(:,1);
disp(['The median orbit length for TF was ' num2str(median(x)) '; the 1st quartile is ' ...
     num2str(prctile(x,25))])
disp(['; the 3rd quartile is ' num2str(prctile(x,75))])
disp(' ')

OM_TF=x;

x=ST_orbit_length_span(:,1);
disp(['The median orbit length for ST was ' num2str(median(x)) '; the 1st quartile is ' ...
     num2str(prctile(x,25))])
disp(['; the 3rd quartile is ' num2str(prctile(x,75))])
 
OM_ST=x;

save([BIGEYEROOT,'fig03_socket_length/FinnedDigitedOrbitLength.mat'],'onlyElpistoOrb','noElpistoOrb','noSecAqOrb','SecAqOrb')

save([BIGEYEROOT,'fig03_socket_length/OM_TF_ST.mat'],'OM_TF','OM_ST')


% Format stat result text for paper
formattedstringA = [ ...
    'The ' ...
    'finned tetrapods had eye sockets that were $' num2str(mean(TF_orbit_length_span(:,1)),2) ...
    ' \pm ' num2str(std(TF_orbit_length_span(:,1)),2) '$~mm long (all numbers mean $\pm$ standard deviation), ' ...
    'while the digited tetrapods had eye sockets that were ' num2str(mean(ST_orbit_length_span(:,1)),2) ...
    ' $\pm$ ' num2str(std(ST_orbit_length_span(:,1)),2) '~mm long. Relative eye socket length, expressed ' ...
    'as a percentage of skull length, ' ...
    'was ' num2str([100.*mean(datTF)],2) ' $\pm$ ' num2str(std(100.*datTF),2)  ...
    '\% ($N=' num2str(length(datTF)) '$) for the' ...
    ' finned tetrapods, and ' num2str([100.*mean(datST)],2) ' $\pm$  ' num2str(std(100.*datST),2)  ...
    '\% ($N=' num2str(length(datST)) '$) for the digited tetrapods.'];

formattedstringB = ['The finned tetrapods ' ... 
    'exclusive of the transitional elpistostegalians had eye sockets that were $' num2str(mean(noElpistoOrb),2) ' \pm ' ...
    num2str(std(noElpistoOrb),2) '$~mm, and  $' ...
    num2str(mean([100.*(noElpistoOrb./noElpistoSkl)]),2) ' \pm ' num2str(std([100.*(noElpistoOrb./noElpistoSkl)]),2) '$\%. ' ... 
    'Elpistostegalians ($N=' num2str(length(onlyElpistoOrb)) '$) alone were $' num2str(mean(onlyElpistoOrb),2) ...
    ' \pm ' num2str(std(onlyElpistoOrb),2) '$~mm, and $' ...
    num2str(mean([100.*(onlyElpistoOrb./onlyElpistoSkl)]),2) ' \pm ' num2str(std([100.*(onlyElpistoOrb./onlyElpistoSkl)]),2) '$\%. ' ... 
    'The digited tetrapods without the secondarily aquatic colosteid-adelospondyls had socket lengths of $' num2str(mean(noSecAqOrb),2) ...
    ' \pm ' num2str(std(noSecAqOrb),2) '$~mm, while the colosteid-adelospondyls ' ...
    '($N=' num2str(length(SecAqOrb)) '$) alone had sockets that were $' ...
  num2str(mean(SecAqOrb),2) ...
  ' \pm ' num2str(std(SecAqOrb),2) '$~mm. ' ...
  'The corresponding relative socket lengths were $' ... 
    num2str(mean([100.*(noSecAqOrb./noSecAqSkl)]),2) ' \pm ' num2str(std([100.*(noSecAqOrb./noSecAqSkl)]),2) '$\% and $' ...
    num2str(mean([100.*(SecAqOrb./SecAqSkl)]),2) ...
    ' \pm ' num2str(std([100.*(SecAqOrb./SecAqSkl)]),2) '$\%, respectively.']
% CONTINUE WRITING

fid = fopen('stat_string1.tex','w');
fprintf(fid, '%s\n\n', formattedstringA);
fprintf(fid, '%s\n\n', formattedstringB);
fclose(fid);


disp('Done!!');

end



function [orbit_length_span, gs] = loc_ConvData(gdat,speciesIndex,lengthCol, refkeyCol,eyeCol,genusCol,speciesCol,stageBegin,stageEnd,filename)
global numTF

refkey=[];
validInd=[];
count=1;
for i = 1:numel(speciesIndex)
    tempkey=gdat.RefKey{speciesIndex(i)};
    
%     if isstr(tempkey)
%         refkey{count} = strrep(gdat.RefKey{speciesIndex(i)},'*','');
%         count=count+1;
%         validInd=[validInd; speciesIndex(i)];
%     end

    if isstr(tempkey)
        refkey{speciesIndex(i)} = strrep(gdat.RefKey{speciesIndex(i)},'*','');
        count=count+1;
        validInd=[validInd; speciesIndex(i)];
    end

end


%% Convert RefData to Cited Reference
citeRef = refkey;
for i = 1:numel(citeRef)
    if length(citeRef{i}) > 4
        citeRef{i} = ['\cite{',citeRef{i},'}'];
    end
end

%% Assemble Format Genus and Species

gs = gdat.GS(validInd,:);
gsFormatted = cell(size(gs,1),1);

startend = gdat.startEnd(validInd,:); % N x 2, with start stage and end stage

span = [];
for i = 1: length(startend)
    MYaS = era2MYr_ICS2015_v01(cell2mat(startend(i,1))); % Start and end of start stage
    MYaE = era2MYr_ICS2015_v01(cell2mat(startend(i,2))); % Start and end of end stage
    % era2MYr returns the
    span = [span; MYaS(1) MYaE(2)]; % Start of start stage and end of end stage
    
end


for i = 1:size(gs,1)
    if iscellstr(gs(i,1)) && iscellstr(gs(i,2)) % species AND genus
        gsFormatted{i,:} = strcat('\textit{',strjoin(gs(i,:)),'}');
    else
        gsFormatted{i,:} = strcat('\textit{',strjoin(gs(i,1)),'}'); % genus only
    end
end

osf = cell(size(gs,1),1);
%% Get Average AP Data
orb2skull = cell2mat(gdat.Eye(validInd))./cell2mat(gdat.Length(validInd));
strOrbSkull = num2str(100.*orb2skull,'%-3.0f');

for i= 1:size(gs,1)
    osf{i,:}=strOrbSkull(i,:); % orbit to skull fraction
end

% hack for full table
% sort the TF and ST separately

if length(startend)>40
    [A,sortedIDX]=sort(gsFormatted(1:numTF,:)) ;
    orbits=gdat.Eye(validInd(1:numTF));
    skulls=gdat.Length(validInd(1:numTF));
    cites=[citeRef(validInd(1:numTF))]';
    osfA=[osf(1:numTF)];
    
    AveAPTF=[A, ...
        repmat({'Finned'},length(A),1), ...
        orbits(sortedIDX),...    % Average AP (mm)
        skulls(sortedIDX),... % Skull Length (mm)
        osfA(sortedIDX),... % orbit skull fraction
        cites(sortedIDX)]; % AP-src BibTeX Reference Key
    
    [B,sortedIDX]=sort(gsFormatted(numTF+1:end,:)) ;
    orbits=gdat.Eye(validInd(numTF+1:end));
    skulls=gdat.Length(validInd(numTF+1:end));
    cites=[citeRef(validInd(numTF+1:end))]';
     osfB=[osf(numTF+1:end)];
    
    AveAPST=[B, ...
        repmat({'Digited'},length(B),1), ...
        orbits(sortedIDX),...    % Average AP (mm)
        skulls(sortedIDX),... % Skull Length (mm)
        osfB(sortedIDX),... % orbit skull fraction
        cites(sortedIDX)]; % AP-src BibTeX Reference Key
    
    AveAP =[AveAPTF; AveAPST];
    
    %% Convert to LaTeX Table Format
    columnLabels = {'Taxon'; 'Group'; 'Orbit Length (mm)'; 'Skull Length (mm)'; 'Orbit/Skull, \%'; 'Reference'};
    
else
    
    [B,sortedIDX]=sort(gsFormatted);
    orbits=gdat.Eye(validInd);
    skulls=gdat.Length(validInd);
    cites=citeRef';
    
    AveAP = [B,... % Genus, species
        orbits(sortedIDX),...    % Average AP (mm)
        skulls(sortedIDX),... % Skull Length (mm)
        osf(sortedIDX),... % Skull Length (mm)
        cites(sortedIDX)]; % AP-src BibTeX Reference Key
    
    %% Convert to LaTeX Table Format
    columnLabels = {'Taxon'; 'Orbit Length (mm)'; 'Skull Length (mm)'; 'Orbit/Skull, \%'; 'Reference'};
end
orbit_length = [gdat.Eye(validInd), gdat.Length(validInd)];

% (Optional) Convert output from string to double
orbit_length = cell2mat(orbit_length);

[orbit_length_span] = [orbit_length span];


end

function matrix2latex(matrix, filename, varargin)
% function: matrix2latex(...)
% Author:   M. Koehler
% Contact:  koehler@in.tum.de
% Version:  1.1
% Date:     May 09, 2004

% This software is published under the GNU GPL, by the free software
% foundation. For further reading see: http://www.gnu.org/licenses/licenses.html#GPL

% Usage:
% matrix2late(matrix, filename, varargs)
% where
%   - matrix is a 2 dimensional numerical or cell array
%   - filename is a valid filename, in which the resulting latex code will
%   be stored
%   - varargs is one ore more of the following (denominator, value) combinations
%      + 'rowLabels', array -> Can be used to label the rows of the
%      resulting latex table
%      + 'columnLabels', array -> Can be used to label the columns of the
%      resulting latex table
%      + 'alignment', 'value' -> Can be used to specify the alginment of
%      the table within the latex document. Valid arguments are: 'l', 'c',
%      and 'r' for left, center, and right, respectively
%      + 'format', 'value' -> Can be used to format the input data. 'value'
%      has to be a valid format string, similar to the ones used in
%      fprintf('format', value);
%      + 'size', 'value' -> One of latex' recognized font-sizes, e.g. tiny,
%      HUGE, Large, large, LARGE, etc.
%
% Example input:
%   matrix = [1.5 1.764; 3.523 0.2];
%   rowLabels = {'row 1', 'row 2'};
%   columnLabels = {'col 1', 'col 2'};
%   matrix2latex(matrix, 'out.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');
%
% The resulting latex file can be included into any latex document by:
% /input{out.tex}
%
% Enjoy life!!!

rowLabels = [];
colLabels = [];
alignment = 'l';
format = [];
textsize = [];
if (rem(nargin,2) == 1 || nargin < 2)
    error('matrix2latex: ', 'Incorrect number of arguments to %s.', mfilename);
end

okargs = {'rowlabels','columnlabels', 'alignment', 'format', 'size'};
for j=1:2:(nargin-2)
    pname = varargin{j};
    pval = varargin{j+1};
    k = strmatch(lower(pname), okargs);
    if isempty(k)
        error('matrix2latex: ', 'Unknown parameter name: %s.', pname);
    elseif length(k)>1
        error('matrix2latex: ', 'Ambiguous parameter name: %s.', pname);
    else
        switch(k)
            case 1  % rowlabels
                rowLabels = pval;
                if isnumeric(rowLabels)
                    rowLabels = cellstr(num2str(rowLabels(:)));
                end
            case 2  % column labels
                colLabels = pval;
                if isnumeric(colLabels)
                    colLabels = cellstr(num2str(colLabels(:)));
                end
            case 3  % alignment
                alignment = lower(pval);
                if alignment == 'right'
                    alignment = 'r';
                end
                if alignment == 'left'
                    alignment = 'l';
                end
                if alignment == 'center'
                    alignment = 'c';
                end
                if alignment ~= 'l' && alignment ~= 'c' && alignment ~= 'r'
                    alignment = 'l';
                    warning('matrix2latex: ', 'Unkown alignment. (Set it to \''left\''.)');
                end
            case 4  % format
                format = lower(pval);
            case 5  % format
                textsize = pval;
        end
    end
end

fid = fopen(filename, 'w');

width = size(matrix, 2);
height = size(matrix, 1);



if isnumeric(matrix)
    matrix = num2cell(matrix);
    for h=1:height
        for w=1:width
            if(~isempty(format))
                matrix{h, w} = num2str(matrix{h, w}, format);
            else
                matrix{h, w} = num2str(matrix{h, w});
            end
        end
    end
end

if(~isempty(textsize))
    fprintf(fid, '\\begin{%s}', textsize);
end

fprintf(fid, '\\begin{tabular}{|');

if(~isempty(rowLabels))
    fprintf(fid, 'l|');
end
% M. MacIver 3/18/2016
% Customized this to be left margin only, fixed width at 4.5 cm
% using package array
fprintf(fid, '%s|', 'L{4.5cm}');
for i=2:width-1
    fprintf(fid, '%c|', alignment);
end
% M. MacIver 3/18/2016
% Customized this to be left margin only, fixed width at 4.5 cm
% using package array
fprintf(fid, '%s|', 'L{4.5cm}')
fprintf(fid, '}\r\n');

fprintf(fid, '\\hline\r\n');

if(~isempty(colLabels))
    if(~isempty(rowLabels))
        fprintf(fid, '&');
    end
    for w=1:width-1
        fprintf(fid, '\\textbf{%s}&', colLabels{w});
    end
    fprintf(fid, '\\textbf{%s}\\\\\\hline\r\n', colLabels{width});
end

for h=1:height
    if(~isempty(rowLabels))
        fprintf(fid, '\\textbf{%s}&', rowLabels{h});
    end
    for w=1:width-1
        
        %     if ~isnan(str2double(matrix{h, w}))
        if ~isstr(matrix{h,w})
            tmpNum = matrix{h, w};
            if(~isempty(format))
                matrix{h, w} = num2str(tmpNum, format);
            end
        end
        fprintf(fid, '%s&', matrix{h, w});
      end
    fprintf(fid, '%s\\\\\\hline\r\n', matrix{h, width});
end

fprintf(fid, '\\end{tabular}\r\n');

if(~isempty(textsize))
    fprintf(fid, '\\end{%s}', textsize);
end

fclose(fid);
end

