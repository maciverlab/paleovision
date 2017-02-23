function figExt06_sensitivity(CONTRASTTHRESH)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2 panel figure depicting aquatic visual range for upward and horizontal viewing under different water conditions
%%
%% Title                : A massive increase in visual range preceded the origin of terrestrial vertebrates
%% Authors              : Ugurcan Mugan, Malcolm A. MacIver
%% Authors' Affiliation : Northwestern University
%% This work is licensed under the Creative Commons Attribution 4.0 International License. 
%% To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
%% January 2017
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global BIGEYEROOT
%% Initialize variables
    close all;   
    load OM_TF_ST.mat
    load FinnedDigitedOrbitLength.mat
    
    pupil_TF = [mean(noElpistoOrb)-std(noElpistoOrb) mean(noElpistoOrb)+std(noElpistoOrb)].*0.449;
    pupil_ST = [mean(noSecAqOrb)-std(noSecAqOrb) mean(noSecAqOrb)+std(noSecAqOrb)].*0.449;
    finnedpupil=mean(noElpistoOrb)*.449;
    digitedpupil=mean(noSecAqOrb)*.449;
%% Fined required data files and ask whether to run again    
    [e,em]=fileExists;  
    while(~all(e))
        notFound=find(e==0);
        warning('Not all *.mat files required are found some are going to be re-run');
        pause(1)
        for i=1:length(notFound)
            fprintf('running %s\n',em{notFound(i)});
            run(em{notFound(i)});
        end
        [e,em]=fileExists;
    end
        
    h=warndlg({'All of the code takes about 1-2hrs to run'},'Warning!');
    waitfor(h);
    choice=questdlg({'All the required *.mat files found!',...
        'Re-run the code?'},'code re-run','yes','no','no');
    if strcmp(choice,'yes')
        for i=1:length(em)
            fprintf('running %s\n',em{i});
            run(em{i})
        end
    end
    
    if nargin<1
        CONTRASTTHRESH=1;
    end
    if CONTRASTTHRESH
        load('visibilityAquaticSensitivity.mat');
        load('visibilityAquatic.mat');
    else
        load('meteoAquaticSensitivity.mat');
        load('meteoAquatic.mat');
    end
    
    visualRangeSensitivity=[smooth(visualRangeAquaticSensitivity(:,1,1),7),smooth(visualRangeAquaticSensitivity(:,2,1),7),smooth(visualRangeAquaticSensitivity(:,3,1),7),...
        smooth(visualRangeAquaticSensitivity(:,4,1),7),smooth(visualRangeAquaticSensitivity(:,1,2),7),smooth(visualRangeAquaticSensitivity(:,2,2),7),...
        smooth(visualRangeAquaticSensitivity(:,3,2),7),smooth(visualRangeAquaticSensitivity(:,4,2),7)];
    linewidthdef=2;

%% Plot figure 6 panel
    fig_props.noYsubplots = 1;
    fig_props.noXsubplots = 2;

    fig_props.figW = 25;   % cm
    fig_props.figH = 10;  % cm

    fig_props.ml = 0.8;
    fig_props.mt = 0.8;
    fig_props.bottom_margin=2;
    create_BE_figure
    fig_props.sub_pW = fig_props.sub_pW-.5;

    colors=get(gca,'colororder'); clf;
    x=17;
% Sensitivity Upward Viewing
    plotnoX= 1;
    plotnoY= 1;
    ha1 = create_BE_axes(plotnoX,plotnoY,fig_props);
    hl1.A=line('XData',DrangeAquatic*10^3,'YData',visualRangeAquaticSensitivity(:,1),...
        'color',colors(1,:),'linewidth',linewidthdef);
    hold on;
    for i=2:4
        str=char(i+'A'-1);
        hl1.(str)=line('XData',DrangeAquatic*10^3,'YData',visualRangeAquaticSensitivity(:,i),...
            'color',colors(i,:),'linewidth',linewidthdef);
    end
    hl1.E=line('XData',DrangeAquatic*10^3,'YData',visualRangeAquatic(:,1,1),...
        'color',colors(5,:),'linewidth',linewidthdef);
    ylim1=get(gca,'ylim');
    line([finnedpupil,finnedpupil],[ylim1(1),ylim1(2)],...
        'linewidth',linewidthdef,'color','r','linestyle',':');
    line([digitedpupil,digitedpupil],[ylim1(1),ylim1(2)],...
            'linewidth',linewidthdef,'color','b','linestyle',':');
    ylabel('\bfvisual range (\itr) \rm\bf(m)','interpreter','tex',...
        'fontsize',12,'fontname','helvetica');
    xlabel('\bfpupil diameter (\itD) \rm\bf(mm)','interpreter','tex',...
        'fontsize',12,'fontname','helvetica');
    text(x,1.8*ylim1(2)/5,'\bfupward viewing','interpreter',...
        'tex','fontsize',13,'fontname','helvetica');
    axis square
        
% Sensitivity Horizontal Viewing
    plotnoX=2;
    plotnoY=1;
    ha2=create_BE_axes(plotnoX,plotnoY,fig_props);
    hl2.A=line('XData',DrangeAquatic*10^3,'YData',visualRangeAquaticSensitivity(:,5),...
        'color',colors(1,:),'linewidth',linewidthdef);
    hold on;
    for i=6:8
        str=char(i+'A'-4);
        hl1.(str)=line('XData',DrangeAquatic*10^3,'YData',visualRangeAquaticSensitivity(:,i),...
            'color',colors(i-4,:),'linewidth',linewidthdef);
    end
    hl1.E=line('XData',DrangeAquatic*10^3,'YData',visualRangeAquatic(:,1,2),...
        'color',colors(5,:),'linewidth',linewidthdef);
    line([finnedpupil,finnedpupil],[ylim1(1),ylim1(2)],...
        'linewidth',linewidthdef,'color','r','linestyle',':');
    line([digitedpupil,digitedpupil],[ylim1(1),ylim1(2)],...
            'linewidth',linewidthdef,'color','b','linestyle',':');
    ylabel('\bfvisual range (\itr) \rm\bf(m)','interpreter','tex',...
        'fontsize',12,'fontname','helvetica');
    xlabel('\bfpupil diameter (\itD) \rm\bf(mm)','interpreter','tex',...
        'fontsize',12,'fontname','helvetica');
    text(x,1.8*ylim1(2)/5,'\bfhorizontal viewing','interpreter','tex',...
        'fontsize',13,'fontname','helvetica');
    axis square
    
hLegend=legend('high turbidity @8 m','clear @8 m',...
    'absorption  dominated @8 m','scattering dominated @8 m','baseline @8 m');
set(hLegend,'box','off'); set(hLegend,'interpreter','tex'); 
set(hLegend,'fontsize',11,'fontname','helvetica'); set(hLegend,'orientation','horizontal')
rect=[0.375 0 0.25 0.1]; set(hLegend,'Position',rect)

filename=[BIGEYEROOT 'figExt06_sensitivity/figures/water_sensitivity.pdf'];
print(filename,'-painters','-dpdf','-r600');
    
function [e,em]=fileExists
    e1={exist('meteoAquaticSensitivity.mat','file')==2, 'Aquatic_firingThreshSensitivity.m'};
    e2={exist('visibilityAquaticSensitivity.mat','file')==2, 'Aquatic_contrastThreshSensitivity.m'};
    e3={exist('meteoAquatic.mat','file')==2, 'Aquatic_firingThresh.m'};
    e4={exist('visibilityAquatic.mat','file')==2, 'Aquatic_contrastThreshold.m'};
    e=[e1{1},e2{1},e3{1},e4{1}];
    em={e1{2},e2{2},e3{2},e4{2}};