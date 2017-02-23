function ha = create_axes(plotnoX,plotnoY,fig_props)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Custom subplot function
%%
%% Title                : A massive increase in visual range preceded the origin of terrestrial vertebrates
%% Authors              : Chen Chen, Ugurcan Mugan, Malcolm A. MacIver
%% Authors' Affiliation : Northwestern University
%% This work is licensed under the Creative Commons Attribution 4.0 International License. 
%% To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
%% January 2017
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ha = axes('unit','centimeters','Position',[fig_props.left_margin+(fig_props.noXsubplots-plotnoX)*fig_props.sub_pW+...
%     (fig_props.noXsubplots-plotnoX)*fig_props.ml,fig_props.bottom_margin+(fig_props.noYsubplots-plotnoY)*fig_props.sub_pH+...
%     (fig_props.noYsubplots-plotnoY)*fig_props.mt, fig_props.sub_pW, fig_props.sub_pH],'FontSize',8,'Box','off');

ha = axes('unit','centimeters','Position',[fig_props.left_margin+(plotnoX-1)*fig_props.sub_pW+...
    (plotnoX-1)*fig_props.ml,fig_props.bottom_margin+(fig_props.noYsubplots-plotnoY)*fig_props.sub_pH+...
    (fig_props.noYsubplots-plotnoY)*fig_props.mt, fig_props.sub_pW, fig_props.sub_pH],'FontSize',8,'Box','off',...
    'fontname','helvetica');