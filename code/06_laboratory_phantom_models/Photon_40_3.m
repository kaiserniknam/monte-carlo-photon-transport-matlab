function [] = Photon_40_3 ()
% Repository group: 06_laboratory_phantom_models
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: Replication of Arnab's Experiment: Dissolving x grams of TiO₂ in a solution of 10 g Sodium Polyacrylate (PAS) in 440 mL of water.
% general analysis: Arnab data

clc
close all

% size and code properties
set_of_Cs = [0  1   4   7];
set_of_ds = [1	2	3	4	5	6	7	8];
set_of_V = [...
    nan      1.76	    1.18        0.94    	0.66    	0.53        0.56    	0.52;
    nan      8.29       2.32    	0.80        0.40        0.28    	0.12    	0.12;
    nan      5.88	    0.96	    0.32    	0.21    	nan	        nan     	nan;
    nan      2.96	    0.72	    nan	        nan	        nan	        nan	       nan];

for i_C = 2:length(set_of_Cs)
    % read dbase
    n_C = set_of_Cs(i_C);
    disp(['C = ',num2str(n_C),' g/L'])
    TheLegend = ['C = ',num2str(n_C),' g/L'];

    % I vs. d
    TheLineStyle = '-'; TheMarker = 's'; x_label = 'd (cm)'; y_label = '\DeltaOD (a.u.)';
    subplot(1,2,1), xlabel(x_label), ylabel(y_label), title('OD vs d'),   axis square
    plot(set_of_ds,-log(set_of_V(i_C,:)./set_of_V(2,:)),           'DisplayName',TheLegend,'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','Color',get_color(set_of_Cs(i_C)),'MarkerSize',16,'LineWidth',2), hold on
    set(gca,'fontsize',24), legend('show','Location','northeast')
    y_label = '\DeltaOD/d (cm^{-1})';
    subplot(1,2,2), xlabel(x_label), ylabel(y_label), title('OD/d vs d'), axis square
    plot(set_of_ds,-log(set_of_V(i_C,:)./set_of_V(2,:))./set_of_ds,'DisplayName',TheLegend,'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','Color',get_color(set_of_Cs(i_C)),'MarkerSize',16,'LineWidth',2), hold on
    set(gca,'fontsize',24), legend('show','Location','northeast')
    clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label Therloess blnFit y_ref





    clearvars t_db the_filename mua mus ind_diff ind_trns n_g TheColor
end
end

function [out] = get_color(index)
if     index==1
    out = [0.0000 0.4470 0.7410];
elseif index==2
    out = [0.8500 0.3250 0.0980];
elseif index==3
    out = [0.9290 0.6940 0.1250];
elseif index==4
    out = [0.4940 0.1840 0.5560];
elseif index==5
    out = [0.4660 0.6740 0.1880];
elseif index==6
    out = [0.3010 0.7450 0.9330];
elseif index==7
    out = [0.6350 0.0780 0.1840];
else
    out = [0.0000 0.0000 0.0000];
end
end
