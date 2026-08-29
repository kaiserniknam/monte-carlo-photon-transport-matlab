function [] = Photon_49_2 ()
% Repository group: 06_laboratory_phantom_models
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: Replication of Arnab's Experiment: Dissolving x grams of TiO₂ in a solution of 10 g Sodium Polyacrylate (PAS) in 440 mL of water.
% the same as Photon_40, but between 1 & 3
% general analysis: Arnab data for Optics Letter

clc
close all
format long
set(0,'DefaultFigureWindowStyle','docked')

% size and code properties
set_of_Ceq = [3,6];
% a = 1.4; set_of_g = [1.4 2.5];
a = 1.5; set_of_g = [1.5 2.6]; % equal to 4 & 7 grams
% a = 1.6; set_of_g = [1.6 2.8];

N_bins = 151;
pq_db = load('DB/Photon_49_1.mat');
Lx = 29.1; Ly = 29.1; Lz = 10.0; % The size of the computational domain (in cm) to prevent reflections from the boundaries.

set_of_ds = [1	2	3	4	5	6	7	8];
set_of_V = [...
            nan     1.76	    1.18        0.94    	0.66    	0.53        0.56    	0.52;
            nan     8.29        2.32    	0.80        0.40        0.28    	0.12    	0.12;
            nan     5.88	    0.96	    0.32    	0.21    	nan	        nan     	nan;
            nan     2.96	    0.72	    nan	        nan	        nan	        nan         nan];
set_of_Is = nan(length(set_of_g),N_bins);

for i_g = 1:length(set_of_g)
    n_g = set_of_g(i_g);
    % read dbase
    the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_49_g_',sprintf('%.2f',set_of_g(i_g)),'.mat'];
    t_db = load(the_filename); clearvars the_filename
    % read refc = 1 gram
    the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_49_g_',sprintf('%.2f',1.00000000000),'.mat'];
    t_rf = load(the_filename); clearvars the_filename

    TheLegend = ['\DeltaC = ',num2str(set_of_Ceq(i_g)),' g/L: True DPF'];
    disp(['g = ',num2str(n_g)])
    disp(['mu_a=',num2str(t_db.mua),', mu_s=',num2str(t_db.mus),', n=',num2str(t_db.n),', & g = ',num2str(t_db.g)])
    disp(['mu_a=',num2str(0.01),', mu_s=',num2str(39.835.*n_g),', n=',num2str(1.333),', & g = ',num2str(0.9)])

    % removing noisy points from db
    u = round(-log(t_db.w)./t_db.s,4); u_unique = unique(u); freq = nan(size(u_unique));
    for i_f = 1:length(u_unique), freq(i_f) = sum(u==u_unique(i_f)); end
    [~,i_fx] = max(freq);
    t_db.x = t_db.x(u==u_unique(i_fx));
    t_db.y = t_db.y(u==u_unique(i_fx));
    t_db.z = t_db.z(u==u_unique(i_fx));
    t_db.d = t_db.d(u==u_unique(i_fx));
    t_db.s = t_db.s(u==u_unique(i_fx));
    t_db.w = t_db.w(u==u_unique(i_fx));
    t_db.c = t_db.c(u==u_unique(i_fx));
    t_db.a = t_db.a(u==u_unique(i_fx));
    disp(num2str([length(u) sum(u==u_unique(i_fx)) sum(u==u_unique(i_fx))/length(u)*100]))
    t_db.no_of_photons = t_db.no_of_photons-sum(u~=u_unique(i_fx));
    clearvars u u_unique freq i_f i_fx

    % 1-D sorting
    d_diff_edges = linspace(0,sqrt((Lx/2).^2+(Ly/2).^2),N_bins+1); [~,~,ind_diff] = histcounts(sqrt(t_db.x(t_db.c==0).^2+t_db.y(t_db.c==0).^2+t_db.z(t_db.c==0).^2),d_diff_edges); % d_diffuse bins
    clearvars d_trns_edges d_diff_edges



    % removing noisy points from rf
    u = round(-log(t_rf.w)./t_rf.s,4); u_unique = unique(u); freq = nan(size(u_unique));
    for i_f = 1:length(u_unique), freq(i_f) = sum(u==u_unique(i_f)); end
    [~,i_fx] = max(freq);
    t_rf.x = t_rf.x(u==u_unique(i_fx));
    t_rf.y = t_rf.y(u==u_unique(i_fx));
    t_rf.z = t_rf.z(u==u_unique(i_fx));
    t_rf.d = t_rf.d(u==u_unique(i_fx));
    t_rf.s = t_rf.s(u==u_unique(i_fx));
    t_rf.w = t_rf.w(u==u_unique(i_fx));
    t_rf.c = t_rf.c(u==u_unique(i_fx));
    t_rf.a = t_rf.a(u==u_unique(i_fx));
    disp(num2str([length(u) sum(u==u_unique(i_fx)) sum(u==u_unique(i_fx))/length(u)*100]))
    t_rf.no_of_photons = t_rf.no_of_photons-sum(u~=u_unique(i_fx));
    clearvars u u_unique freq i_f i_fx

    % 1-D sorting
    d_diff_edges = linspace(0,sqrt((Lx/2).^2+(Ly/2).^2),N_bins+1); [~,~,ind_refc] = histcounts(sqrt(t_rf.x(t_rf.c==0).^2+t_rf.y(t_rf.c==0).^2+t_rf.z(t_rf.c==0).^2),d_diff_edges); % d_diffuse bins for rf
    d_centers = 1/2*(d_diff_edges(1:end-1)+d_diff_edges(2:end-0));
    clearvars d_trns_edges d_diff_edges






    % I vs. d
    i_fig = 1; fun_x = @mean; fun_y = @sum; TheOutFun = @(x)(-log(x)); TheLineStyle = '-'; TheMarker = 'none'; x_label = 'd (cm)'; y_label = '\DeltaOD (a.u.)';
    index_in = ind_diff; TheCode = 0; TheTitle = '\DeltaOD vs. d';
    x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
    y_temp = t_db.w(t_db.c==TheCode)./t_db.no_of_photons;                                            y_bind = accumarray(index_in,y_temp,[],fun_y,nan);

    index_in = ind_refc; TheCode = 0;
    x_temp = sqrt(t_rf.x(t_rf.c==TheCode).^2+t_rf.y(t_rf.c==TheCode).^2+t_rf.z(t_rf.c==TheCode).^2); x_refc = accumarray(index_in,x_temp,[],fun_x,nan);
    y_temp = t_rf.w(t_rf.c==TheCode)./t_rf.no_of_photons;                                            y_refc = accumarray(index_in,y_temp,[],fun_y,nan);

    set_of_Is(i_g,1:min([length(y_bind),length(y_refc)])) = y_bind(1:min([length(y_bind),length(y_refc)]))./y_refc(1:min([length(y_bind),length(y_refc)]));



    subplot(1,2,i_fig), xlabel(x_label), ylabel(y_label), title(TheTitle)
    plot(d_centers,TheOutFun(set_of_Is(i_g,:)),'DisplayName',TheLegend,'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','w','Color',get_color((i_g)),'MarkerSize',8,'LineWidth',2), hold on
    set(gca,'fontsize',20), grid on, hold on
    clearvars index_in x_temp y_temp TheCode y_bind
    clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label Therloess blnFit y_ref



    % corrected DPF vs. d
    i_fig = 2; fun_x = @mean; fun_y = @sum; TheOutFun = @(x)(-log(x)); TheLineStyle = '-'; TheMarker = 'none'; x_label = 'd (cm)'; y_label = '\DeltaOD/d (cm^{-1})';
    subplot(1,2,i_fig), xlabel(x_label), ylabel(y_label), title(TheTitle)
    plot(d_centers,TheOutFun(set_of_Is(i_g,:))./d_centers,'DisplayName',TheLegend,'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','w','Color',get_color((i_g)),'MarkerSize',8,'LineWidth',2), hold on
    set(gca,'fontsize',20), grid on, hold on
    clearvars index_in x_temp y_temp TheCode TheTitle x_bind y_bind
    clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label Therloess blnFit y_ref

    clearvars t_db t_rf the_filename mua mus ind_diff ind_trns n_g TheColor ind_refc
    clearvars x_refc y_refc TheLegend a
end

subplot(1,2,1)
d_min = 1; d_max = 6; % max & min of signal
axis([d_min d_max 0 30]),
subplot(1,2,2)
axis([d_min d_max 0.1 10]),



% experimental data
set_of_Cs = [0  1   4   7];
for i_C = 3:length(set_of_Cs)
    subplot(1,2,1), plot(set_of_ds,-log(set_of_V(i_C,:)./set_of_V(2,:)),           'LineStyle','none','LineWidth',2,'DisplayName',['\DeltaC = ',num2str(set_of_Cs(i_C)-1),' g/L: data'],'Marker','square','MarkerSize',24,'MarkerFaceColor',get_color((i_C-2)),'MarkerEdgeColor','k'), hold on
    legend('show','Location','northwest','NumColumns',2), axis square
    subplot(1,2,2), plot(set_of_ds,-log(set_of_V(i_C,:)./set_of_V(2,:))./set_of_ds,'LineStyle','none','LineWidth',2,'DisplayName',['\DeltaC = ',num2str(set_of_Cs(i_C)-1),' g/L: data'],'Marker','square','MarkerSize',24,'MarkerFaceColor',get_color((i_C-2)),'MarkerEdgeColor','k'), hold on
    legend('show','Location','northeast','NumColumns',2), axis square

    [~,idx] = min(abs(pq_db.set_of_data(:,1)-set_of_g(i_C-2))); p_est = pq_db.set_of_data(idx,4); q_est = pq_db.set_of_data(idx,5); clearvars idx
    [~,idx] = min(abs(pq_db.set_of_data(:,1)-1.0            )); p_bas = pq_db.set_of_data(idx,4); q_bas = pq_db.set_of_data(idx,5); clearvars idx
    subplot(1,2,1), plot(set_of_ds,(p_est.*set_of_ds+q_est)-(p_bas.*set_of_ds+q_bas),'LineStyle','-.','LineWidth',2,'DisplayName',['\DeltaC = ',num2str(set_of_Cs(i_C)-1),' g/L: inv. dist.'],'Marker','none','Color',get_color((i_C-2))), hold on
    subplot(1,2,2), plot(set_of_ds,(p_est+q_est./set_of_ds)-(p_bas+q_bas./set_of_ds),'LineStyle','-.','LineWidth',2,'DisplayName',['\DeltaC = ',num2str(set_of_Cs(i_C)-1),' g/L: inv. dist.'],'Marker','none','Color',get_color((i_C-2))), hold on

    [~,idx] = min(abs(pq_db.set_of_data(:,1)-set_of_g(i_C-2))); c_est = pq_db.set_of_cnst_DPF(idx,1); clearvars idx
    [~,idx] = min(abs(pq_db.set_of_data(:,1)-1.0            )); c_bas = pq_db.set_of_cnst_DPF(idx,1); clearvars idx
    subplot(1,2,1), plot(set_of_ds,c_est.*ones(size(set_of_ds))-c_bas.*ones(size(set_of_ds)),'LineStyle',':','LineWidth',2,'DisplayName',['\DeltaC = ',num2str(set_of_Cs(i_C)-1),' g/L: constant'],'Marker','none','Color',get_color((i_C-2))), hold on
    subplot(1,2,2), plot(set_of_ds,c_est./set_of_ds - c_bas./set_of_ds,'LineStyle',':','LineWidth',2,'DisplayName',['\DeltaC = ',num2str(set_of_Cs(i_C)-1),' g/L: constant'],'Marker','none','Color',get_color((i_C-2))), hold on

    % [~,idx] = min(abs(pq_db.set_of_data(:,1)-set_of_g(i_C-2))); mua_est = pq_db.set_of_data(idx,2); mus_est = pq_db.set_of_data(idx,3); clearvars idx
    % [~,idx] = min(abs(pq_db.set_of_data(:,1)-1.0            )); mua_bas = pq_db.set_of_data(idx,2); mus_bas = pq_db.set_of_data(idx,3); clearvars idx
    % p_est = 0.493.*((mua_est).^(-0.507)).*(mus_est.^(+0.506)); q_est = 3.112.*((mua_est).^(-1.370)).*(mus_est.^(-0.198));
    % p_bas = 0.493.*((mua_bas).^(-0.507)).*(mus_bas.^(+0.506)); q_bas = 3.112.*((mua_bas).^(-1.370)).*(mus_bas.^(-0.198));
    % subplot(1,2,1), plot(set_of_ds,(p_est.*set_of_ds+q_est)-(p_bas.*set_of_ds+q_bas),'LineStyle','..','LineWidth',2,'DisplayName',['\DeltaC = ',num2str(set_of_Cs(i_C)-1),' g/L: inv. dist.'],'Marker','none','Color',get_color((i_C-2))), hold on
    % subplot(1,2,2), plot(set_of_ds,(p_est+q_est./set_of_ds)-(p_bas+q_bas./set_of_ds),'LineStyle','..','LineWidth',2,'DisplayName',['\DeltaC = ',num2str(set_of_Cs(i_C)-1),' g/L: inv. dist.'],'Marker','none','Color',get_color((i_C-2))), hold on

    clearvars p_est q_est p_bas q_bas c_bas c_est
    clearvars mua_est mus_est mua_bas mus_bas
end
clearvars i_C set_of_Cs i_g



clc
set_of_Cs = [0  1   4   7];
set_of_ind = nan(size(set_of_ds));
for i_d = 1:length(set_of_ds)
    [~,set_of_ind(i_d)] = min(abs(set_of_ds(i_d)-d_centers));
end
clearvars i_d
set_of_true = nan(2,length(set_of_ds));
set_of_estimate = nan(2,4,length(set_of_ds));
for i_C = 3:length(set_of_Cs)
    [~,i_ref] = min(abs(pq_db.set_of_data(:,1)-1              ));
    [~,idx  ] = min(abs(pq_db.set_of_data(:,1)-set_of_g(i_C-2))); mua_est = pq_db.set_of_data(idx,2); mus_est = pq_db.set_of_data(idx,3); clearvars idx

    [~,idx] = min(abs(pq_db.set_of_data(:,1)-set_of_g(i_C-2)));
    set_of_estimate(i_C-2,1,:) = -log(set_of_V(i_C,:)./set_of_V(2,:))./set_of_ds./(pq_db.set_of_cnst_DPF(idx,1) -            pq_db.set_of_cnst_DPF(i_ref,1));            % const
    set_of_estimate(i_C-2,2,:) = -log(set_of_V(i_C,:)./set_of_V(2,:))./set_of_ds./(pq_db.set_of_savg_DPF(idx,set_of_ind) -   pq_db.set_of_savg_DPF(i_ref,set_of_ind));   % savg
    set_of_estimate(i_C-2,3,:) = -log(set_of_V(i_C,:)./set_of_V(2,:))./set_of_ds./(pq_db.set_of_true_DPF(idx,set_of_ind) -   pq_db.set_of_true_DPF(i_ref,set_of_ind));   % true
    set_of_estimate(i_C-2,4,:) = -log(set_of_V(i_C,:)./set_of_V(2,:))./set_of_ds./(pq_db.set_of_idst_DPF(idx,2)./set_of_ds - pq_db.set_of_idst_DPF(i_ref,2)./set_of_ds); % inv. dist.

    set_of_true(i_C-2,:) = -log(set_of_Is(i_C-2,set_of_ind))./set_of_ds./(pq_db.set_of_true_DPF(idx,set_of_ind) - pq_db.set_of_true_DPF(i_ref,set_of_ind));
    clearvars idx i_ref

    disp(['DC = ',num2str(set_of_Cs(i_C)-1),' g/L'])
    disp(round(array2table(    squeeze(set_of_estimate(i_C-2,:,:))),2))
    disp(round(array2table(abs(squeeze(set_of_estimate(i_C-2,:,:))./mua_est-1)*100),1))
end

end

function [out] = get_color(index)
if     index==1
    out = [0.0000 0.4470 0.7410];
    out = 'b';
elseif index==2
    out = [0.8500 0.3250 0.0980];
    out = 'r';
elseif index==3
    out = [0.9290 0.6940 0.1250];set_of_colors = [4 7];

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
