function [] = Photon_40_2 ()
% Repository group: 06_laboratory_phantom_models
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: Replication of Arnab's Experiment: Dissolving x grams of TiO₂ in a solution of 10 g Sodium Polyacrylate (PAS) in 440 mL of water.
% n_g = 2, 3 ,4, 7, ...
% general analysis: Arnab data

clc
close all

% size and code properties
set_of_g = [0,1,2,3];
N_bins = 250; d_ref = 2;
Lx = 29.1; Ly = 29.1; Lz = 10.0; % The size of the computational domain (in cm) to prevent reflections from the boundaries.

set_of_Cs = [0  1   4   7];
set_of_ds = [1	2	3	4	5	6	7	8];
set_of_V = [...
    2.73	1.76	    1.186666667	0.94	0.666666667	0.6	    0.56	0.52;
    nan     8.293333333	2.32    	0.8	    0.4	        0.28	0.12	0.12;
    nan     5.88	    0.96	    0.32	0.213333333	nan	    nan 	nan;
    nan     3.26	    0.72	    nan	    nan	        nan	    nan	    nan];

for i_g = 3:length(set_of_g)
    % read dbase
    the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_40_g_',sprintf('%.1f',set_of_g(i_g)),'.mat'];
    t_db = load(the_filename); clearvars the_filename
    % read refc
    the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_40_g_',sprintf('%.1f',set_of_g(2  )),'.mat'];
    t_rf = load(the_filename); clearvars the_filename
    n_g = set_of_g(i_g);
    TheLegend = ['g = ',num2str(n_g),': model'];
    disp(['g = ',num2str(n_g)])
    disp(['n = ',num2str(t_db.n)])

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
    clearvars d_trns_edges d_diff_edges






    % I vs. d
    i_fig = 1; fun_x = @mean; fun_y = @sum; TheOutFun = @(x)(-log(x)); TheLineStyle = '-'; TheMarker = 'none'; x_label = '\Deltad (cm)'; y_label = '\DeltaOD (a.u.)';
    index_in = ind_diff; TheCode = 0; TheTitle = '\deltaOD vs. \Deltad';
    x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
    y_temp = t_db.w(t_db.c==TheCode)./t_db.no_of_photons;                                            y_bind = accumarray(index_in,y_temp,[],fun_y,nan);
    [~,i_ref] = min(abs(x_bind-d_ref)); y_ref = y_bind(i_ref); clearvars i_ref
    y_bind = y_bind(3<=x_bind&x_bind<=8); x_bind = x_bind(3<=x_bind&x_bind<=8);

    figure(i_fig), xlabel(x_label), ylabel(y_label), title(TheTitle)
    plot(x_bind-d_ref,TheOutFun(y_bind./y_ref),'DisplayName',TheLegend,'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','w','Color',get_color((i_g)),'MarkerSize',8,'LineWidth',2), hold on
    set(gca,'fontsize',20), grid on, axis([0 7 0 4.5]), hold on
    clearvars index_in x_temp y_temp TheCode TheTitle x_bind y_bind
    clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label Therloess blnFit y_ref



    % corrected DPF vs. d
    i_fig = 2; fun_x = @mean; fun_y = @sum; TheOutFun = @(x)(-log(x)); TheLineStyle = '-'; TheMarker = 'none'; x_label = '\Deltad (cm)'; y_label = '\DeltaOD/d (cm^{-1})';
    index_in = ind_diff; TheCode = 0; TheTitle = '\deltaOD/d vs. \Deltad';
    x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
    y_temp = t_db.w(t_db.c==TheCode)./t_db.no_of_photons;                                            y_bind = accumarray(index_in,y_temp,[],fun_y,nan);
    [~,i_ref] = min(abs(x_bind-d_ref)); y_ref = y_bind(i_ref); clearvars i_ref
    y_bind = y_bind(3<=x_bind&x_bind<=8); x_bind = x_bind(3<=x_bind&x_bind<=8);

    figure(i_fig), xlabel(x_label), ylabel(y_label), title(TheTitle)
    plot(x_bind-d_ref,TheOutFun(y_bind./y_ref)./(x_bind-d_ref),'DisplayName',TheLegend,'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','w','Color',get_color((i_g)),'MarkerSize',8,'LineWidth',2), hold on
    set(gca,'fontsize',20), grid on, axis([0 7 0 2]), hold on
    clearvars index_in x_temp y_temp TheCode TheTitle x_bind y_bind
    clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label Therloess blnFit y_ref



    % corrected DPF vs. d
    i_fig = 3; fun_x = @mean; fun_y = @sum; TheOutFun = @(x)(-log(x)); TheLineStyle = '-'; TheMarker = 'none'; x_label = '\Deltad (cm)'; y_label = '\DeltaOD/d (cm^{-1})';
    index_in = ind_diff; TheCode = 0; TheTitle = '\deltaOD/d vs. \Deltad';
    x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
    y_temp = t_db.w(t_db.c==TheCode)./t_db.no_of_photons;                                            y_bind = accumarray(index_in,y_temp,[],fun_y,nan);
    y_bind = y_bind(1<=x_bind&x_bind<=5); x_bind = x_bind(1<=x_bind&x_bind<=5);

    index_in = ind_refc; TheCode = 0; TheTitle = '\deltaOD/d vs. \Deltad';
    x_temp = sqrt(t_rf.x(t_rf.c==TheCode).^2+t_rf.y(t_rf.c==TheCode).^2+t_rf.z(t_rf.c==TheCode).^2); x_refc = accumarray(index_in,x_temp,[],fun_x,nan);
    y_temp = t_rf.w(t_rf.c==TheCode)./t_rf.no_of_photons;                                            y_refc = accumarray(index_in,y_temp,[],fun_y,nan);
    y_refc = y_refc(1<=x_refc&x_refc<=5); x_refc = x_refc(1<=x_refc&x_refc<=5);

    figure(i_fig), xlabel(x_label), ylabel(y_label), title(TheTitle)
    plot(x_bind,TheOutFun(y_bind./y_refc)./x_bind,'DisplayName',TheLegend,'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','w','Color',get_color((i_g)),'MarkerSize',8,'LineWidth',2), hold on
    set(gca,'fontsize',20), grid on, axis([1 5 0 1]), hold on
    clearvars index_in x_temp y_temp TheCode TheTitle x_bind y_bind
    clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label Therloess blnFit y_ref



    clearvars t_db the_filename mua mus ind_diff ind_trns n_g TheColor
end

for i_C = 3:length(set_of_Cs)
    figure(1), plot(set_of_ds(3:end)-set_of_ds(2),-log(set_of_V(i_C,3:end)/set_of_V(i_C,2)),'LineStyle','none','LineWidth',2,'DisplayName',['g = ',num2str(set_of_Cs(i_C)),': data'],'Marker','square','MarkerSize',24,'MarkerFaceColor',get_color((i_C)),'MarkerEdgeColor','k'), hold on
    legend('show','Location','southeast'),
    figure(2), plot(set_of_ds(3:end)-set_of_ds(2),-log(set_of_V(i_C,3:end)/set_of_V(i_C,2))./(set_of_ds(3:end)-set_of_ds(2)),'LineStyle','none','LineWidth',2,'DisplayName',['g = ',num2str(set_of_Cs(i_C)),': data'],'Marker','square','MarkerSize',24,'MarkerFaceColor',get_color((i_C)),'MarkerEdgeColor','k'), hold on
    legend('show','Location','northeast'),
    figure(3), plot(set_of_ds,-log(set_of_V(i_C,:)./set_of_V(2,:))./set_of_ds,'LineStyle','none','LineWidth',2,'DisplayName',['C = ',num2str(set_of_Cs(i_C)),': data'],'Marker','square','MarkerSize',24,'MarkerFaceColor',get_color((i_C)),'MarkerEdgeColor','k'), hold on
    legend('show','Location','northeast'),
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
