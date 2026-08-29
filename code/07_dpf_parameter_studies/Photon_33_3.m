function [] = Photon_33_3 ()
% Repository group: 07_dpf_parameter_studies
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: w & s vs. d over mu_s & mu_a: more photons, more combinations
% the analysis: population stat - more consice compared to 33_2

clc
close all
set(0,'DefaultFigureWindowStyle','docked')

% optical & geometry properties
% z_air = 0.0; % the thickness of air layer
Lx = 29.1; Ly = 29.1; Lz = 05.0; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
% beam_phi = 0;     % polar angle of beam in cm
% beam_tht = 0;     % azimuthal angle of beam in cm
% beam_X = 0.0;     % X deviation of beam in cm
% beam_Y = 0.0;     % Y deviation of beam in cm

% n = 1.4; g = 0.95;
bln_diffuse = true;
set_of_muamus = [
    0.15 150
    0.15 350
    0.35 150
    0.35 350];



% arrays to save variables
N_bins = 250; the_percent = 5;
set_of_______OD = nan(length(set_of_muamus),2,N_bins);
set_of_conv_DPF = nan(length(set_of_muamus),3);
set_of_dfnd_DPF = nan(length(set_of_muamus),2,N_bins);
set_of_true_DPF = nan(length(set_of_muamus),2,N_bins);

for i_case = 1:size(set_of_muamus,1)
    % read dbase
    mua = set_of_muamus(i_case,1);
    mus = set_of_muamus(i_case,2);
    the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_33_mua_',sprintf('%.2f',mua),'_mus_',sprintf('%.2f',mus),'.mat'];
    t_db = load(the_filename); clearvars the_filename
    disp(['mua = ',num2str(mua,'%.2f'),', mus = ',num2str(mus,'%.0f')])

    % removing noisy points
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
    d_diff_edges = linspace(0, sqrt((Lx/2).^2+(Ly/2).^2      ),N_bins+1); [~,~,ind_diff] = histcounts(sqrt(t_db.x(t_db.c==0).^2+t_db.y(t_db.c==0).^2+t_db.z(t_db.c==0).^2),d_diff_edges); % d_diffuse bins
    d_trns_edges = linspace(Lz,sqrt((Lx/2).^2+(Ly/2).^2+Lz.^2),N_bins+1); [~,~,ind_trns] = histcounts(sqrt(t_db.x(t_db.c==1).^2+t_db.y(t_db.c==1).^2+t_db.z(t_db.c==1).^2),d_trns_edges); % d_Transmittance bins
    % figure(200+i_case)
    % subplot(1,2,1), histogram(sqrt(t_db.x(t_db.c==0).^2+t_db.y(t_db.c==0).^2+t_db.z(t_db.c==0).^2),N_bins,'FaceColor','r','EdgeColor','none'), axis tight, y_lim = ylim; hold on, plot([0,   0],[min(y_lim) max(y_lim)],'k-.'), plot([sqrt((Lx/2).^2+(Ly/2).^2),       sqrt((Lx/2).^2+(Ly/2).^2)],      [min(y_lim) max(y_lim)],'k-.'), hold off
    % subplot(1,2,2), histogram(sqrt(t_db.x(t_db.c==1).^2+t_db.y(t_db.c==1).^2+t_db.z(t_db.c==1).^2),N_bins,'FaceColor','r','EdgeColor','none'), axis tight, y_lim = ylim; hold on, plot([Lz, Lz],[min(y_lim) max(y_lim)],'k-.'), plot([sqrt((Lx/2).^2+(Ly/2).^2+Lz.^2), sqrt((Lx/2).^2+(Ly/2).^2+Lz.^2)],[min(y_lim) max(y_lim)],'k-.'), hold off
    clearvars d_trns_edges d_diff_edges



    % I vs. s (scatterplot)
    if bln_diffuse
        % (diffuse)
        index_in = ind_diff; TheCode = 0;
    else
        % (transmittance)
        index_in = ind_trns; TheCode = 1;
    end
    x_temp = t_db.s(t_db.c==TheCode);
    d_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2);
    x_bind = accumarray(index_in,d_temp,[],@(x)(mean(x)),nan);
    y_bind = accumarray(index_in,x_temp,[],@(x)(sum(exp(-mua.*x))./t_db.no_of_photons),nan);
    set_of_true_DPF(i_case,1,1:length(x_bind)) = x_bind;
    set_of_true_DPF(i_case,2,1:length(x_bind)) = -log(y_bind)./mua./x_bind;
    clearvars index_in TheCode x_temp d_temp x_bind y_bind



    % I vs. d
    i_fig = 1; fun_x = @mean; fun_y = @sum; TheOutFun = @(x)(-log(x)); TheLineStyle = 'none'; TheMarker = 'o'; x_label = 'd (cm)'; y_label = 'OD/\mu_a (cm^{-1})'; Therloess = 0.0; blnFit = false; blnShow = false;
    TheTitle = ['\mu_a = ',num2str(mua,'%.2f'),', \mu_s = ',num2str(mus,'%.0f')];
    if bln_diffuse
        % (diffuse)
        index_in = ind_diff; TheCode = 0; TheLegend = 'diffuse';
    else
        % (transmittance)
        index_in = ind_trns; TheCode = 1; TheLegend = 'transmittance';
    end
    x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.w(t_db.c==TheCode)./t_db.no_of_photons;
    [m,b,m_s,x_bind,y_bind] = make_a_subplot(i_fig,i_case,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,0.5,500,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit,TheLegend,'k',blnShow);
    set_of_conv_DPF(i_case,1) = m/mua;
    set_of_conv_DPF(i_case,2) = b/mua;
    set_of_conv_DPF(i_case,3) = m_s/mua;
    disp(['p = ',num2str(set_of_conv_DPF(i_case,1)),', q = ',num2str(set_of_conv_DPF(i_case,2))])
    disp(['-----------------------------------------------------------------------------------'])
    set_of_______OD(i_case,1,1:length(x_bind)) = x_bind;
    set_of_______OD(i_case,2,1:length(y_bind)) = TheOutFun(y_bind);
    figure(i_fig)
    subplot(2,2,i_case)
    plot(x_bind,TheOutFun(y_bind)/mua,'DisplayName','data', 'Marker','none','LineStyle','-','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6,'LineWidth',2.0,'Color','k'), hold on
    plot(x_bind,(x_bind*m_s)./mua,    'DisplayName','m*d',  'Marker','none','LineStyle','-','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',6,'LineWidth',2.0,'Color','r'), hold on
    plot(x_bind,(x_bind*m+b)./mua,    'DisplayName','p*d+q','Marker','none','LineStyle','-','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',6,'LineWidth',2.0,'Color','g'), hold on
    xlim([0  006]), xlabel(x_label), set(gca,'xtick',0:2:20)
    ylim([0  120]), ylabel(y_label), set(gca,'ytick',0:30:150)
    title(TheTitle), set(gca,'fontsize',14), grid on, hold on, legend('show','Location','southeast')
    clearvars index_in TheCode TheLegend x_temp y_temp x_bind y_bind m m_s b
    clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label Therloess blnFit blnShow



    % DPF vs. d
    i_fig = 2; fun_x = @mean; fun_y = @mean; TheOutFun = @(x)(x); TheLineStyle = '-'; TheMarker = 'none'; x_label = 'd (cm)'; y_label = 'DPF (a.u.)'; Therloess = 0.5; blnFit = false; blnShow = false;
    TheTitle = ['\mu_a = ',num2str(mua,'%.2f'),', \mu_s = ',num2str(mus,'%.0f')];
    if bln_diffuse
        % (diffuse)
        index_in = ind_diff; TheCode = 0; TheLegend = 'diffuse';
    else
        % (diffuse)
        index_in = ind_trns; TheCode = 1; TheLegend = 'transmittance';
    end
    x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.s(t_db.c==TheCode)./x_temp;
    [~,~,~,x_bind,y_bind] = make_a_subplot(i_fig,i_case,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,0.5,500,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit,TheLegend,'k',blnShow);
    set_of_dfnd_DPF(i_case,1,1:length(x_bind)) = x_bind; clearvars x_bind
    set_of_dfnd_DPF(i_case,2,1:length(y_bind)) = y_bind; clearvars y_bind
    figure(i_fig)
    subplot(2,2,i_case)
    plot(squeeze(set_of_dfnd_DPF(i_case,1,:)),squeeze(set_of_dfnd_DPF(i_case,2,:)),                                                     'DisplayName','<s>/d',   'Marker','none','LineStyle','-', 'MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',6,'LineWidth',2.0,'Color','b'), hold on
    plot(squeeze(set_of_dfnd_DPF(i_case,1,:)),ones(size(squeeze(set_of_dfnd_DPF(i_case,1,:)))).*set_of_conv_DPF(i_case,3),              'DisplayName','constant','Marker','none','LineStyle','-', 'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',6,'LineWidth',2.0,'Color','r'), hold on
    plot(squeeze(set_of_dfnd_DPF(i_case,1,:)),set_of_conv_DPF(i_case,1)+set_of_conv_DPF(i_case,2)./squeeze(set_of_dfnd_DPF(i_case,1,:)),'DisplayName','p+q/d',   'Marker','none','LineStyle','-', 'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',6,'LineWidth',2.0,'Color','g'), hold on
    plot(squeeze(set_of_true_DPF(i_case,1,:)),squeeze(set_of_true_DPF(i_case,2,:)),                                                     'DisplayName','true',    'Marker','none','LineStyle','-', 'MarkerEdgeColor','k','MarkerFaceColor','m','MarkerSize',6,'LineWidth',2.0,'Color','m'), hold on
    xlim([0  006]), xlabel(x_label), set(gca,'xtick',0:2:20)
    ylim([0  150]), ylabel(y_label), set(gca,'ytick',0:30:150)
    title(TheTitle), set(gca,'fontsize',14), grid on, hold on, legend('show','Location','northeast')
    clearvars index_in TheCode TheLegend x_temp y_temp x_bind y_bind m m_s b
    clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label Therloess blnFit blnShow



    % estimation of mu_a vs. d
    i_fig = 3; x_label = 'd (cm)'; y_label = '\mu_a (cm^{-1})';
    TheTitle = ['\mu_a = ',num2str(mua,'%.2f'),', \mu_s = ',num2str(mus,'%.0f')];
    if bln_diffuse
        % (diffuse)
        index_in = ind_diff; TheCode = 0; TheLegend = 'diffuse';
    else
        % (diffuse)
        index_in = ind_trns; TheCode = 1; TheLegend = 'transmittance';
    end
    figure(i_fig)
    subplot(2,2,i_case)
    plot(squeeze(set_of_______OD(i_case,1,:)),ones(size(squeeze(set_of_______OD(i_case,2,:)))).*mua,                                                                                                                  'DisplayName','true value',  'Marker','none','LineStyle','-', 'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6,'LineWidth',4.0,'Color','k'), hold on
    plot(squeeze(set_of_______OD(i_case,1,:)),squeeze(set_of_______OD(i_case,2,:))./squeeze(set_of_______OD(i_case,1,:))./(squeeze(set_of_dfnd_DPF(i_case,2,:))),                                                     'DisplayName','DPF = <s>/d', 'Marker','none','LineStyle','-', 'MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',6,'LineWidth',2.0,'Color','b'), hold on
    plot(squeeze(set_of_______OD(i_case,1,:)),squeeze(set_of_______OD(i_case,2,:))./squeeze(set_of_______OD(i_case,1,:))./(ones(size(squeeze(set_of_dfnd_DPF(i_case,1,:)))).*set_of_conv_DPF(i_case,3)),              'DisplayName','constant DPF','Marker','none','LineStyle','-', 'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',6,'LineWidth',2.0,'Color','r'), hold on
    plot(squeeze(set_of_______OD(i_case,1,:)),squeeze(set_of_______OD(i_case,2,:))./squeeze(set_of_______OD(i_case,1,:))./(set_of_conv_DPF(i_case,1)+set_of_conv_DPF(i_case,2)./squeeze(set_of_dfnd_DPF(i_case,1,:))),'DisplayName','DPF = p+q/d', 'Marker','none','LineStyle','-', 'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',6,'LineWidth',2.0,'Color','g'), hold on
    plot(squeeze(set_of_______OD(i_case,1,:)),squeeze(set_of_______OD(i_case,2,:))./squeeze(set_of_______OD(i_case,1,:))./(squeeze(set_of_true_DPF(i_case,2,:))),                                                     'DisplayName','true DPF',    'Marker','none','LineStyle','-', 'MarkerEdgeColor','k','MarkerFaceColor','m','MarkerSize',6,'LineWidth',2.0,'Color','m'), hold on

    xlim([0  006]), xlabel(x_label), set(gca,'xtick',0:2:20)
    ylim([mua-0.1  mua+0.1]), ylabel(y_label), set(gca,'ytick',0:0.1:2)
    title(TheTitle), set(gca,'fontsize',14), grid on, hold on, legend('show','Location','northeast')
    clearvars index_in TheCode TheLegend x_temp y_temp x_bind y_bind m m_s b
    clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label Therloess blnFit blnShow



    clearvars t_db the_filename mua mus ind_diff ind_trns
end
end

function [m,b,s_m,x_bind,y_bind] = make_a_subplot(i_fig,i_subplot,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,set_of_mua,set_of_mus,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit,TheLegend,TheEdgeColor,blnShow)
x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
y_bind = accumarray(index_in,y_temp,[],fun_y,nan);
TheColor = get_color(mua/max(set_of_mua),mus/max(set_of_mus));
mdl = fitlm(x_bind(~isnan(x_bind)),TheOutFun(y_bind(~isnan(x_bind))));
p = mdl.Coefficients.pValue(2);
if p<=0.05
    b = mdl.Coefficients.Estimate(1);
    m = mdl.Coefficients.Estimate(2);
    s_m = x_bind(~isnan(x_bind))\TheOutFun(y_bind(~isnan(x_bind)));
else
    b = nan;
    m = nan;
    s_m = nan;
end
if blnShow
    figure(i_fig)
    subplot(2,2,i_subplot)
    plot(x_bind,TheOutFun(y_bind),'Color',TheColor,'DisplayName',TheLegend,'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor',TheEdgeColor,'MarkerFaceColor',TheColor,'MarkerSize',6,'LineWidth',1.0), hold on
end
% if blnFit
%     plot(x_bind,m*x_bind+b,'Color',TheColor,'HandleVisibility','off'), hold on
% end
end
function [out] = get_color(mua_norm,mus_norm)
col_mua_0_mus_0 = [1 1 0];
col_mua_0_mus_1 = [1 0 0];
col_mua_1_mus_0 = [0 0 1];
col_mua_1_mus_1 = [1 0.0 1];
% out = [mua_norm 0 mus_norm];
out = ...
    (1-mua_norm)*(1-mus_norm)*col_mua_0_mus_0 + ...
    (  mua_norm)*(1-mus_norm)*col_mua_1_mus_0 + ...
    (1-mua_norm)*(  mus_norm)*col_mua_0_mus_1 + ...
    (  mua_norm)*(  mus_norm)*col_mua_1_mus_1;
end
