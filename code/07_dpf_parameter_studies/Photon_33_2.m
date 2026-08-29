function [] = Photon_33_2 ()
% Repository group: 07_dpf_parameter_studies
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: w & s vs. d over mu_s & mu_a: more photons, more combinations
% the analysis: population stat - more consice compared to 33_1

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
set_of_muamus = [
    0.05 050
    0.05 250
    0.05 450
    0.25 050
    0.25 250
    0.25 450
    0.45 050
    0.45 250
    0.45 450];

% arrays to save variables
N_bins = 250; the_percent = 5;
set_of_______OD = nan(2,length(set_of_muamus),2,N_bins);
set_of_monv_DPF = nan(2,length(set_of_muamus),3);
set_of_conv_DPF = nan(2,length(set_of_muamus),3);
set_of_dfnd_DPF = nan(2,length(set_of_muamus),2,N_bins);
set_of_true_DPF = nan(2,length(set_of_muamus),2,N_bins);

for i_case = 1:length(set_of_muamus)
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
    clearvars d_trns_edges d_diff_edges



    % I vs. s (scatterplot)
    i_fig = 1; TheOutFun = @(x)(-log(x)); TheLineStyle = 'none'; TheMarker = 'o'; x_label = 's (cm)'; y_label = 'OD (a.u.)';
    TheTitle = ['\mu_a = ',num2str(mua,'%.2f'),', \mu_s = ',num2str(mus,'%.0f')];
    TheColor = get_color(mua/0.5,mus/500);
    % (diffuse)
    index_in = ind_diff; TheCode = 0; TheLegend = 'diffuse';
    x_temp = t_db.s(t_db.c==TheCode); y_temp = t_db.w(t_db.c==TheCode);
    %%
    d_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2);
    x_bind = accumarray(index_in,d_temp,[],@(x)(mean(x)),nan);
    y_bind = accumarray(index_in,x_temp,[],@(x)(sum(exp(-mua.*x))./t_db.no_of_photons),nan);
    set_of_true_DPF(1,i_case,1,1:length(x_bind)) = x_bind;
    set_of_true_DPF(1,i_case,2,1:length(x_bind)) = -log(y_bind)./mua./x_bind;
    clearvars d_temp x_bind y_bind
    %%
    x_temp = x_temp; y_temp = TheOutFun(y_temp);
    figure(i_fig), subplot(3,3,i_case)
    plot(x_temp,y_temp,'Color',TheColor,'DisplayName',TheLegend,'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor',TheColor,'MarkerSize',6), hold on
    xlim([0 100]), xlabel(x_label), set(gca,'xtick',0:50:150)
    ylim([0   5]), ylabel(y_label), set(gca,'ytick',0:1:5)
    title(TheTitle), set(gca,'fontsize',12), grid on, hold on, legend('show','Location','best')
    clearvars i_subplot index_in x_temp y_temp TheCode TheLegend mdl
    % (transmittance)
    index_in = ind_trns; TheCode = 1; TheLegend = 'transmittance';
    x_temp = t_db.s(t_db.c==TheCode); y_temp = t_db.w(t_db.c==TheCode);
    %%
    d_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2);
    x_bind = accumarray(index_in,d_temp,[],@(x)(mean(x)),nan);
    y_bind = accumarray(index_in,x_temp,[],@(x)(sum(exp(-mua.*x))./t_db.no_of_photons),nan);
    set_of_true_DPF(2,i_case,1,1:length(x_bind)) = x_bind;
    set_of_true_DPF(2,i_case,2,1:length(x_bind)) = -log(y_bind)./mua./x_bind;
    clearvars d_temp x_bind y_bind
    %%
    x_temp = x_temp; y_temp = TheOutFun(y_temp);
    figure(i_fig), subplot(3,3,i_case)
    plot(x_temp,y_temp,'Color',TheColor,'DisplayName',TheLegend,'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','w','MarkerFaceColor',TheColor,'MarkerSize',6), hold on
    xlim([0 100]), xlabel(x_label), set(gca,'xtick',0:50:150)
    ylim([0   5]), ylabel(y_label), set(gca,'ytick',0:1:5)
    title(TheTitle), set(gca,'fontsize',12), grid on, hold on, legend('show','Location','best')
    set(gca,'fontsize',12), grid on, hold on
    clearvars i_subplot index_in x_temp y_temp TheCode TheTitle mdl
    clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label TheLegend TheColor



    % I vs. d
    i_fig = 2; fun_x = @mean; fun_y = @sum; TheOutFun = @(x)(-log(x)); TheLineStyle = 'none'; TheMarker = 'o'; x_label = 'd (cm)'; y_label = 'OD (a.u.)'; Therloess = 0.0; blnFit = false; blnShow = true;
    TheTitle = ['\mu_a = ',num2str(mua,'%.2f'),', \mu_s = ',num2str(mus,'%.0f')];
    % (diffuse)
    index_in = ind_diff; TheCode = 0; TheLegend = 'diffuse';
    x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.w(t_db.c==TheCode)./t_db.no_of_photons;
    [m,b,m_s,x_bind,y_bind] = make_a_subplot(i_fig,i_case,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,0.5,500,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit,TheLegend,'k',blnShow);
    set_of_conv_DPF(1,i_case,1) = m/mua;
    set_of_conv_DPF(1,i_case,2) = b/mua;
    set_of_conv_DPF(1,i_case,3) = m_s/mua;
    set_of_______OD(1,i_case,1,1:length(x_bind)) = x_bind;
    set_of_______OD(1,i_case,2,1:length(y_bind)) = TheOutFun(y_bind);
    xlim([0  15]), xlabel(x_label), set(gca,'xtick',0:5:20)
    ylim([0  20]), ylabel(y_label), set(gca,'ytick',0:5:20)
    title(TheTitle), set(gca,'fontsize',12), grid on, hold on, legend('show','Location','best')
    figure(100+i_case), subplot(1,2,1), plot(x_bind,TheOutFun(y_bind),'DisplayName','data','linewidth',2), hold on, plot(x_bind,x_bind*m+b,'-','linewidth',2,'DisplayName','a*x+b'), plot(x_bind,x_bind*m_s,'-.','DisplayName','a*x','linewidth',2), hold off, axis tight, legend show, xlabel('d (cm)'), ylabel('OD/\mu'), set(gca,'fontsize',12), axis square
    clearvars x_bind y_bind m m_s b
    clearvars i_subplot index_in x_temp y_temp TheCode TheLegend
    % (transmittance)
    index_in = ind_trns; TheCode = 1; TheLegend = 'transmittance';
    x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.w(t_db.c==TheCode)./t_db.no_of_photons;
    [m,b,m_s,x_bind,y_bind] = make_a_subplot(i_fig,i_case,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,0.5,500,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit,TheLegend,'w',blnShow);
    set_of_conv_DPF(2,i_case,1) = m/mua; clearvars m
    set_of_conv_DPF(2,i_case,2) = b/mua; clearvars b
    set_of_conv_DPF(2,i_case,3) = m_s/mua; clearvars m_s
    set_of_______OD(2,i_case,1,1:length(x_bind)) = x_bind; clearvars x_bind
    set_of_______OD(2,i_case,2,1:length(y_bind)) = TheOutFun(y_bind); clearvars y_bind
    xlim([0  15]), xlabel(x_label), set(gca,'xtick',0:5:20)
    ylim([0  20]), ylabel(y_label), set(gca,'ytick',0:5:20)
    title(TheTitle), set(gca,'fontsize',12), grid on, hold on, legend('show','Location','best')
    clearvars i_subplot index_in x_temp y_temp TheCode TheTitle
    clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label Therloess blnFit TheLegend blnShow



    % min_{s} vs. d
    i_fig = 3; fun_x = @mean; fun_y = @(x)(prctile(x,the_percent)); TheOutFun = @(x)(x); TheLineStyle = '-'; TheMarker = '.'; x_label = 'd (cm)'; y_label = 's_{min} (cm)'; Therloess = 0.5; blnFit = false; blnShow = true;
    TheTitle = ['\mu_a = ',num2str(mua,'%.2f'),', \mu_s = ',num2str(mus,'%.0f')];
    % (diffuse)
    index_in = ind_diff; TheCode = 0; TheLegend = 'diffuse';
    x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.s(t_db.c==TheCode);
    [~,~,~,~,~] = make_a_subplot(i_fig,i_case,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,0.5,500,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit,TheLegend,'k',blnShow);
    xlim([0   20]), xlabel(x_label), set(gca,'xtick',0:5:20)
    ylim([0  100]), ylabel(y_label), set(gca,'ytick',0:20:100)
    title(TheTitle), set(gca,'fontsize',12), grid on, hold on, legend('show','Location','best')
    clearvars i_subplot index_in x_temp y_temp TheCode TheLegend
    % (transmittance)
    index_in = ind_trns; TheCode = 1; TheLegend = 'transmittance';
    x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.s(t_db.c==TheCode);
    [~,~,~,~,~] = make_a_subplot(i_fig,i_case,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,0.5,500,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit,TheLegend,'w',blnShow);
    xlim([0   20]), xlabel(x_label), set(gca,'xtick',0:5:20)
    ylim([0  100]), ylabel(y_label), set(gca,'ytick',0:20:100)
    title(TheTitle), set(gca,'fontsize',12), grid on, hold on, legend('show','Location','best')
    clearvars i_subplot index_in x_temp y_temp TheCode TheTitle
    clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label Therloess blnFit TheLegend blnShow

    % max_{s} vs. d
    i_fig = 3; fun_x = @mean; fun_y = @(x)(prctile(x,100-the_percent)); TheOutFun = @(x)(x); TheLineStyle = '-'; TheMarker = '.'; x_label = 'd (cm)'; y_label = 's_{max} (cm)'; Therloess = 0.5; blnFit = false; blnShow = true;
    TheTitle = ['\mu_a = ',num2str(mua,'%.2f'),', \mu_s = ',num2str(mus,'%.0f')];
    % (diffuse)
    index_in = ind_diff; TheCode = 0; TheLegend = 'diffuse';
    x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.s(t_db.c==TheCode);
    [~,~,~,~,~] = make_a_subplot(i_fig,i_case,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,0.5,500,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit,TheLegend,'k',blnShow);
    xlim([0   20]), xlabel(x_label), set(gca,'xtick',0:5:20)
    ylim([0  100]), ylabel(y_label), set(gca,'ytick',0:20:100)
    title(TheTitle), set(gca,'fontsize',12), grid on, hold on, legend('show','Location','best')
    clearvars i_subplot index_in x_temp y_temp TheCode TheLegend
    % (transmittance)
    index_in = ind_trns; TheCode = 1; TheLegend = 'transmittance';
    x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.s(t_db.c==TheCode);
    [~,~,~,~,~] = make_a_subplot(i_fig,i_case,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,0.5,500,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit,TheLegend,'w',blnShow);
    xlim([0   20]), xlabel(x_label), set(gca,'xtick',0:5:20)
    ylim([0  100]), ylabel(y_label), set(gca,'ytick',0:20:100)
    title(TheTitle), set(gca,'fontsize',12), grid on, hold on, legend('show','Location','best')
    clearvars i_subplot index_in x_temp y_temp TheCode TheTitle
    clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label Therloess blnFit TheLegend blnShow



    % s vs. d
    i_fig = 4; fun_x = @mean; fun_y = @mean; TheOutFun = @(x)(x); TheLineStyle = '-'; TheMarker = '.'; x_label = 'd (cm)'; y_label = 's (cm)'; Therloess = 0.5; blnFit = false; blnShow = true;
    TheTitle = ['\mu_a = ',num2str(mua,'%.2f'),', \mu_s = ',num2str(mus,'%.0f')];
    % (diffuse)
    index_in = ind_diff; TheCode = 0; TheLegend = 'diffuse';
    x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.s(t_db.c==TheCode);
    [m,b,m_s,~,~] = make_a_subplot(i_fig,i_case,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,0.5,500,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit,TheLegend,'k',blnShow);
    set_of_monv_DPF(1,i_case,1) = m; clearvars m
    set_of_monv_DPF(1,i_case,2) = b; clearvars b
    set_of_monv_DPF(1,i_case,3) = m_s; clearvars m_s
    xlim([0   20]), xlabel(x_label), set(gca,'xtick',0:5:20)
    ylim([0  100]), ylabel(y_label), set(gca,'ytick',0:20:100)
    title(TheTitle), set(gca,'fontsize',12), grid on, hold on, legend('show','Location','best')
    clearvars i_subplot index_in x_temp y_temp TheCode TheLegend
    % (transmittance)
    index_in = ind_trns; TheCode = 1; TheLegend = 'transmittance';
    x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.s(t_db.c==TheCode);
    [m,b,m_s,~,~] = make_a_subplot(i_fig,i_case,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,0.5,500,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit,TheLegend,'w',blnShow);
    set_of_monv_DPF(2,i_case,1) = m; clearvars m
    set_of_monv_DPF(2,i_case,2) = b; clearvars b
    set_of_monv_DPF(2,i_case,3) = m_s; clearvars m_s
    xlim([0   20]), xlabel(x_label), set(gca,'xtick',0:5:20)
    ylim([0  100]), ylabel(y_label), set(gca,'ytick',0:20:100)
    title(TheTitle), set(gca,'fontsize',12), grid on, hold on, legend('show','Location','best')
    clearvars i_subplot index_in x_temp y_temp TheCode TheTitle
    clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label Therloess blnFit TheLegend blnShow



    % DPF vs. d
    i_fig = 5; fun_x = @mean; fun_y = @mean; TheOutFun = @(x)(x); TheLineStyle = '-'; TheMarker = '.'; x_label = 'd (cm)'; y_label = '<s>/d (a.u.)'; Therloess = 0.5; blnFit = false; blnShow = true;
    TheTitle = ['\mu_a = ',num2str(mua,'%.2f'),', \mu_s = ',num2str(mus,'%.0f')];
    % (diffuse)
    index_in = ind_diff; TheCode = 0; TheLegend = 'diffuse';
    x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.s(t_db.c==TheCode)./x_temp;
    [~,~,~,x_bind,y_bind] = make_a_subplot(i_fig,i_case,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,0.5,500,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit,TheLegend,'k',blnShow);
    set_of_dfnd_DPF(1,i_case,1,1:length(x_bind)) = x_bind; clearvars x_bind
    set_of_dfnd_DPF(1,i_case,2,1:length(y_bind)) = y_bind; clearvars y_bind
    xlim([0   10]), xlabel(x_label), set(gca,'xtick',0:5:20)
    ylim([0   20]), ylabel(y_label), set(gca,'ytick',0:5:100)
    title(TheTitle), set(gca,'fontsize',12), grid on, hold on, legend('show','Location','best')
    clearvars i_subplot index_in x_temp y_temp TheCode TheLegend
    % (transmittance)
    index_in = ind_trns; TheCode = 1; TheLegend = 'transmittance';
    x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.s(t_db.c==TheCode)./x_temp;
    [~,~,~,x_bind,y_bind] = make_a_subplot(i_fig,i_case,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,0.5,500,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit,TheLegend,'w',blnShow);
    set_of_dfnd_DPF(2,i_case,1,1:length(x_bind)) = x_bind; clearvars x_bind
    set_of_dfnd_DPF(2,i_case,2,1:length(y_bind)) = y_bind; clearvars y_bind
    xlim([0   10]), xlabel(x_label), set(gca,'xtick',0:5:20)
    ylim([0   25]), ylabel(y_label), set(gca,'ytick',0:5:100)
    title(TheTitle), set(gca,'fontsize',12), grid on, hold on, legend('show','Location','best')
    clearvars i_subplot index_in x_temp y_temp TheCode TheTitle
    clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label Therloess blnFit TheLegend blnShow



    % histogram of s
    i_fig = 6;
    TheColor = get_color(mua/0.5,mus/500);
    % (diffuse)
    TheCode = 0; TheTitle = ['\mu_a = ',num2str(mua,'%.2f'),', \mu_s = ',num2str(mus,'%.0f')];
    y_temp = t_db.s(t_db.c==TheCode);
    [n_temp,x_temp] = histcounts(y_temp,'BinEdges',linspace(0,150)); x_temp = 1/2*(x_temp(1:end-1)+x_temp(2:end-0)); n_temp = n_temp./t_db.no_of_photons*100;
    if ~isempty(y_temp)>0
        figure(i_fig), subplot(3,3,i_case), semilogy(x_temp,n_temp,'Color',TheColor,'LineWidth',1,'DisplayName','diffuse','Marker','.','MarkerEdgeColor','k'), hold on
        xlim([0  100]), xlabel('s (cm)'),   set(gca,'xtick',0:25:100)
        ylim([0  100]), ylabel('freq (%)'), set(gca,'ytick',0:25:100)
        title(TheTitle), set(gca,'fontsize',12), grid on, hold on, legend('show','Location','best')
        set(gca,'fontsize',12), grid on, hold on
    end
    clearvars i_subplot index_in x_temp y_temp n_temp TheCode
    % (transmittance)
    TheCode = 1;
    y_temp = t_db.s(t_db.c==TheCode);
    [n_temp,x_temp] = histcounts(y_temp,'BinEdges',linspace(0,150)); x_temp = 1/2*(x_temp(1:end-1)+x_temp(2:end-0)); n_temp = n_temp./t_db.no_of_photons*100;
    if ~isempty(y_temp)>0
        figure(i_fig), subplot(3,3,i_case), semilogy(x_temp,n_temp,'Color',TheColor,'LineWidth',1,'DisplayName','transmittance','Marker','.','MarkerEdgeColor','w'), hold on
        xlabel('s (cm)'), ylabel('freq (%)'), title(TheTitle), axis([0 100 0 100])
        set(gca,'fontsize',12), grid on, hold on
    end
    clearvars i_subplot index_in x_temp y_temp n_temp TheCode TheTitle TheColor
    clearvars i_fig



    % N vs. d
    i_fig = 7; fun_x = @mean; fun_y = @sum; TheOutFun = @(x)(x); TheLineStyle = '-'; TheMarker = '.'; x_label = 'd (cm)'; y_label = 'N (#)'; Therloess = 0.0; blnFit = false; blnShow = true;
    TheTitle = ['\mu_a = ',num2str(mua,'%.2f'),', \mu_s = ',num2str(mus,'%.0f')];
    % (diffuse)
    index_in = ind_diff; TheCode = 0; TheLegend = 'diffuse';
    x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = ones(size(t_db.s(t_db.c==TheCode)))./t_db.no_of_photons*100;
    [~,~,~,~,~] = make_a_subplot(i_fig,i_case,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,0.5,500,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit,TheLegend,'k',blnShow);
    xlim([0   20]), xlabel(x_label), set(gca,'xtick',0:5:20)
    ylim([0   5e1]), ylabel(y_label), yscale("log")
    title(TheTitle), set(gca,'fontsize',12), grid on, hold on, legend('show','Location','best')
    clearvars i_subplot index_in x_temp y_temp TheCode TheLegend
    % (transmittance)
    index_in = ind_trns; TheCode = 1; TheLegend = 'transmittance';
    x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = ones(size(t_db.s(t_db.c==TheCode)))./t_db.no_of_photons*100;
    [~,~,~,~,~] = make_a_subplot(i_fig,i_case,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,0.5,500,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit,TheLegend,'w',blnShow);
    xlim([0   20]), xlabel(x_label), set(gca,'xtick',0:5:20)
    ylim([0   5e1]), ylabel(y_label), yscale("log")
    title(TheTitle), set(gca,'fontsize',12), grid on, hold on, legend('show','Location','best')
    clearvars i_subplot index_in x_temp y_temp TheCode TheTitle
    clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label Therloess blnFit TheLegend blnShow



    % w_bar vs. d
    i_fig = 8; fun_x = @mean; fun_y = @mean; TheOutFun = @(x)(-log(x)); TheLineStyle = '-'; TheMarker = '.'; x_label = 'd (cm)'; y_label = 'w_{bar} (a.u.)'; Therloess = 0.5; blnFit = false; blnShow = true;
    TheTitle = ['\mu_a = ',num2str(mua,'%.2f'),', \mu_s = ',num2str(mus,'%.0f')];
    % (diffuse)
    index_in = ind_diff; TheCode = 0; TheLegend = 'diffuse';
    x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.w(t_db.c==TheCode);
    [~,~,~,~,~] = make_a_subplot(i_fig,i_case,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,0.5,500,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit,TheLegend,'k',blnShow);
    xlim([0   10]), xlabel(x_label), set(gca,'xtick',0:5:20)
    ylim([0    5]), ylabel(y_label), set(gca,'ytick',0:1:5)
    title(TheTitle), set(gca,'fontsize',12), grid on, hold on, legend('show','Location','best')
    clearvars i_subplot index_in x_temp y_temp TheCode TheLegend
    % (transmittance)
    index_in = ind_trns; TheCode = 1; TheLegend = 'transmittance';
    x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.w(t_db.c==TheCode);
    [~,~,~,~,~] = make_a_subplot(i_fig,i_case,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,0.5,500,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit,TheLegend,'w',blnShow);
    xlim([0   10]), xlabel(x_label), set(gca,'xtick',0:5:20)
    ylim([0    5]), ylabel(y_label), set(gca,'ytick',0:1:50)
    title(TheTitle), set(gca,'fontsize',12), grid on, hold on, legend('show','Location','best')
    clearvars i_subplot index_in x_temp y_temp TheCode TheTitle
    clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label Therloess blnFit TheLegend blnShow






    % delta(I)/delta(d) vs. d
    i_fig = 9; fun_x = @mean; fun_y = @sum; TheOutFun = @(x)(-log(x)); TheLineStyle = '-'; TheMarker = 'none'; x_label = 'd (cm)'; y_label = 'OD/d (cm^{-1})'; Therloess = 0.0; blnFit = false; blnShow = false;
    TheTitle = ['\mu_a = ',num2str(mua,'%.2f'),', \mu_s = ',num2str(mus,'%.0f')];
    TheColor = get_color(mua/0.5,mus/500);
    % (diffuse)
    figure(i_fig), subplot(3,3,i_case)
    index_in = ind_diff; TheCode = 0; TheLegend = 'diffuse';
    x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.w(t_db.c==TheCode)./t_db.no_of_photons;
    [~,~,~,x_bind,y_bind] = make_a_subplot(i_fig,i_case,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,0.5,500,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit,TheLegend,'k',blnShow);
    plot(x_bind,TheOutFun(y_bind)./x_bind,'Color',TheColor,'DisplayName',TheLegend,'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor',TheColor,'MarkerSize',6,'LineWidth',2.0), hold on
    xlim([0  10]), xlabel(x_label), set(gca,'xtick',0:5:20)
    ylim([0  20]), ylabel(y_label), set(gca,'ytick',0:5:20)
    clearvars i_subplot index_in x_temp y_temp TheCode TheLegend x_bind y_bind
    % (transmittance)
    figure(i_fig), subplot(3,3,i_case)
    index_in = ind_trns; TheCode = 1; TheLegend = 'transmittance'; TheLineStyle = '-.';
    x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.w(t_db.c==TheCode)./t_db.no_of_photons;
    [~,~,~,x_bind,y_bind] = make_a_subplot(i_fig,i_case,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,0.5,500,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit,TheLegend,'k',blnShow);
    plot(x_bind,TheOutFun(y_bind)./x_bind,'Color',TheColor,'DisplayName',TheLegend,'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','w','MarkerFaceColor',TheColor,'MarkerSize',6,'LineWidth',2.0), hold on
    xlim([0  10]), xlabel(x_label), set(gca,'xtick',0:5:20)
    ylim([0  20]), ylabel(y_label), set(gca,'ytick',0:5:20)
    title(TheTitle), set(gca,'fontsize',12), grid on, hold on, legend('show','Location','best')
    clearvars i_subplot index_in x_temp y_temp TheCode TheLegend x_bind y_bind
    clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label Therloess blnFit TheLegend blnShow TheColor TheTitle



    % DPF vs. d
    TheTitle = ['\mu_a = ',num2str(mua,'%.2f'),', \mu_s = ',num2str(mus,'%.0f')]; TheColor = get_color(mua/0.5,mus/500); x_label = 'd (cm)'; y_label = 'DPF (a.u.)';
    % (diffuse)
    i_fig = 10; the_propagation_index = 1; figure(i_fig), subplot(3,3,i_case); ThePreLegend = 'diffuse';
    TheLegend = 'conv. DPF'; TheLineStyle = ':'; TheMarker = 'none';
    x_temp = squeeze(set_of_dfnd_DPF(the_propagation_index,i_case,1,:));
    y_temp = ones(size(squeeze(set_of_dfnd_DPF(the_propagation_index,i_case,1,:)))).*squeeze(set_of_conv_DPF(the_propagation_index,i_case,3));
    semilogy(x_temp,y_temp,'Color',TheColor,'DisplayName',[ThePreLegend ,': ',TheLegend],'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor',TheColor,'MarkerSize',6,'LineWidth',2), hold on
    clearvars TheLegend TheLineStyle TheMarker x_temp y_temp
    TheLegend = 'Imp. conv. DPF'; TheLineStyle = '-.'; TheMarker = 'none';
    x_temp = squeeze(set_of_dfnd_DPF(the_propagation_index,i_case,1,:));
    y_temp = squeeze(set_of_conv_DPF(the_propagation_index,i_case,1))+squeeze(set_of_conv_DPF(the_propagation_index,i_case,2))./squeeze(set_of_dfnd_DPF(the_propagation_index,i_case,1,:));
    semilogy(x_temp,y_temp,'Color',TheColor,'DisplayName',[ThePreLegend ,': ',TheLegend],'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor',TheColor,'MarkerSize',6,'LineWidth',2), hold on
    clearvars TheLegend TheLineStyle TheMarker x_temp y_temp
    TheLegend = '<s>/d'; TheLineStyle = '--'; TheMarker = 'none';
    x_temp = squeeze(set_of_dfnd_DPF(the_propagation_index,i_case,1,:));
    y_temp = squeeze(set_of_dfnd_DPF(the_propagation_index,i_case,2,:));
    semilogy(x_temp,y_temp,'Color',TheColor,'DisplayName',[ThePreLegend ,': ',TheLegend],'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor',TheColor,'MarkerSize',6,'LineWidth',2), hold on
    clearvars TheLegend TheLineStyle TheMarker x_temp y_temp
    TheLegend = 'true DPF'; TheLineStyle = '-'; TheMarker = 'none';
    x_temp = squeeze(set_of_true_DPF(the_propagation_index,i_case,1,:));
    y_temp = squeeze(set_of_true_DPF(the_propagation_index,i_case,2,:));
    semilogy(x_temp,y_temp,'Color',TheColor,'DisplayName',[ThePreLegend ,': ',TheLegend],'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor',TheColor,'MarkerSize',6,'LineWidth',2), hold on
    clearvars TheLegend TheLineStyle TheMarker x_temp y_temp
    clearvars the_propagation_index ThePreLegend
    % (transmittance)
    i_fig = 10; the_propagation_index = 2; figure(i_fig), subplot(3,3,i_case); ThePreLegend = 'transmittance';
    TheLegend = 'conv. DPF'; TheLineStyle = ':'; TheMarker = '.';
    x_temp = squeeze(set_of_dfnd_DPF(the_propagation_index,i_case,1,:));
    y_temp = ones(size(squeeze(set_of_dfnd_DPF(the_propagation_index,i_case,1,:)))).*squeeze(set_of_conv_DPF(the_propagation_index,i_case,3));
    semilogy(x_temp,y_temp,'Color',TheColor,'DisplayName',[ThePreLegend ,': ',TheLegend],'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6,'LineWidth',2), hold on
    clearvars TheLegend TheLineStyle TheMarker x_temp y_temp
    TheLegend = 'Imp. conv. DPF'; TheLineStyle = '-.'; TheMarker = '.';
    x_temp = squeeze(set_of_dfnd_DPF(the_propagation_index,i_case,1,:));
    y_temp = squeeze(set_of_conv_DPF(the_propagation_index,i_case,1))+squeeze(set_of_conv_DPF(the_propagation_index,i_case,2))./squeeze(set_of_dfnd_DPF(the_propagation_index,i_case,1,:));
    semilogy(x_temp,y_temp,'Color',TheColor,'DisplayName',[ThePreLegend ,': ',TheLegend],'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6,'LineWidth',2), hold on
    clearvars TheLegend TheLineStyle TheMarker x_temp y_temp
    TheLegend = '<s>/d'; TheLineStyle = '--'; TheMarker = '.';
    x_temp = squeeze(set_of_dfnd_DPF(the_propagation_index,i_case,1,:));
    y_temp = squeeze(set_of_dfnd_DPF(the_propagation_index,i_case,2,:));
    semilogy(x_temp,y_temp,'Color',TheColor,'DisplayName',[ThePreLegend ,': ',TheLegend],'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6,'LineWidth',2), hold on
    clearvars TheLegend TheLineStyle TheMarker x_temp y_temp
    TheLegend = 'true DPF'; TheLineStyle = '-'; TheMarker = '.';
    x_temp = squeeze(set_of_true_DPF(the_propagation_index,i_case,1,:));
    y_temp = squeeze(set_of_true_DPF(the_propagation_index,i_case,2,:));
    semilogy(x_temp,y_temp,'Color',TheColor,'DisplayName',[ThePreLegend ,': ',TheLegend],'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6,'LineWidth',2), hold on
    clearvars TheLegend TheLineStyle TheMarker x_temp y_temp
    xlim([0    10]), xlabel(x_label), set(gca,'xtick',0:5:10)
    ylim([0  1000]), ylabel(y_label), set(gca,'ytick',[0,1,10,100,1000])
    title(TheTitle), set(gca,'fontsize',12), grid on, hold on, % legend('show','Location','best')
    clearvars the_propagation_index ThePreLegend
    clearvars i_fig TheColor TheTitle x_label y_label



    % estimation of mu_a vs. d
    TheTitle = ['\mu_a = ',num2str(mua,'%.2f'),', \mu_s = ',num2str(mus,'%.0f')]; TheColor = get_color(mua/0.5,mus/500); x_label = 'd (cm)'; y_label = '\mu_a (cm^{-1})';
    i_fig = 11; the_propagation_index = 1; figure(i_fig), subplot(3,3,i_case);
    TheLegend = 'true \mu_a'; TheLineStyle = 'none'; TheMarker = 'o';
    d_temp = squeeze(set_of_______OD(the_propagation_index,i_case,1,:));
    m_temp = ones(size(squeeze(set_of_______OD(the_propagation_index,i_case,2,:)))).*mua;
    plot(d_temp,m_temp,'Color','k','DisplayName',TheLegend,'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor',TheColor,'MarkerSize',6), hold on
    clearvars TheLegend TheLineStyle TheMarker d_temp m_temp the_propagation_index i_fig
    % (diffuse)
    i_fig = 11; the_propagation_index = 1; figure(i_fig), subplot(3,3,i_case); ThePreLegend = 'diffuse';
    TheLegend = 'conv. \mu_a'; TheLineStyle = ':'; TheMarker = 'none';
    d_temp = squeeze(set_of_______OD(the_propagation_index,i_case,1,:));
    m_temp = squeeze(set_of_______OD(the_propagation_index,i_case,2,:))./squeeze(set_of_______OD(the_propagation_index,i_case,1,:))./(set_of_conv_DPF(the_propagation_index,i_case,3)*squeeze(set_of_______OD(the_propagation_index,i_case,1,:)));
    plot(d_temp,m_temp,'Color',TheColor,'DisplayName',[ThePreLegend ,': ',TheLegend],'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6,'LineWidth',2), hold on
    clearvars TheLegend TheLineStyle TheMarker d_temp m_temp
    TheLegend = 'Imp. conv. \mu_a'; TheLineStyle = '-.'; TheMarker = 'none';
    d_temp = squeeze(set_of_______OD(the_propagation_index,i_case,1,:));
    m_temp = squeeze(set_of_______OD(the_propagation_index,i_case,2,:))./(set_of_conv_DPF(the_propagation_index,i_case,1)*squeeze(set_of_______OD(the_propagation_index,i_case,1,:))+set_of_conv_DPF(the_propagation_index,i_case,2));
    plot(d_temp,m_temp,'Color',TheColor,'DisplayName',[ThePreLegend ,': ',TheLegend],'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6,'LineWidth',2), hold on
    clearvars TheLegend TheLineStyle TheMarker d_temp m_temp
    TheLegend = '<s>/d \mu_a'; TheLineStyle = '--'; TheMarker = 'none';
    d_temp = squeeze(set_of_______OD(the_propagation_index,i_case,1,:));
    m_temp = squeeze(set_of_______OD(the_propagation_index,i_case,2,:))./squeeze(set_of_______OD(the_propagation_index,i_case,1,:))./squeeze(set_of_dfnd_DPF(the_propagation_index,i_case,2,:));
    plot(d_temp,m_temp,'Color',TheColor,'DisplayName',[ThePreLegend ,': ',TheLegend],'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6,'LineWidth',2), hold on
    clearvars TheLegend TheLineStyle TheMarker d_temp m_temp
    TheLegend = 'true \mu_a'; TheLineStyle = '-'; TheMarker = 'none';
    d_temp = squeeze(set_of_______OD(the_propagation_index,i_case,1,:));
    m_temp = squeeze(set_of_______OD(the_propagation_index,i_case,2,:))./squeeze(set_of_______OD(the_propagation_index,i_case,1,:))./squeeze(set_of_true_DPF(the_propagation_index,i_case,2,:));
    plot(d_temp,m_temp,'Color',TheColor,'DisplayName',[ThePreLegend ,': ',TheLegend],'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6,'LineWidth',2), hold on
    clearvars TheLegend TheLineStyle TheMarker d_temp m_temp
    clearvars the_propagation_index ThePreLegend
    % (transmittance)
    i_fig = 11; the_propagation_index = 2; figure(i_fig), subplot(3,3,i_case); ThePreLegend = 'transmittance';
    TheLegend = 'conv. \mu_a'; TheLineStyle = ':'; TheMarker = '.';
    d_temp = squeeze(set_of_______OD(the_propagation_index,i_case,1,:));
    m_temp = squeeze(set_of_______OD(the_propagation_index,i_case,2,:))./squeeze(set_of_______OD(the_propagation_index,i_case,1,:))./(set_of_conv_DPF(the_propagation_index,i_case,3)*squeeze(set_of_______OD(the_propagation_index,i_case,1,:)));
    plot(d_temp,m_temp,'Color',TheColor,'DisplayName',[ThePreLegend ,': ',TheLegend],'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6,'LineWidth',2), hold on
    clearvars TheLegend TheLineStyle TheMarker d_temp m_temp
    TheLegend = 'Imp. conv. \mu_a'; TheLineStyle = '-.'; TheMarker = '.';
    d_temp = squeeze(set_of_______OD(the_propagation_index,i_case,1,:));
    m_temp = squeeze(set_of_______OD(the_propagation_index,i_case,2,:))./(set_of_conv_DPF(the_propagation_index,i_case,1)*squeeze(set_of_______OD(the_propagation_index,i_case,1,:))+set_of_conv_DPF(the_propagation_index,i_case,2));
    plot(d_temp,m_temp,'Color',TheColor,'DisplayName',[ThePreLegend ,': ',TheLegend],'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6,'LineWidth',2), hold on
    clearvars TheLegend TheLineStyle TheMarker d_temp m_temp
    TheLegend = '<s>/d \mu_a'; TheLineStyle = '--'; TheMarker = '.';
    d_temp = squeeze(set_of_______OD(the_propagation_index,i_case,1,:));
    m_temp = squeeze(set_of_______OD(the_propagation_index,i_case,2,:))./squeeze(set_of_______OD(the_propagation_index,i_case,1,:))./squeeze(set_of_dfnd_DPF(the_propagation_index,i_case,2,:));
    plot(d_temp,m_temp,'Color',TheColor,'DisplayName',[ThePreLegend ,': ',TheLegend],'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6,'LineWidth',2), hold on
    clearvars TheLegend TheLineStyle TheMarker d_temp m_temp
    TheLegend = 'true \mu_a'; TheLineStyle = '-'; TheMarker = '.';
    d_temp = squeeze(set_of_______OD(the_propagation_index,i_case,1,:));
    m_temp = squeeze(set_of_______OD(the_propagation_index,i_case,2,:))./squeeze(set_of_______OD(the_propagation_index,i_case,1,:))./squeeze(set_of_true_DPF(the_propagation_index,i_case,2,:));
    plot(d_temp,m_temp,'Color',TheColor,'DisplayName',[ThePreLegend ,': ',TheLegend],'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6,'LineWidth',2), hold on
    clearvars TheLegend TheLineStyle TheMarker d_temp m_temp
    xlim([0    10]), xlabel(x_label), set(gca,'xtick',0:1:10)
    ylim([mua-0.1   mua+0.1]), ylabel(y_label), set(gca,'ytick',[0:0.1:1.0])
    title(TheTitle), set(gca,'fontsize',12), grid on, hold on, % legend('show','Location','best')
    clearvars the_propagation_index ThePreLegend
    clearvars i_fig TheColor TheTitle x_label y_label



    % error in estimation of mu_a vs. d
    TheTitle = ['\mu_a = ',num2str(mua,'%.2f'),', \mu_s = ',num2str(mus,'%.0f')]; TheColor = get_color(mua/0.5,mus/500); x_label = 'd (cm)'; y_label = 'error in \mu_a estimation';
    % (diffuse)
    i_fig = 12; the_propagation_index = 1; figure(i_fig), subplot(3,3,i_case); ThePreLegend = 'diffuse';
    TheLegend = 'conv. \mu_a'; TheLineStyle = ':'; TheMarker = 'none';
    d_temp = squeeze(set_of_______OD(the_propagation_index,i_case,1,:));
    m_temp = squeeze(set_of_______OD(the_propagation_index,i_case,2,:))./squeeze(set_of_______OD(the_propagation_index,i_case,1,:))./set_of_conv_DPF(the_propagation_index,i_case,3);
    m_temp = (m_temp-mua)./mua*100;
    plot(d_temp,m_temp,'Color',TheColor,'DisplayName',[ThePreLegend ,': ',TheLegend],'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6,'LineWidth',2), hold on
    figure(100+i_case), subplot(1,2,2), plot(d_temp,m_temp,'Color',TheColor,'DisplayName',[ThePreLegend ,': ',TheLegend],'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6,'LineWidth',2), hold on, legend show, axis square, axis tight, ylim([-50 250]), xlabel('d (cm)'), ylabel('est. error (%)'), set(gca,'fontsize',12), axis square
    figure(i_fig), subplot(3,3,i_case);
    TheLegend = 'Imp. conv. \mu_a'; TheLineStyle = '-.'; TheMarker = 'none';
    d_temp = squeeze(set_of_______OD(the_propagation_index,i_case,1,:));
    m_temp = squeeze(set_of_______OD(the_propagation_index,i_case,2,:))./(set_of_conv_DPF(the_propagation_index,i_case,1)*squeeze(set_of_______OD(the_propagation_index,i_case,1,:))+set_of_conv_DPF(the_propagation_index,i_case,2));
    m_temp = (m_temp-mua)./mua*100;
    plot(d_temp,m_temp,'Color',TheColor,'DisplayName',[ThePreLegend ,': ',TheLegend],'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6,'LineWidth',2), hold on
    figure(100+i_case), subplot(1,2,2), plot(d_temp,m_temp,'Color',TheColor,'DisplayName',[ThePreLegend ,': ',TheLegend],'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6,'LineWidth',2), hold on, legend show, axis square, axis tight, ylim([-100 250]), xlabel('d (cm)'), ylabel('est. error (%)'), set(gca,'fontsize',12), axis square
    figure(i_fig), subplot(3,3,i_case);
    clearvars TheLegend TheLineStyle TheMarker d_temp m_temp
    TheLegend = '<s>/d \mu_a'; TheLineStyle = '--'; TheMarker = 'none';
    d_temp = squeeze(set_of_______OD(the_propagation_index,i_case,1,:));
    m_temp = squeeze(set_of_______OD(the_propagation_index,i_case,2,:))./squeeze(set_of_______OD(the_propagation_index,i_case,1,:))./squeeze(set_of_dfnd_DPF(the_propagation_index,i_case,2,:));
    m_temp = (m_temp-mua)./mua*100;
    plot(d_temp,m_temp,'Color',TheColor,'DisplayName',[ThePreLegend ,': ',TheLegend],'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6,'LineWidth',2), hold on
    clearvars TheLegend TheLineStyle TheMarker d_temp m_temp
    TheLegend = 'true \mu_a'; TheLineStyle = '-'; TheMarker = 'none';
    d_temp = squeeze(set_of_______OD(the_propagation_index,i_case,1,:));
    m_temp = squeeze(set_of_______OD(the_propagation_index,i_case,2,:))./squeeze(set_of_______OD(the_propagation_index,i_case,1,:))./squeeze(set_of_true_DPF(the_propagation_index,i_case,2,:));
    m_temp = (m_temp-mua)./mua*100;
    plot(d_temp,m_temp,'Color',TheColor,'DisplayName',[ThePreLegend ,': ',TheLegend],'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6,'LineWidth',2), hold on
    clearvars TheLegend TheLineStyle TheMarker d_temp m_temp
    clearvars the_propagation_index ThePreLegend
    % (transmittance)
    i_fig = 12; the_propagation_index = 2; figure(i_fig), subplot(3,3,i_case); ThePreLegend = 'transmittance';
    TheLegend = 'conv. \mu_a'; TheLineStyle = ':'; TheMarker = '.';
    d_temp = squeeze(set_of_______OD(the_propagation_index,i_case,1,:));
    m_temp = squeeze(set_of_______OD(the_propagation_index,i_case,2,:))./squeeze(set_of_______OD(the_propagation_index,i_case,1,:))./set_of_conv_DPF(the_propagation_index,i_case,3);
    m_temp = (m_temp-mua)./mua*100;
    plot(d_temp,m_temp,'Color',TheColor,'DisplayName',[ThePreLegend ,': ',TheLegend],'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6,'LineWidth',2), hold on
    TheLegend = 'Imp. conv. \mu_a'; TheLineStyle = '-.'; TheMarker = '.';
    d_temp = squeeze(set_of_______OD(the_propagation_index,i_case,1,:));
    m_temp = squeeze(set_of_______OD(the_propagation_index,i_case,2,:))./(set_of_conv_DPF(the_propagation_index,i_case,1)*squeeze(set_of_______OD(the_propagation_index,i_case,1,:))+set_of_conv_DPF(the_propagation_index,i_case,2));
    m_temp = (m_temp-mua)./mua*100;
    plot(d_temp,m_temp,'Color',TheColor,'DisplayName',[ThePreLegend ,': ',TheLegend],'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6,'LineWidth',2), hold on
    clearvars TheLegend TheLineStyle TheMarker d_temp m_temp
    TheLegend = '<s>/d \mu_a'; TheLineStyle = '--'; TheMarker = '.';
    d_temp = squeeze(set_of_______OD(the_propagation_index,i_case,1,:));
    m_temp = squeeze(set_of_______OD(the_propagation_index,i_case,2,:))./squeeze(set_of_______OD(the_propagation_index,i_case,1,:))./squeeze(set_of_dfnd_DPF(the_propagation_index,i_case,2,:));
    m_temp = (m_temp-mua)./mua*100;
    plot(d_temp,m_temp,'Color',TheColor,'DisplayName',[ThePreLegend ,': ',TheLegend],'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6,'LineWidth',2), hold on
    clearvars TheLegend TheLineStyle TheMarker d_temp m_temp
    TheLegend = 'true \mu_a'; TheLineStyle = '-'; TheMarker = '.';
    d_temp = squeeze(set_of_______OD(the_propagation_index,i_case,1,:));
    m_temp = squeeze(set_of_______OD(the_propagation_index,i_case,2,:))./squeeze(set_of_______OD(the_propagation_index,i_case,1,:))./squeeze(set_of_true_DPF(the_propagation_index,i_case,2,:));
    m_temp = (m_temp-mua)./mua*100;
    plot(d_temp,m_temp,'Color',TheColor,'DisplayName',[ThePreLegend ,': ',TheLegend],'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6,'LineWidth',2), hold on
    clearvars TheLegend TheLineStyle TheMarker d_temp m_temp
    xlim([0    10]), xlabel(x_label), set(gca,'xtick',0:1:10)
    ylim([-25 +25]), ylabel(y_label), set(gca,'ytick',-100:5:100)
    title(TheTitle), set(gca,'fontsize',12), grid on, hold on, % legend('show','Location','best')
    clearvars the_propagation_index ThePreLegend
    clearvars i_fig TheColor TheTitle x_label y_label



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
if Therloess<0 & length(x_bind)>3
    y_bind = smooth(x_bind,y_bind,Therloess,'rloess');
end
if blnShow
    figure(i_fig)
    subplot(3,3,i_subplot)
    plot(x_bind,TheOutFun(y_bind),'Color',TheColor,'DisplayName',TheLegend,'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor',TheEdgeColor,'MarkerFaceColor',TheColor,'MarkerSize',6,'LineWidth',1.0), hold on
end
if blnFit
    plot(x_bind,m*x_bind+b,'Color',TheColor,'HandleVisibility','off'), hold on
end
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
