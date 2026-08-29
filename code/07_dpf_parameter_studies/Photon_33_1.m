function [] = Photon_33_1 ()
% Repository group: 07_dpf_parameter_studies
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: w & s vs. d over mu_s & mu_a: more photons, more combinations
% the analysis: population stat

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
set_of_mua = (0.0:.05:0.5)+eps;
set_of_mus = (000:050:500)+eps;

% arrays to save variables
N_bins = 250; the_percent = 5;

N_diff = nan(length(set_of_mua),length(set_of_mus));
N_trns = nan(length(set_of_mua),length(set_of_mus));
N_absp = nan(length(set_of_mua),length(set_of_mus));
R_of_s_diff = nan(length(set_of_mua),length(set_of_mus),2);
R_of_s_trns = nan(length(set_of_mua),length(set_of_mus),2);

mIs_diff = nan(length(set_of_mua),length(set_of_mus),4);
mIs_trns = nan(length(set_of_mua),length(set_of_mus),4);
mId_diff = nan(length(set_of_mua),length(set_of_mus),4);
mId_trns = nan(length(set_of_mua),length(set_of_mus),4);
msd_diff = nan(length(set_of_mua),length(set_of_mus),3);
msd_trns = nan(length(set_of_mua),length(set_of_mus),3);

sum_N_diff = nan(length(set_of_mua),length(set_of_mus),2,N_bins);
sum_N_trns = nan(length(set_of_mua),length(set_of_mus),2,N_bins);
avg_w_diff = nan(length(set_of_mua),length(set_of_mus),2,N_bins);
avg_w_trns = nan(length(set_of_mua),length(set_of_mus),2,N_bins);
sum_I_diff = nan(length(set_of_mua),length(set_of_mus),2,N_bins);
sum_I_trns = nan(length(set_of_mua),length(set_of_mus),2,N_bins);
my_colormap = nan(length(set_of_mua),length(set_of_mus),3);

for i_a = 1:1:length(set_of_mua)
    for i_s = 1:1:length(set_of_mus)
        % read dbase
        mua = set_of_mua(i_a);
        mus = set_of_mus(i_s);
        the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_33_mua_',sprintf('%.2f',mua),'_mus_',sprintf('%.2f',mus),'.mat'];
        t_db = load(the_filename); clearvars the_filename
        disp(['mua = ',num2str(mua,'%.2f'),', mus = ',num2str(mus,'%.0f')])
        my_colormap(i_a,i_s,:) = get_color(mua/max(set_of_mua),mus/max(set_of_mus)); % make my own colormap

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
        d_diff_edges = linspace(0,sqrt((Lx/2).^2+(Ly/2).^2      ),N_bins+1); [~,~,ind_diff] = histcounts(sqrt(t_db.x(t_db.c==0).^2+t_db.y(t_db.c==0).^2+t_db.z(t_db.c==0).^2),d_diff_edges); % d_diffuse bins
        d_trns_edges = linspace(0,sqrt((Lx/2).^2+(Ly/2).^2+Lz.^2),N_bins+1); [~,~,ind_trns] = histcounts(sqrt(t_db.x(t_db.c==1).^2+t_db.y(t_db.c==1).^2+t_db.z(t_db.c==1).^2),d_trns_edges); % d_Transmittance bins
        clearvars d_trns_edges d_diff_edges



        % I vs. s (scatterplot)
        i_fig = 1; TheOutFun = @(x)(-log(x)); TheLineStyle = 'none'; TheMarker = 'o'; x_label = 's (cm)'; y_label = 'OD (a.u.)';
        TheLegend = ['\mu_a = ',num2str(mua,'%.2f'),', \mu_s = ',num2str(mus,'%.0f')];
        TheColor = get_color(mua/max(set_of_mua),mus/max(set_of_mus));
        % (diffuse)
        i_subplot = 1; index_in = ind_diff; TheCode = 0; TheTitle = 'diffuse';
        x_temp = t_db.s(t_db.c==TheCode); y_temp = t_db.w(t_db.c==TheCode);
        x_temp = x_temp; y_temp = TheOutFun(y_temp);
        if get_plot(mua,mus)
            figure(i_fig), subplot(1,2,i_subplot)
            plot(x_temp,y_temp,'Color',TheColor,'DisplayName',TheLegend,'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','w','MarkerFaceColor',TheColor,'MarkerSize',8), hold on
            xlim([0 150]), xlabel(x_label), ylabel(y_label), title(TheTitle), axis([0 50 0 5])
            set(gca,'fontsize',16), axis square, grid on, hold on
        end
        mdl = fitlm(x_temp,y_temp);
        mIs_diff(i_a,i_s,3) = mdl.Coefficients.pValue(2);
        if mIs_diff(i_a,i_s,3)<=0.05
            mIs_diff(i_a,i_s,2) = mdl.Coefficients.Estimate(1);
            mIs_diff(i_a,i_s,1) = mdl.Coefficients.Estimate(2);
            mIs_diff(i_a,i_s,4) = mIs_diff(i_a,i_s,1)/mua;
        end
        clearvars i_subplot index_in x_temp y_temp TheCode TheTitle mdl
        % (transmittance)
        i_subplot = 2; index_in = ind_trns; TheCode = 1; TheTitle = 'transmittance';
        x_temp = t_db.s(t_db.c==TheCode); y_temp = t_db.w(t_db.c==TheCode);
        x_temp = x_temp; y_temp = TheOutFun(y_temp);
        if get_plot(mua,mus)
            figure(i_fig), subplot(1,2,i_subplot)
            plot(x_temp,y_temp,'Color',TheColor,'DisplayName',TheLegend,'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','w','MarkerFaceColor',TheColor,'MarkerSize',8), hold on
            xlim([0 150]), xlabel(x_label), ylabel(y_label), title(TheTitle), axis([0 50 0 5])
            set(gca,'fontsize',16), axis square, grid on, hold on
        end
        mdl = fitlm(x_temp,y_temp);
        mIs_trns(i_a,i_s,3) = mdl.Coefficients.pValue(2);
        if mIs_trns(i_a,i_s,3)<=0.05
            mIs_trns(i_a,i_s,2) = mdl.Coefficients.Estimate(1);
            mIs_trns(i_a,i_s,1) = mdl.Coefficients.Estimate(2);
            mIs_trns(i_a,i_s,4) = mIs_trns(i_a,i_s,1)/mua;
        end
        clearvars i_subplot index_in x_temp y_temp TheCode TheTitle mdl
        clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label TheLegend TheColor



        % I vs. d
        i_fig = 2; fun_x = @mean; fun_y = @sum; TheOutFun = @(x)(-log(x)); TheLineStyle = 'none'; TheMarker = 'o'; x_label = 'd (cm)'; y_label = 'OD (a.u.)'; Therloess = 0.0; blnFit = false;
        % (diffuse)
        i_subplot = 1; index_in = ind_diff; TheCode = 0; TheTitle = 'diffuse';
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.w(t_db.c==TheCode)./t_db.no_of_photons;
        [mId_diff(i_a,i_s,1),mId_diff(i_a,i_s,2),mId_diff(i_a,i_s,3)] = make_a_subplot(i_fig,i_subplot,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,set_of_mua,set_of_mus,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit);
        mId_diff(i_a,i_s,4) = mId_diff(i_a,i_s,1)/mua;
        if ~isempty(index_in), sum_I_diff(i_a,i_s,1,1:max(index_in)) = accumarray(index_in,x_temp,[],fun_x,nan); end
        if ~isempty(index_in), sum_I_diff(i_a,i_s,2,1:max(index_in)) = accumarray(index_in,y_temp,[],fun_y,nan); end
        if get_plot(mua,mus)
            figure(i_fig), subplot(1,2,i_subplot), xlabel(x_label), ylabel(y_label), title(TheTitle), axis([0 20 0 20])
            set(gca,'fontsize',16), axis square, grid on, hold on
        end
        clearvars i_subplot index_in x_temp y_temp TheCode TheTitle
        % (transmittance)
        i_subplot = 2; index_in = ind_trns; TheCode = 1; TheTitle = 'transmittance';
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.w(t_db.c==TheCode)./t_db.no_of_photons;
        [mId_trns(i_a,i_s,1),mId_trns(i_a,i_s,2),mId_trns(i_a,i_s,3)] = make_a_subplot(i_fig,i_subplot,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,set_of_mua,set_of_mus,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit);
        mId_trns(i_a,i_s,4) = mId_trns(i_a,i_s,1)/mua;
        if ~isempty(index_in), sum_I_trns(i_a,i_s,1,1:max(index_in)) = accumarray(index_in,x_temp,[],fun_x,nan); end
        if ~isempty(index_in), sum_I_trns(i_a,i_s,2,1:max(index_in)) = accumarray(index_in,y_temp,[],fun_y,nan); end
        if get_plot(mua,mus)
            figure(i_fig), subplot(1,2,i_subplot), xlabel(x_label), ylabel(y_label), title(TheTitle), axis([0 20 0 20])
            set(gca,'fontsize',16), axis square, grid on, hold on
        end
        clearvars i_subplot index_in x_temp y_temp TheCode TheTitle
        clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label Therloess blnFit



        % min_{s} vs. d
        i_fig = 3; fun_x = @mean; fun_y = @(x)(prctile(x,the_percent)); TheOutFun = @(x)(x); TheLineStyle = '-'; TheMarker = 'none'; x_label = 'd (cm)'; y_label = 's_{min} (cm)'; Therloess = 0.5; blnFit = false;
        % (diffuse)
        i_subplot = 1; index_in = ind_diff; TheCode = 0; TheTitle = 'diffuse';
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.s(t_db.c==TheCode);
        [~,~,~] = make_a_subplot(i_fig,i_subplot,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,set_of_mua,set_of_mus,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit);
        if get_plot(mua,mus)
            figure(i_fig), subplot(1,2,i_subplot), xlabel(x_label), ylabel(y_label), title(TheTitle), axis([0 20 0 100])
            set(gca,'fontsize',16), axis square, grid on, hold on
        end
        clearvars i_subplot index_in x_temp y_temp TheCode TheTitle
        % (transmittance)
        i_subplot = 2; index_in = ind_trns; TheCode = 1; TheTitle = 'transmittance';
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.s(t_db.c==TheCode);
        [~,~,~] = make_a_subplot(i_fig,i_subplot,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,set_of_mua,set_of_mus,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit);
        if get_plot(mua,mus)
            figure(i_fig), subplot(1,2,i_subplot), xlabel(x_label), ylabel(y_label), title(TheTitle), axis([0 20 0 100])
            set(gca,'fontsize',16), axis square, grid on, hold on
        end
        clearvars i_subplot index_in x_temp y_temp TheCode TheTitle
        clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label Therloess blnFit

        % max_{s} vs. d
        i_fig = 3; fun_x = @mean; fun_y = @(x)(prctile(x,100-the_percent)); TheOutFun = @(x)(x); TheLineStyle = '-'; TheMarker = 'none'; x_label = 'd (cm)'; y_label = 's_{max} (cm)'; Therloess = 0.5; blnFit = false;
        % (diffuse)
        i_subplot = 1; index_in = ind_diff; TheCode = 0; TheTitle = 'diffuse';
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.s(t_db.c==TheCode);
        [~,~,~] = make_a_subplot(i_fig,i_subplot,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,set_of_mua,set_of_mus,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit);
        if get_plot(mua,mus)
            figure(i_fig), subplot(1,2,i_subplot), xlabel(x_label), ylabel(y_label), title(TheTitle), axis([0 20 0 100])
            set(gca,'fontsize',16), axis square, grid on, hold on
        end
        clearvars i_subplot index_in x_temp y_temp TheCode TheTitle
        % (transmittance)
        i_subplot = 2; index_in = ind_trns; TheCode = 1; TheTitle = 'transmittance';
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.s(t_db.c==TheCode);
        [~,~,~] = make_a_subplot(i_fig,i_subplot,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,set_of_mua,set_of_mus,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit);
        if get_plot(mua,mus)
            figure(i_fig), subplot(1,2,i_subplot), xlabel(x_label), ylabel(y_label), title(TheTitle), axis([0 20 0 100])
            set(gca,'fontsize',16), axis square, grid on, hold on
        end
        clearvars i_subplot index_in x_temp y_temp TheCode TheTitle
        clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label Therloess blnFit



        % s vs. d
        i_fig = 4; fun_x = @mean; fun_y = @mean; TheOutFun = @(x)(x); TheLineStyle = '-'; TheMarker = 'none'; x_label = 'd (cm)'; y_label = 's (cm)'; Therloess = 0.5; blnFit = false;
        % (diffuse)
        i_subplot = 1; index_in = ind_diff; TheCode = 0; TheTitle = 'diffuse';
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.s(t_db.c==TheCode);
        [msd_diff(i_a,i_s,1),msd_diff(i_a,i_s,2),msd_diff(i_a,i_s,3)] = make_a_subplot(i_fig,i_subplot,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,set_of_mua,set_of_mus,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit);
        if get_plot(mua,mus)
            figure(i_fig), subplot(1,2,i_subplot), xlabel(x_label), ylabel(y_label), title(TheTitle), axis([0 20 0 100])
            set(gca,'fontsize',16), axis square, grid on, hold on
        end
        clearvars i_subplot index_in x_temp y_temp TheCode TheTitle
        % (transmittance)
        i_subplot = 2; index_in = ind_trns; TheCode = 1; TheTitle = 'transmittance';
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.s(t_db.c==TheCode);
        [msd_trns(i_a,i_s,1),msd_trns(i_a,i_s,2),msd_trns(i_a,i_s,3)] = make_a_subplot(i_fig,i_subplot,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,set_of_mua,set_of_mus,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit);
        if get_plot(mua,mus)
            figure(i_fig), subplot(1,2,i_subplot), xlabel(x_label), ylabel(y_label), title(TheTitle), axis([0 20 0 100])
            set(gca,'fontsize',16), axis square, grid on, hold on
        end
        clearvars i_subplot index_in x_temp y_temp TheCode TheTitle
        clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label Therloess blnFit



        % DPF vs. d
        i_fig = 5; fun_x = @mean; fun_y = @mean; TheOutFun = @(x)(x); TheLineStyle = '-'; TheMarker = 'none'; x_label = 'd (cm)'; y_label = 'DPF (a.u.)'; Therloess = 0.5; blnFit = false;
        % (diffuse)
        i_subplot = 1; index_in = ind_diff; TheCode = 0; TheTitle = 'diffuse';
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.s(t_db.c==TheCode)./x_temp;
        [~,~,~] = make_a_subplot(i_fig,i_subplot,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,set_of_mua,set_of_mus,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit);
        if get_plot(mua,mus) && ~isnan(msd_diff(i_a,i_s,3))
            % x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
            % plot(x_bind,msd_diff(i_a,i_s,1)*ones(size(x_bind)),'Color',get_color(mua/max(set_of_mua),mus/max(set_of_mus)),'HandleVisibility','off','LineStyle',':'), hold on
            % clearvars x_bind
            figure(i_fig), subplot(1,2,i_subplot), xlabel(x_label), ylabel(y_label), title(TheTitle), axis([0 20 0 25])
            set(gca,'fontsize',16), axis square, grid on, hold on
        end
        clearvars i_subplot index_in x_temp y_temp TheCode TheTitle
        % (transmittance)
        i_subplot = 2; index_in = ind_trns; TheCode = 1; TheTitle = 'transmittance';
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.s(t_db.c==TheCode)./x_temp;
        [~,~,~] = make_a_subplot(i_fig,i_subplot,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,set_of_mua,set_of_mus,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit);
        if get_plot(mua,mus) && ~isnan(msd_trns(i_a,i_s,3))
            % x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
            % plot(x_bind,msd_trns(i_a,i_s,1)*ones(size(x_bind)),'Color',get_color(mua/max(set_of_mua),mus/max(set_of_mus)),'HandleVisibility','off','LineStyle',':'), hold on
            % clearvars x_bind
            figure(i_fig), subplot(1,2,i_subplot), xlabel(x_label), ylabel(y_label), title(TheTitle), axis([0 20 0 25])
            set(gca,'fontsize',16), axis square, grid on, hold on
        end
        clearvars i_subplot index_in x_temp y_temp TheCode TheTitle
        clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label Therloess blnFit



        % histogram of s
        i_fig = 6;
        % (diffuse)
        i_subplot = 1; TheCode = 0; TheTitle = 'diffuse';
        y_temp = t_db.s(t_db.c==TheCode);
        R_of_s_diff(i_a,i_s,:) = [prctile(y_temp,the_percent) , prctile(y_temp,100-the_percent)];
        [n_temp,x_temp] = histcounts(y_temp,'BinEdges',linspace(0,150)); x_temp = 1/2*(x_temp(1:end-1)+x_temp(2:end-0)); n_temp = n_temp./t_db.no_of_photons*100;
        if ~isempty(y_temp)>0 & get_plot(mua,mus)
            figure(i_fig), subplot(1,2,i_subplot), semilogy(x_temp,n_temp,'Color',get_color(mua/max(set_of_mua),mus/max(set_of_mus)),'LineWidth',1,'DisplayName',['\mu_a = ',num2str(mua,'%.2f'),', \mu_s = ',num2str(mus,'%.0f')]), hold on
            xlabel('s (cm)'), ylabel('freq (%)'), title(TheTitle), axis([0 100 0 100])
            set(gca,'fontsize',16), axis square, grid on, hold on
        end
        clearvars i_subplot index_in x_temp y_temp n_temp TheCode TheTitle
        % (transmittance)
        i_subplot = 2; TheCode = 1; TheTitle = 'transmittance';
        y_temp = t_db.s(t_db.c==TheCode);
        R_of_s_trns(i_a,i_s,:) = [prctile(y_temp,the_percent) , prctile(y_temp,100-the_percent)];
        [n_temp,x_temp] = histcounts(y_temp,'BinEdges',linspace(0,150)); x_temp = 1/2*(x_temp(1:end-1)+x_temp(2:end-0)); n_temp = n_temp./t_db.no_of_photons*100;
        if ~isempty(y_temp)>0 & get_plot(mua,mus)
            figure(i_fig), subplot(1,2,i_subplot), semilogy(x_temp,n_temp,'Color',get_color(mua/max(set_of_mua),mus/max(set_of_mus)),'LineWidth',1,'DisplayName',['\mu_a = ',num2str(mua,'%.2f'),', \mu_s = ',num2str(mus,'%.0f')]), hold on
            xlabel('s (cm)'), ylabel('freq (%)'), title(TheTitle), axis([0 100 0 100])
            set(gca,'fontsize',16), axis square, grid on, hold on
        end
        clearvars i_subplot index_in x_temp y_temp n_temp TheCode TheTitle
        clearvars i_fig



        % N vs. d
        i_fig = 7; fun_x = @mean; fun_y = @sum; TheOutFun = @(x)(x); TheLineStyle = '-'; TheMarker = 'none'; x_label = 'd (cm)'; y_label = 'N (#)'; Therloess = 0.0; blnFit = false;
        % (diffuse)
        i_subplot = 1; index_in = ind_diff; TheCode = 0; TheTitle = 'diffuse';
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = ones(size(t_db.s(t_db.c==TheCode)))./t_db.no_of_photons*100;
        [~,~,~] = make_a_subplot(i_fig,i_subplot,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,set_of_mua,set_of_mus,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit);
        if ~isempty(index_in), sum_N_diff(i_a,i_s,1,1:max(index_in)) = accumarray(index_in,x_temp,[],fun_x,nan); end
        if ~isempty(index_in), sum_N_diff(i_a,i_s,2,1:max(index_in)) = accumarray(index_in,y_temp,[],fun_y,nan)/100; end
        if get_plot(mua,mus)
            figure(i_fig), subplot(1,2,i_subplot), xlabel(x_label), ylabel(y_label), yscale("log"), title(TheTitle), axis([0 20 0 5e1])
            set(gca,'fontsize',16), axis square, grid on, hold on
        end
        clearvars i_subplot index_in x_temp y_temp TheCode TheTitle
        % (transmittance)
        i_subplot = 2; index_in = ind_trns; TheCode = 1; TheTitle = 'transmittance';
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = ones(size(t_db.s(t_db.c==TheCode)))./t_db.no_of_photons*100;
        [~,~,~] = make_a_subplot(i_fig,i_subplot,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,set_of_mua,set_of_mus,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit);
        if ~isempty(index_in), sum_N_trns(i_a,i_s,1,1:max(index_in)) = accumarray(index_in,x_temp,[],fun_x,nan); end
        if ~isempty(index_in), sum_N_trns(i_a,i_s,2,1:max(index_in)) = accumarray(index_in,y_temp,[],fun_y,nan)/100; end
        if get_plot(mua,mus)
            figure(i_fig), subplot(1,2,i_subplot), xlabel(x_label), ylabel(y_label), yscale("log"), title(TheTitle), axis([0 20 0 5e1])
            set(gca,'fontsize',16), axis square, grid on, hold on
        end
        clearvars i_subplot index_in x_temp y_temp TheCode TheTitle
        clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label Therloess blnFit



        % w_bar vs. d
        i_fig = 8; fun_x = @mean; fun_y = @mean; TheOutFun = @(x)(-log(x)); TheLineStyle = '-'; TheMarker = 'none'; x_label = 'd (cm)'; y_label = 'w_{bar} (a.u.)'; Therloess = 0.5; blnFit = false;
        % (diffuse)
        i_subplot = 1; index_in = ind_diff; TheCode = 0; TheTitle = 'diffuse';
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.w(t_db.c==TheCode);
        [~,~,~] = make_a_subplot(i_fig,i_subplot,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,set_of_mua,set_of_mus,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit);
        if ~isempty(index_in), avg_w_diff(i_a,i_s,1,1:max(index_in)) = accumarray(index_in,x_temp,[],fun_x,nan); end
        if ~isempty(index_in), avg_w_diff(i_a,i_s,2,1:max(index_in)) = accumarray(index_in,y_temp,[],fun_y,nan); end
        if get_plot(mua,mus)
            figure(i_fig), subplot(1,2,i_subplot), xlabel(x_label), ylabel(y_label), title(TheTitle), axis([0 20 0 5])
            set(gca,'fontsize',16), axis square, grid on, hold on
        end
        clearvars i_subplot index_in x_temp y_temp TheCode TheTitle
        % (transmittance)
        i_subplot = 2; index_in = ind_trns; TheCode = 1; TheTitle = 'transmittance';
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.w(t_db.c==TheCode);
        [~,~,~] = make_a_subplot(i_fig,i_subplot,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,set_of_mua,set_of_mus,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit);
        if ~isempty(index_in), avg_w_trns(i_a,i_s,1,1:max(index_in)) = accumarray(index_in,x_temp,[],fun_x,nan); end
        if ~isempty(index_in), avg_w_trns(i_a,i_s,2,1:max(index_in)) = accumarray(index_in,y_temp,[],fun_y,nan); end
        if get_plot(mua,mus)
            figure(i_fig), subplot(1,2,i_subplot), xlabel(x_label), ylabel(y_label), title(TheTitle), axis([0 20 0 5])
            set(gca,'fontsize',16), axis square, grid on, hold on
        end
        clearvars i_subplot index_in x_temp y_temp TheCode TheTitle
        clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label Therloess blnFit



        % corrected DPF vs. d
        i_fig = 9; fun_x = @mean; fun_y = @(x)(sum(exp(-mua.*x))./t_db.no_of_photons); TheOutFun = @(x)(log(x)./mua); TheLineStyle = '-'; TheMarker = 'none'; x_label = 'd (cm)'; y_label = 'corrected DPF (a.u.)'; Therloess = 0.5; blnFit = false;
        % (diffuse)
        i_subplot = 1; index_in = ind_diff; TheCode = 0; TheTitle = 'diffuse';
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.s(t_db.c==TheCode);
        x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
        y_bind = accumarray(index_in,y_temp,[],fun_y,nan); y_bind = -TheOutFun(y_bind)./x_bind;
        TheColor = get_color(mua/max(set_of_mua),mus/max(set_of_mus));
        TheLegend = ['\mu_a = ',num2str(mua,'%.2f'),', \mu_s = ',num2str(mus,'%.0f')];
        if get_plot(mua,mus)
            figure(i_fig), subplot(1,2,i_subplot), xlabel(x_label), ylabel(y_label), yscale("log"), title(TheTitle), axis([0 10 0 150])
            plot(x_bind,TheOutFun(y_bind),'Color',TheColor,'DisplayName',TheLegend,'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','w','MarkerFaceColor',TheColor,'MarkerSize',8,'LineWidth',1.5), hold on
            set(gca,'fontsize',16), axis square, grid on, hold on
        end
        clearvars i_subplot index_in x_temp y_temp x_bind y_bind TheCode TheTitle TheColor TheLegend
        % (transmittance)
        i_subplot = 2; index_in = ind_trns; TheCode = 1; TheTitle = 'transmittance';
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.s(t_db.c==TheCode);
        x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
        y_bind = accumarray(index_in,y_temp,[],fun_y,nan); y_bind = -TheOutFun(y_bind)./x_bind;
        TheColor = get_color(mua/max(set_of_mua),mus/max(set_of_mus));
        TheLegend = ['\mu_a = ',num2str(mua,'%.2f'),', \mu_s = ',num2str(mus,'%.0f')];
        if get_plot(mua,mus)
            figure(i_fig), subplot(1,2,i_subplot), xlabel(x_label), ylabel(y_label), yscale("log"), title(TheTitle), axis([0 10 0 150])
            plot(x_bind,TheOutFun(y_bind),'Color',TheColor,'DisplayName',TheLegend,'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','w','MarkerFaceColor',TheColor,'MarkerSize',8,'LineWidth',1.5), hold on
            set(gca,'fontsize',16), axis square, grid on, hold on
        end
        clearvars i_subplot index_in x_temp y_temp x_bind y_bind TheCode TheTitle TheColor TheLegend
        clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label Therloess blnFit



        N_diff(i_a,i_s) = sum(t_db.c==0)./t_db.no_of_photons;
        N_trns(i_a,i_s) = sum(t_db.c==1)./t_db.no_of_photons;
        N_absp(i_a,i_s) = 1 - N_diff(i_a,i_s) - N_trns(i_a,i_s);
        clearvars t_db the_filename mua mus ind_diff ind_trns
    end
end

figure(10)
[X_lgnd,Y_lgnd] = ndgrid(set_of_mua*1000,set_of_mus);
image(set_of_mua,set_of_mus,permute(my_colormap,[2,1,3])), shading flat
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
title('legend'), set(gca,'fontsize',16), axis square, axis tight, axis xy

figure(11)
subplot(1,3,1), pcolor(set_of_mua,set_of_mus,(N_diff+eps).'), shading flat
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
title('percent of diffused photons'),
set(gca,'fontsize',16), axis square, axis tight,
colormap jet; clim([0 1]); % set(gca,'colorscale','log');
subplot(1,3,2), pcolor(set_of_mua,set_of_mus,(N_absp+eps).'), shading flat
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
title('percent of absorbed photons'),
set(gca,'fontsize',16), axis square, axis tight,
colormap jet; clim([0 1]); % set(gca,'colorscale','log');
subplot(1,3,3), pcolor(set_of_mua,set_of_mus,(N_trns+eps).'), shading flat
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
title('percent of transmitted photons'),
set(gca,'fontsize',16), axis square, axis tight,
colormap jet; h = colorbar; clim([0 1]); % set(gca,'colorscale','log');
ylabel(h,'% (a.u.)','FontSize',16), clearvars h

figure(12)
subplot(1,2,1), pcolor(set_of_mua,set_of_mus,squeeze(msd_diff(:,:,1)).'), shading flat
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
title('s/d for diffused photons'),
set(gca,'fontsize',16), axis square, axis tight,
colormap jet;
clim([0 10])
subplot(1,2,2), pcolor(set_of_mua,set_of_mus,squeeze(msd_trns(:,:,1)).'), shading flat
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
title('s/d for transmitted photons'),
set(gca,'fontsize',16), axis square, axis tight,
colormap jet; h = colorbar;
ylabel(h,'s/d (a.u.)','FontSize',16), clim([0 10]), clearvars h

figure(13)
subplot(1,2,1), pcolor(set_of_mua,set_of_mus,squeeze(mIs_diff(:,:,1)).'), shading flat
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
title('OD/s for diffused photons'),
set(gca,'fontsize',16), axis square, axis tight,
colormap jet; clim([min(set_of_mua) max(set_of_mua)])
subplot(1,2,2), pcolor(set_of_mua,set_of_mus,squeeze(mIs_trns(:,:,1)).'), shading flat
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
title('OD/s for transmitted photons'),
set(gca,'fontsize',16), axis square, axis tight,
colormap jet; h = colorbar; clim([min(set_of_mua) max(set_of_mua)])
ylabel(h,'OD/s (cm^{-1})','FontSize',16), clearvars h

figure(14)
subplot(1,2,1), pcolor(set_of_mua,set_of_mus,squeeze(mId_diff(:,:,1)).'), shading flat
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
title('OD/d for diffused photons'),
set(gca,'fontsize',16), axis square, axis tight,
colormap jet;
subplot(1,2,2), pcolor(set_of_mua,set_of_mus,squeeze(mId_trns(:,:,1)).'), shading flat
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
title('OD/d for transmitted photons'),
set(gca,'fontsize',16), axis square, axis tight,
colormap jet; h = colorbar; clim([0 8])
ylabel(h,'OD/d (cm^{-1})','FontSize',16), clearvars h

% for i_a = 2:3:length(set_of_mua)
%     for i_s = 2:3:length(set_of_mus)
%         figure
%         TheColor = get_color(set_of_mua(i_a)/max(set_of_mua),set_of_mus(i_s)/max(set_of_mus));
%         TheLegend = ['\mu_a = ',num2str(set_of_mua(i_a),'%.2f'),', \mu_s = ',num2str(set_of_mus(i_s),'%.0f')];
%         subplot(1,2,1)
%         x_N = squeeze(sum_N_diff(i_a,i_s,1,:));
%         y_N = squeeze(sum_N_diff(i_a,i_s,2,:));
%         x_w = squeeze(avg_w_diff(i_a,i_s,1,:));
%         y_w = squeeze(avg_w_diff(i_a,i_s,2,:));
%         x_I = squeeze(sum_I_diff(i_a,i_s,1,:));
%         y_I = squeeze(sum_I_diff(i_a,i_s,2,:));
%         plot(1/2*(x_N+x_w),-log(y_N.*y_w),'Color',TheColor,'DisplayName',['indirect: ',TheLegend],'Marker','none','LineStyle','-','LineWidth',1), hold on
%         plot(x_I,          -log(y_I),     'Color',TheColor,'DisplayName',['direct: ',TheLegend],  'Marker','^',   'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor',TheColor), hold on
%         xlabel('d (cm)'), ylabel('OD (a.u.)'), title('diffuse'), axis([0 20 0 20])
%         set(gca,'fontsize',16), axis square, grid on, hold off, legend show
%         clearvars x_N y_N x_w y_w x_I y_I
%         subplot(1,2,2)
%         x_N = squeeze(sum_N_trns(i_a,i_s,1,:));
%         y_N = squeeze(sum_N_trns(i_a,i_s,2,:));
%         x_w = squeeze(avg_w_trns(i_a,i_s,1,:));
%         y_w = squeeze(avg_w_trns(i_a,i_s,2,:));
%         x_I = squeeze(sum_I_trns(i_a,i_s,1,:));
%         y_I = squeeze(sum_I_trns(i_a,i_s,2,:));
%         plot(1/2*(x_N+x_w),-log(y_N.*y_w),'Color',TheColor,'DisplayName',['indirect: ',TheLegend],'Marker','none','LineStyle','-','LineWidth',1), hold on
%         plot(x_I,          -log(y_I),     'Color',TheColor,'DisplayName',['direct: ',TheLegend],  'Marker','^',   'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor',TheColor), hold on
%         xlabel('d (cm)'), ylabel('OD (a.u.)'), title('transmittance'), axis([0 20 0 20])
%         set(gca,'fontsize',16), axis square, grid on, hold off, legend show
%         clearvars x_N y_N x_w y_w x_I y_I
%     end
% end
for i_fig = 1:14
    figure(i_fig)
    set(gcf,'Renderer','painters');
    saveas(gcf, ['Fig-',sprintf('%2.0f',i_fig),'.eps'],'epsc'); % Save as EPS with 300 DPI
    saveas(gcf, ['Fig-',sprintf('%2.0f',i_fig),'.jpg'],'jpeg'); % Save as EPS with 300 DPI
end
end

function [m,b,p] = make_a_subplot(i_fig,i_subplot,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,set_of_mua,set_of_mus,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit)
x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
y_bind = accumarray(index_in,y_temp,[],fun_y,nan);
TheColor = get_color(mua/max(set_of_mua),mus/max(set_of_mus));
TheLegend = ['\mu_a = ',num2str(mua,'%.2f'),', \mu_s = ',num2str(mus,'%.0f')];
mdl = fitlm(x_bind,TheOutFun(y_bind));
p = mdl.Coefficients.pValue(2);
if p<=0.05
    b = mdl.Coefficients.Estimate(1);
    m = mdl.Coefficients.Estimate(2);
else
    b = nan;
    m = nan;
end
if get_plot(mua,mus)
    figure(i_fig)
    subplot(1,2,i_subplot)
    if Therloess<0 & length(x_bind)>3
        y_bind = smooth(x_bind,y_bind,Therloess,'rloess');
    end
    plot(x_bind,TheOutFun(y_bind),'Color',TheColor,'DisplayName',TheLegend,'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','w','MarkerFaceColor',TheColor,'MarkerSize',8,'LineWidth',1.5), hold on
    if blnFit
        plot(x_bind,m*x_bind+b,'Color',TheColor,'HandleVisibility','off'), hold on
    end
end
end% out = (mua==0.10||mua==0.30||mua==0.50)&&(mus==010||mus==300||mus==500);

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
function [out] = get_plot(mua,mus)
mua = round(mua,2);
mus = round(mus,0);
out = (mua==0.05||mua==0.15||mua==0.25||mua==0.35||mua==0.45)&&(mus==50||mus==150||mus==250||mus==350||mus==450);
% out = (mua==0.00||mua==0.25||mua==0.50)&&(mus==000||mus==250||mus==500);
% out = true;
% out = (mua==0.10||mua==0.30||mua==0.50)&&(mus==100||mus==300||mus==500);
% out = (mua==0.05||mua==0.25||mua==0.45)&&(mus== 50||mus==250||mus==450);
% out = (mua==0.10||mua==0.30)&&(mus==100||mus==300);
% out = (mua==0.4)&&(mus==200||mus==500);
% out = (mua==0.05)&&(mus==50||mus==100);
end
