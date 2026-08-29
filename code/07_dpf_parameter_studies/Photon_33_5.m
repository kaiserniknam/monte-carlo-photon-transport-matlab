function [] = Photon_33_5 ()
% Repository group: 07_dpf_parameter_studies
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: w & s vs. d over mu_s & mu_a: more photons, more combinations
% the analysis: population stat - to plot p and q for mu_a's and mu_s's
% Optics Letters DPF Test (based on dataset from 39)
% figures for DPF manuscript

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

n = 1.4; g = 0.95;
set_of_mua = (.00:.05:0.5);
set_of_mus = (000:050:500);
number_of_all_photons = [];
N_bins = 200;
d_edges = linspace(0,sqrt((Lx/2).^2+(Ly/2).^2),N_bins+1); clearvars Lx Ly Lz
d_cntrs = 1/2*(d_edges(1:end-1)+d_edges(2:end-0));
set_of_cnst_DPF = nan(length(set_of_mua),length(set_of_mus),1);
set_of_idst_DPF = nan(length(set_of_mua),length(set_of_mus),2);
set_of_savg_DPF = nan(length(set_of_mua),length(set_of_mus),N_bins);
set_of_slop_DPF = nan(length(set_of_mua),length(set_of_mus),N_bins);
set_of_true_DPF = nan(length(set_of_mua),length(set_of_mus),N_bins);
set_of________I = nan(length(set_of_mua),length(set_of_mus),N_bins);
set_of_dpfs_err = nan(length(set_of_mua),length(set_of_mus),6,N_bins);
clearvars N_bins

p_est = nan(size(squeeze(set_of_idst_DPF(:,:,1))));
q_est = nan(size(squeeze(set_of_idst_DPF(:,:,2))));

p_anc = nan(size(squeeze(set_of_idst_DPF(:,:,1))));
q_anc = nan(size(squeeze(set_of_idst_DPF(:,:,2))));

for i_a = 1:length(set_of_mua)
    for i_s = 1:length(set_of_mus)
        % p_est(i_a,i_s) = exp(-0.707 - 0.507*log(set_of_mua(i_a)) + 0.506*log(set_of_mus(i_s)));
        % q_est(i_a,i_s) = exp(+1.135 - 1.370*log(set_of_mua(i_a)) - 0.198*log(set_of_mus(i_s)));

        p_est(i_a,i_s) = 0.493.*((set_of_mua(i_a)).^(-0.507)).*((set_of_mus(i_s)).^(+0.506));
        q_est(i_a,i_s) = 3.112.*((set_of_mua(i_a)).^(-1.370)).*((set_of_mus(i_s)).^(-0.198));

        mua  = set_of_mua(i_a);
        musp = set_of_mus(i_s)*(1-g);
        D = 1./(3.*(mua+musp));
        mueff = sqrt(3.*mua.*(mua+musp));
        z0 = 1./(mua+musp);
        Reff = -1.440./n.^2 + 0.710./n + 0.668 + 0.0636.*n;
        A = (1 + Reff)./(1 - Reff);
        zb = 2.*A.*D;

        p_anc(i_a,i_s) = mueff/mua;
        q_anc(i_a,i_s) = -log(1/4/pi/D*z0*(z0+2*zb))/mua;
        clearvars mua musp D mueff z0 Reff A zb
    end
end
clearvars i_a i_s

for i_a = 1:length(set_of_mua)
    for i_s = 1:length(set_of_mus)
        % read dbase
        mua = set_of_mua(i_a);
        mus = set_of_mus(i_s);
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
        number_of_all_photons = [number_of_all_photons ; sum(u==u_unique(i_fx))];
        t_db.no_of_photons = t_db.no_of_photons-sum(u~=u_unique(i_fx));
        clearvars u u_unique freq i_f i_fx

        % 1-D sorting
        [~,~,ind_diff] = histcounts(sqrt(t_db.x(t_db.c==0).^2+t_db.y(t_db.c==0).^2+t_db.z(t_db.c==0).^2),d_edges); % d_diffuse bins
        clearvars d_trns_edges d_diff_edges

        % s vs. d
        fun_x = @mean; fun_y = @mean;
        index_in = ind_diff; TheCode = 0;
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2);
        y_temp = t_db.s(t_db.c==TheCode);
        x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
        y_bind = accumarray(index_in,y_temp,[],fun_y,nan);
        set_of_savg_DPF(i_a,i_s,1:length(y_bind)) = y_bind./x_bind;
        clearvars fun_x fun_y index_in TheCode x_temp y_temp x_bind y_bind

        % slope DPF
        fun_x = @mean; fun_s = @(x)(sum(x.*exp(-mua.*x))./sum(exp(-mua.*x)));
        index_in = ind_diff; TheCode = 0;
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2);
        s_temp = t_db.s(t_db.c==TheCode);
        x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
        s_bind = accumarray(index_in,s_temp,[],fun_s,nan); s_bind = s_bind./x_bind;
        set_of_slop_DPF(i_a,i_s,1:length(s_bind)) = s_bind;
        clearvars fun_x fun_s index_in TheCode x_temp s_temp x_bind s_bind

        % I vs. d
        fun_x = @mean; fun_y = @sum; fun_s = @(x)(sum(exp(-mua.*x))./t_db.no_of_photons); TheOutFun = @(x)(-log(x));
        index_in = ind_diff; TheCode = 0;
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2);
        y_temp = t_db.w(t_db.c==TheCode)./t_db.no_of_photons;
        s_temp = t_db.s(t_db.c==TheCode);
        x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
        y_bind = accumarray(index_in,y_temp,[],fun_y,nan); y_bind = TheOutFun(y_bind);
        s_bind = accumarray(index_in,s_temp,[],fun_s,nan); s_bind = TheOutFun(s_bind)./mua./x_bind;

        mdl = fitlm(x_bind(~isnan(x_bind)),y_bind(~isnan(x_bind)));
        p = mdl.Coefficients.pValue(2);
        if p<=0.05
            set_of_cnst_DPF(i_a,i_s,1) = (x_bind(~isnan(x_bind))\y_bind(~isnan(x_bind)))/mua;
            set_of_idst_DPF(i_a,i_s,1) = mdl.Coefficients.Estimate(2)/mua;
            set_of_idst_DPF(i_a,i_s,2) = mdl.Coefficients.Estimate(1)/mua;
        else
            set_of_cnst_DPF(i_a,i_s,1) = nan;
            set_of_idst_DPF(i_a,i_s,1) = nan;
            set_of_idst_DPF(i_a,i_s,2) = nan;
        end
        set_of________I(i_a,i_s,1:length(y_bind)) = y_bind;
        set_of_true_DPF(i_a,i_s,1:length(s_bind)) = s_bind;

        % error
        set_of_dpfs_err(i_a,i_s,1,1:length(x_bind)) = ((y_bind./x_bind./set_of_cnst_DPF(i_a,i_s,1))-mua)./mua; % constant
        set_of_dpfs_err(i_a,i_s,2,1:length(x_bind)) = ((y_bind./x_bind./squeeze(set_of_savg_DPF(i_a,i_s,1:length(x_bind))))-mua)./mua; % <s>/d
        set_of_dpfs_err(i_a,i_s,3,1:length(x_bind)) = ((y_bind./x_bind./squeeze(set_of_true_DPF(i_a,i_s,1:length(x_bind))))-mua)./mua; % true
        set_of_dpfs_err(i_a,i_s,4,1:length(x_bind)) = ((y_bind./x_bind./(set_of_idst_DPF(i_a,i_s,1)+set_of_idst_DPF(i_a,i_s,2)./x_bind))-mua)./mua; % iverse-distance
        set_of_dpfs_err(i_a,i_s,5,1:length(x_bind)) = ((y_bind./x_bind./(p_est(i_a,i_s)            +q_est(i_a,i_s)            ./x_bind))-mua)./mua; % empirical
        set_of_dpfs_err(i_a,i_s,6,1:length(x_bind)) = ((y_bind./x_bind./(p_anc(i_a,i_s)            +q_anc(i_a,i_s)            ./x_bind))-mua)./mua; % analytical
        set_of_dpfs_err(i_a,i_s,7,1:length(x_bind)) = ((y_bind./x_bind./squeeze(set_of_slop_DPF(i_a,i_s,1:length(x_bind))))-mua)./mua; % slope DPF

        clearvars mdl p blnFit blnShow fun_x fun_y fun_s i_fig index_in TheCode TheOutFun x_temp y_temp s_temp x_bind y_bind s_bind
        clearvars t_db the_filename mua mus ind_diff ind_trns
    end
end
clearvars i_a i_s

figure(1)
i_a = 4; i_s = 4;
subplot(1,2,1)
plot(d_cntrs,squeeze(set_of_cnst_DPF(i_a,i_s,1)).*ones(size(d_cntrs)),      'DisplayName','constant','LineWidth',2), hold on
plot(d_cntrs,squeeze(set_of_savg_DPF(i_a,i_s,:)),                           'DisplayName','average', 'LineWidth',2), hold on
plot(d_cntrs,squeeze(set_of_slop_DPF(i_a,i_s,:)),                           'DisplayName','slope',   'LineWidth',2), hold on
plot(d_cntrs,squeeze(set_of_true_DPF(i_a,i_s,:)),                           'DisplayName','true',    'LineWidth',2), hold on
plot(d_cntrs,set_of_idst_DPF(i_a,i_s,1)+set_of_idst_DPF(i_a,i_s,2)./d_cntrs,'DisplayName','inverse', 'LineWidth',2), hold on
xlabel('d (cm)'),     xlim([0 6]),   set(gca,'xtick',0:1:6)
ylabel('DPF (a.u.)'), ylim([0 150]), set(gca,'ytick',0:25:150)
title(['Diff. DPF definitions for \mu_a = ',num2str(set_of_mua(i_a)),', \mu_s = ',num2str(set_of_mus(i_s))])
hold off, legend('show','Location','northeast'), set(gca,'fontsize',16)
subplot(1,2,2)
plot(d_cntrs,squeeze(set_of_dpfs_err(i_a,i_s,1,:)).*100,                    'DisplayName','constant','LineWidth',2), hold on
plot(d_cntrs,squeeze(set_of_dpfs_err(i_a,i_s,2,:)).*100,                    'DisplayName','average', 'LineWidth',2), hold on
plot(d_cntrs,squeeze(set_of_dpfs_err(i_a,i_s,7,:)).*100,                    'DisplayName','slope',   'LineWidth',2), hold on
plot(d_cntrs,squeeze(set_of_dpfs_err(i_a,i_s,3,:)).*100,                    'DisplayName','true',    'LineWidth',2), hold on
plot(d_cntrs,squeeze(set_of_dpfs_err(i_a,i_s,4,:)).*100,                    'DisplayName','inverse', 'LineWidth',2), hold on
xlabel('d (cm)'),     xlim([0 6]),   set(gca,'xtick',0:1:6)
ylabel('estimation error (%)'), ylim([-50 +250]), set(gca,'ytick',-50:50:250)
title(['Diff. DPF error generation for \mu_a = ',num2str(set_of_mua(i_a)),', \mu_s = ',num2str(set_of_mus(i_s))])
hold off, legend('show','Location','northeast'), set(gca,'fontsize',16)

figure(2)
subplot(2,1,1), surfc(set_of_mua,set_of_mus,squeeze(set_of_idst_DPF(:,:,1)).'), shading flat
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
zlabel('a (a.u.)'), set(gca,'ztick',0:25:75), title('p')
set(gca,'fontsize',16), axis square, axis tight, colormap jet; view([145 22.5])
axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus) 0 75])
subplot(2,1,2), surfc(set_of_mua,set_of_mus,squeeze(set_of_idst_DPF(:,:,2)).'), shading flat
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
zlabel('b (a.u.)'), set(gca,'ztick',0:25:75), title('q')
set(gca,'fontsize',16), axis square, axis tight, colormap jet; view([145 22.5])
axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus) 0 75])
clearvars h

figure(3)
idx = 1<=d_cntrs&d_cntrs<=6; icol = -1:0.5:+2;
subplot(2,2,1)
contourf(set_of_mua,set_of_mus,squeeze(nanmean(set_of_dpfs_err(:,:,1,idx)*100,4)).',icol,'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
h = colorbar; ylabel(h,'relative (%)'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',icol), title('constant DPF'); clim([min(icol) max(icol)])
set(gca,'fontsize',16), axis square, axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus)]); colormap jet
subplot(2,2,3)
contourf(set_of_mua,set_of_mus,squeeze(nanmean(set_of_dpfs_err(:,:,4,idx)*100,4)).',icol,'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
h = colorbar; ylabel(h,'relative (%)'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',icol), title('inverse DPF'); clim([min(icol) max(icol)])
set(gca,'fontsize',16), axis square, axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus)]); colormap jet
subplot(2,2,4)
contourf(set_of_mua,set_of_mus,squeeze(nanmean(set_of_dpfs_err(:,:,5,idx)*100,4)).',icol,'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
h = colorbar; ylabel(h,'relative (%)'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',icol), title('empirical DPF'); clim([min(icol) max(icol)])
set(gca,'fontsize',16), axis square, axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus)]); colormap jet
clearvars h idx icol

figure(4)
idx = 1<=d_cntrs&d_cntrs<=8;
subplot(2,2,1)
icol = -1:3:50;
contourf(set_of_mua,set_of_mus,squeeze(nanmean(set_of_dpfs_err(:,:,1,idx)*100,4)).',icol,'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
h = colorbar; ylabel(h,'relative (%)'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',[min(icol) mean(icol) max(icol)]), title('constant DPF'); clim([min(icol) max(icol)])
set(gca,'fontsize',16), axis square, axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus)]); colormap jet

subplot(2,2,2)
icol = 10:5:110;
contourf(set_of_mua,set_of_mus,squeeze(nanmean(set_of_dpfs_err(:,:,6,idx)*100,4)).','EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
h = colorbar; ylabel(h,'relative (%)'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',[min(icol) mean(icol) max(icol)]), title('empirical DPF'); clim([min(icol) max(icol)])
set(gca,'fontsize',16), axis square, axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus)]); colormap(gca,'jet')
clearvars icol set_of_lines my_colormap

subplot(2,2,3)
icol = -1.5:0.25:+2;
contourf(set_of_mua,set_of_mus,squeeze(nanmean(set_of_dpfs_err(:,:,4,idx)*100,4)).',icol,'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
h = colorbar; ylabel(h,'relative (%)'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',[min(icol) mean(icol) max(icol)]), title('inverse DPF'); clim([min(icol) max(icol)])
set(gca,'fontsize',16), axis square, axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus)]); colormap jet
subplot(2,2,4)
icol = -12:1:+4;
contourf(set_of_mua,set_of_mus,squeeze(nanmean(set_of_dpfs_err(:,:,5,idx)*100,4)).',icol,'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
h = colorbar; ylabel(h,'relative (%)'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',[min(icol) mean(icol) max(icol)]), title('empirical DPF'); clim([min(icol) max(icol)])
set(gca,'fontsize',16), axis square, axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus)]); colormap jet
clearvars h idx icol

figure(5)
idx = 1<=d_cntrs&d_cntrs<=8;
subplot(2,2,1)
icol = -1:3:50; [set_of_lines,my_colormap] = make_colormap(min(icol),max(icol),length(icol));
contourf(set_of_mua,set_of_mus,squeeze(nanmean(set_of_dpfs_err(:,:,1,idx)*100,4)).',set_of_lines,'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
h = colorbar; ylabel(h,'relative (%)'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',[min(icol) mean(icol) max(icol)]), title('constant DPF'); clim([min(icol) max(icol)])
set(gca,'fontsize',16), axis square, axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus)]); colormap(gca,my_colormap)
clearvars icol set_of_lines my_colormap
subplot(2,2,3)
icol = -1.5:0.25:+2; icol = -12:1:+4; [set_of_lines,my_colormap] = make_colormap(min(icol),max(icol),length(icol));
contourf(set_of_mua,set_of_mus,squeeze(nanmean(set_of_dpfs_err(:,:,4,idx)*100,4)).',set_of_lines,'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
h = colorbar; ylabel(h,'relative (%)'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',[min(icol) mean(icol) max(icol)]), title('inverse DPF'); clim([min(icol) max(icol)])
set(gca,'fontsize',16), axis square, axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus)]); colormap(gca,my_colormap)
clearvars icol set_of_lines my_colormap
subplot(2,2,4)
icol = -12:1:+4; [set_of_lines,my_colormap] = make_colormap(min(icol),max(icol),length(icol));
contourf(set_of_mua,set_of_mus,squeeze(nanmean(set_of_dpfs_err(:,:,5,idx)*100,4)).',set_of_lines,'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
h = colorbar; ylabel(h,'relative (%)'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',[min(icol) mean(icol) max(icol)]), title('empirical DPF'); clim([min(icol) max(icol)])
set(gca,'fontsize',16), axis square, axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus)]); colormap(gca,my_colormap)
clearvars icol set_of_lines my_colormap
clearvars h idx icol





figure(6)
for i_a = 1:length(set_of_mua)
    for i_s = 1:length(set_of_mus)
        plot(d_cntrs,squeeze(set_of________I(i_a,i_s,:)),'Color',get_color(set_of_mua(i_a)/max(set_of_mua),set_of_mus(i_s)/max(set_of_mus)),'DisplayName',['\mu_a = ',num2str(set_of_mua(i_a)),', \mu_s = ',num2str(set_of_mus(i_s))]), hold on
    end
end
clearvars i_a i_s

figure(7)
for i_a = 1:length(set_of_mua)
    for i_s = 1:length(set_of_mus)
        plot(d_cntrs,squeeze(set_of_dpfs_err(i_a,i_s,1,:)),'Color',get_color(set_of_mua(i_a)/max(set_of_mua),set_of_mus(i_s)/max(set_of_mus)),'DisplayName',['\mu_a = ',num2str(set_of_mua(i_a)),', \mu_s = ',num2str(set_of_mus(i_s))]), hold on
        % plot(d_cntrs,squeeze(set_of_dpfs_err(i_a,i_s,2,:)),'Color',get_color(set_of_mua(i_a)/max(set_of_mua),set_of_mus(i_s)/max(set_of_mus)),'DisplayName',['\mu_a = ',num2str(set_of_mua(i_a)),', \mu_s = ',num2str(set_of_mus(i_s))]), hold on
        plot(d_cntrs,squeeze(set_of_dpfs_err(i_a,i_s,4,:)),'Color',get_color(set_of_mua(i_a)/max(set_of_mua),set_of_mus(i_s)/max(set_of_mus)),'DisplayName',['\mu_a = ',num2str(set_of_mua(i_a)),', \mu_s = ',num2str(set_of_mus(i_s))]), hold on
    end
end
clearvars i_a i_s

figure(8)
subplot(2,1,1), mesh(set_of_mua,set_of_mus,squeeze(set_of_idst_DPF(:,:,1)).','FaceColor','none','EdgeColor','b','DisplayName','fitted'), hold on
subplot(2,1,1), mesh(set_of_mua,set_of_mus,p_est.',                          'FaceColor','none','EdgeColor','r','DisplayName','empirical'), hold on
subplot(2,1,1), mesh(set_of_mua,set_of_mus,p_anc.',                          'FaceColor','none','EdgeColor','g','DisplayName','analytical'), hold on
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
zlabel('p (a.u.)'), set(gca,'ztick',0:25:75), title('p'), hold off; legend
set(gca,'fontsize',16), axis square, axis tight, colormap jet; view([145 22.5])
axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus) 0 75])
subplot(2,1,2), mesh(set_of_mua,set_of_mus,squeeze(set_of_idst_DPF(:,:,2)).','FaceColor','none','EdgeColor','b','DisplayName','fitted'), hold on
subplot(2,1,2), mesh(set_of_mua,set_of_mus,q_est.',                          'FaceColor','none','EdgeColor','r','DisplayName','empirical'), hold on
subplot(2,1,2), mesh(set_of_mua,set_of_mus,q_anc.',                          'FaceColor','none','EdgeColor','g','DisplayName','analytical'), hold on
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
zlabel('b (a.u.)'), set(gca,'ztick',0:25:75), title('q'), hold off; legend
set(gca,'fontsize',16), axis square, axis tight, colormap jet; view([145 22.5])
axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus) 0 75])
clearvars p_est q_est



figure(9)
idx = 1<=d_cntrs&d_cntrs<=8;
subplot(2,2,1)
icol = -1:3:50; icol = -12:1:+12; [set_of_lines,my_colormap] = make_colormap(min(icol),max(icol),length(icol));
contourf(set_of_mua,set_of_mus,squeeze(nanmean(set_of_dpfs_err(:,:,1,idx)*100,4)).',set_of_lines,'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
h = colorbar; ylabel(h,'relative (%)'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',[min(icol):6:max(icol)]), title('constant DPF'); clim([min(icol) max(icol)])
set(gca,'fontsize',16), axis square, axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus)]); colormap(gca,my_colormap)
clearvars icol set_of_lines my_colormap
subplot(2,2,3)
icol = -1.5:0.25:+2; icol = -12:1:+12; [set_of_lines,my_colormap] = make_colormap(min(icol),max(icol),length(icol));
contourf(set_of_mua,set_of_mus,squeeze(nanmean(set_of_dpfs_err(:,:,4,idx)*100,4)).',set_of_lines,'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
h = colorbar; ylabel(h,'relative (%)'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',[min(icol):6:max(icol)]), title('inverse DPF'); clim([min(icol) max(icol)])
set(gca,'fontsize',16), axis square, axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus)]); colormap(gca,my_colormap)
clearvars icol set_of_lines my_colormap
subplot(2,2,4)
icol = -12:1:+12; [set_of_lines,my_colormap] = make_colormap(min(icol),max(icol),length(icol));
contourf(set_of_mua,set_of_mus,squeeze(nanmean(set_of_dpfs_err(:,:,5,idx)*100,4)).',set_of_lines,'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
h = colorbar; ylabel(h,'relative (%)'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',[min(icol):6:max(icol)]), title('empirical DPF'); clim([min(icol) max(icol)])
set(gca,'fontsize',16), axis square, axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus)]); colormap(gca,my_colormap)
clearvars icol set_of_lines my_colormap
clearvars h idx icol



figure(10)
idx = 1<=d_cntrs&d_cntrs<=8;
subplot(2,2,1)
icol = -1:3:50; icol = -12:1:+12; [set_of_lines,my_colormap] = make_colormap(min(icol),max(icol),length(icol));
contourf(set_of_mua,set_of_mus,squeeze(nanmean(set_of_dpfs_err(:,:,1,idx)*100,4)).',set_of_lines,'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
h = colorbar; ylabel(h,'relative (%)'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',[min(icol):6:max(icol)]), title('constant DPF'); clim([min(icol) max(icol)])
set(gca,'fontsize',16), axis square, axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus)]); colormap(gca,'jet')
clearvars icol set_of_lines my_colormap
subplot(2,2,3)
icol = -1.5:0.25:+2; icol = -12:1:+12; [set_of_lines,my_colormap] = make_colormap(min(icol),max(icol),length(icol));
contourf(set_of_mua,set_of_mus,squeeze(nanmean(set_of_dpfs_err(:,:,4,idx)*100,4)).',set_of_lines,'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
h = colorbar; ylabel(h,'relative (%)'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',[min(icol):6:max(icol)]), title('inverse DPF'); clim([min(icol) max(icol)])
set(gca,'fontsize',16), axis square, axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus)]); colormap(gca,'jet')
clearvars icol set_of_lines my_colormap
subplot(2,2,4)
icol = -12:1:+12; [set_of_lines,my_colormap] = make_colormap(min(icol),max(icol),length(icol));
contourf(set_of_mua,set_of_mus,squeeze(nanmean(set_of_dpfs_err(:,:,5,idx)*100,4)).',set_of_lines,'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
h = colorbar; ylabel(h,'relative (%)'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',[min(icol):6:max(icol)]), title('empirical DPF'); clim([min(icol) max(icol)])
set(gca,'fontsize',16), axis square, axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus)]); colormap(gca,'jet')
clearvars icol set_of_lines my_colormap
clearvars h idx icol

% save('Photon_33_5.mat','d_cntrs','d_edges','set_of_mua','set_of_mus','set_of________I','set_of_cnst_DPF','set_of_savg_DPF','set_of_idst_DPF','set_of_dpfs_err','set_of_true_DPF')
disp(['----------------------------------------------------------------------------------------------'])
disp(['----------------------------------------------------------------------------------------------'])
disp(['----------------------------------------------------------------------------------------------'])
disp(['# of photons = ',num2str(mean(number_of_all_photons)),' ± ',num2str(std(number_of_all_photons))])
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
function [set_of_lines,my_colormap] = make_colormap(l_start,l_end,N_lines)
set_of_lines = linspace(l_start,l_end,2*N_lines+1);
my_colormap = [1, 1, 1];
if     l_end<=0
    for idx = 1:2*N_lines
        my_colormap = [...
            [1-idx/2/N_lines, 1-idx/2/N_lines, 1];...
            my_colormap];
    end
elseif l_start>=0
    for idx = 1:2*N_lines
        my_colormap = [...
            my_colormap;...
            [1, 1-idx/2/N_lines, 1-idx/2/N_lines]];
    end
else
    Np = sum(set_of_lines>0);
    Nn = sum(set_of_lines<0);
    for idx = 1:Np
        my_colormap = [...
            my_colormap;...
            [1, 1-idx/Np, 1-idx/Np]];
    end
    for idx = 1:Nn
        my_colormap = [...
            [1-idx/Nn, 1-idx/Nn, 1];...
            my_colormap];
    end
end
end
