function [] = Photon_27_1 ()
% Repository group: 04_compression_beam_and_dpf
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: s vs. d for breast tissue (compression-induced hemodynamic changes in breast tissue / analysis)
% ref. #1: Assessment of the size, position, and optical properties of breast tumors in vivo by noninvasive optical methods
% ref. #2: Quantitative Absorption and Scattering Spectra in Thick Tissues Using Broadband Diffuse Optical Spectroscopy
% ref. #3: Tumor location of the central and nipple portion is associated with impaired survival for women with breast cancer
% ref. #4: Broadband Optical Mammography: Chromophore Concentration and Hemoglobin Saturation Contrast in Breast Cancer

clc
close all

z_air = 0.0; % the thickness of air layer
wvlnt_sample = 800; % sample wavelength (nm)
nPhotonsReq = 1e5; % number of requested photon
nExamplePaths = 50000; % number of example photon paths
% Optical properties of breast and tumor [1-4]
[~,~,~,~,brst_strct,tumr_strct] = give_breast_mua_mus(wvlnt_sample,false); clearvars wvlnt_sample
Lx = 19.1; Ly = 19.1; Lz = 10.0; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
cx = 0.00; cy = 0.00; cz = 1.50; % Tumor depth from the skin (in cm)(0.5 to 5.0 cm)
rx = 1.00; ry = 1.00; rz = 1.00; % Tumor radius within the breast (in cm)(1.0 to 2.5 cm).
do_count = 0;

for idx = 1
    set_of_scal = 1:0.5:3; % set of scales
    for i_scale = 1:length(set_of_scal)
        do_count = do_count + 1;
        mua_brst = brst_strct(idx).mua; g_brst = brst_strct(idx).g; n_brst = brst_strct(idx).n; mus_brst = brst_strct(idx).mus/(1-g_brst);
        mua_tumr = tumr_strct(idx).mua; g_tumr = tumr_strct(idx).g; n_tumr = tumr_strct(idx).n; mus_tumr = tumr_strct(idx).mus/(1-g_tumr); wvlnt = tumr_strct(idx).lambda;
        db = load(['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_27_',sprintf('%.1f',set_of_scal(i_scale)),'.mat']);
        d = db.d; s = db.s; w = db.w; a = db.a; no_of_surf_photons = db.no_of_surf_photons; no_of_dpth_photons = db.no_of_dpth_photons; M_raw = db.M_raw;

        N_bin = 150;
        s_edges = linspace(min(s),max(s),N_bin+1); % # s bins
        s_c = 1/2*(s_edges(1:end-1)+s_edges(2:end-0))';
        [f_s,~,ind_s] = histcounts(s,s_edges); f_s = f_s./sum(f_s); clearvars s_edges
        d_edges = linspace(min(d),max(d),N_bin+1); % # d bins
        d_c = 1/2*(d_edges(1:end-1)+d_edges(2:end-0))';
        [f_d,~,ind_d] = histcounts(d,d_edges); f_d = f_d./sum(f_d); clearvars d_edges

        w_mean = accumarray(ind_d,w,[],@sum); w_mean = w_mean/(no_of_surf_photons+no_of_dpth_photons);
        s_mean = accumarray(ind_d,s,[],@mean);
        s_var  = accumarray(ind_d,s,[],@var);
        s_var(s_mean==0) = nan; s_mean(s_mean==0) = nan;
        clearvars ind_s ind_d

        figure(2),
        imagesc(linspace(-Ly/2,+Ly/2,size(M_raw,2)),-linspace(0,Lz,size(M_raw,3)),squeeze(M_raw(round(size(M_raw,1)/2),:,:)).'),
        xlabel('y (cm)'), ylabel('z (cm)'), title(['model for scale = ',num2str(set_of_scal(i_scale))])
        set(gca,'fontsize',24), axis xy, axis equal, axis tight, grid on, hold off
        saveas(gcf,['Photon_27_',sprintf('%.1f',set_of_scal(i_scale)),'_side.tiff'],'tif')

        figure(3),
        imagesc(linspace(-Lx/2,+Lx/2,size(M_raw,1)),linspace(-Ly/2,+Ly/2,size(M_raw,2)),squeeze(M_raw(:,:,round(cz/set_of_scal(i_scale)/0.1))).'),
        xlabel('x (cm)'), ylabel('y (cm)'), title(['model for scale = ',num2str(set_of_scal(i_scale))])
        set(gca,'fontsize',24), axis xy, axis equal, axis tight, grid on, hold off
        saveas(gcf,['Photon_27_',sprintf('%.1f',set_of_scal(i_scale)),'_top.tiff'],'tif')

        figure(4)
        [N_w,edges_w] = histcounts(w,0:.01:1,"Normalization","probability");
        plot(1/2*(edges_w(1:end-1)+edges_w(2:end-0)),N_w,'Color',get_color(do_count),'DisplayName',['scale = ',num2str(set_of_scal(i_scale))]), hold on
        xlabel('w (a.u.)'), ylabel('freq (a.u.)'), title('hist. of w')
        set(gca,'fontsize',24), axis tight, grid on, legend show, hold on
        clearvars N_w edges_w

        figure(5)
        [N_s,edges_s] = histcounts(s,0:.5:100,"Normalization","probability");
        plot(1/2*(edges_s(1:end-1)+edges_s(2:end-0)),N_s,'Color',get_color(do_count),'DisplayName',['scale = ',num2str(set_of_scal(i_scale))]), hold on
        xlabel('s (cm)'), ylabel('freq (a.u.)'), title('hist. of s')
        set(gca,'fontsize',24), axis tight, grid on, legend show, hold on
        clearvars N_s edges_s

        clearvars d s a w d_c s_c s_mean w_mean i_mean s_var f_d f_s no_of_dpth_photons no_of_surf_photons
        clearvars mua_brst g_brst n_brst mus_brst
        clearvars mua_tumr g_tumr n_tumr mus_tumr
    end
end
end
function [out] = get_color(idx)
if     idx==1
    out = [0.0000 0.4470 0.7410];
elseif idx==2
    out = [0.8500 0.3250 0.0980];
elseif idx==3
    out = [0.9290 0.6940 0.1250];
elseif idx==4
    out = [0.4940 0.1840 0.5560];
elseif idx==5
    out = [0.4660 0.6740 0.1880];
elseif idx==6
    out = [0.3010 0.7450 0.9330];
elseif idx==7
    out = [0.6350 0.0780 0.1840];
else
    out = [0.0000 0.0000 0.0000];
end
end
function [mua_brst,mua_tumr,mus_brst,mus_tumr,brst_strct,tumr_strct] = give_breast_mua_mus(lambda_in,do_fig)
% data
tumr_db = readtable('DB/breast.tumor.csv');
brst_db = readtable('DB/breast.normal.csv');
brst_strct(1).lambda = 690; brst_strct(1).mua = 0.030; brst_strct(1).mus = 12.0; brst_strct(1).g = 0.95; brst_strct(1).n = 1.4;
brst_strct(2).lambda = 825; brst_strct(2).mua = 0.040; brst_strct(2).mus = 11.0; brst_strct(2).g = 0.95; brst_strct(2).n = 1.4;
tumr_strct(1).lambda = 690; tumr_strct(1).mua = 0.084; tumr_strct(1).mus = 15.0; tumr_strct(1).g = 0.95; tumr_strct(1).n = 1.4;
tumr_strct(2).lambda = 825; tumr_strct(2).mua = 0.085; tumr_strct(2).mus = 12.7; tumr_strct(2).g = 0.95; tumr_strct(2).n = 1.4;
% fit
bs_brst = -log(brst_strct(2).mus/brst_strct(1).mus)/log(brst_strct(2).lambda/brst_strct(1).lambda); as_brst = brst_strct(1).mus/brst_strct(1).lambda^(-bs_brst); mus_brst = as_brst*(lambda_in)^(-bs_brst); % clearvars as_brst bs_brst
bs_tumr = -log(tumr_strct(2).mus/tumr_strct(1).mus)/log(tumr_strct(2).lambda/tumr_strct(1).lambda); as_tumr = tumr_strct(1).mus/tumr_strct(1).lambda^(-bs_tumr); mus_tumr = as_tumr*(lambda_in)^(-bs_tumr); % clearvars as_tumr bs_tumr
[x,iv,~] = unique(brst_db.Var1); v = brst_db.Var2(iv)*10; mua_brst = interp1(x,v,lambda_in,'pchip','extrap'); clearvars x v iv
[x,iv,~] = unique(tumr_db.Var1); v = tumr_db.Var2(iv)*10; mua_tumr = interp1(x,v,lambda_in,'pchip','extrap'); clearvars x v iv
if do_fig
    subplot(1,2,1)
    plot(brst_db.Var1,brst_db.Var2*10,'DisplayName','breast: model','LineStyle','-','Color','b','LineWidth',1.5), hold on
    plot([brst_strct(:).lambda],[brst_strct(:).mua],'Marker','s','MarkerFaceColor','b','MarkerEdgeColor','k','DisplayName','breast: data','MarkerSize',16,'LineStyle','none'), hold on
    plot(lambda_in,mua_brst,'Marker','h','MarkerFaceColor','b','MarkerEdgeColor','k','DisplayName','breast: output','MarkerSize',16,'LineStyle','none'), hold on
    plot(tumr_db.Var1,tumr_db.Var2*10,'DisplayName','tumor: model', 'LineStyle','-','Color','r','LineWidth',1.5), hold on
    plot([tumr_strct(:).lambda],[tumr_strct(:).mua],'Marker','s','MarkerFaceColor','r','MarkerEdgeColor','k','DisplayName','tumor: data','MarkerSize',16,'LineStyle','none'), hold on
    plot(lambda_in,mua_tumr,'Marker','h','MarkerFaceColor','r','MarkerEdgeColor','k','DisplayName','tumor: output','MarkerSize',16,'LineStyle','none'), hold on
    xlabel('wavelength (\lambda, nm)'), ylabel('\mu_a (cm^{-1})'), title('\mu_a vs. wavelength')
    axis tight, axis square, set(gca,'fontsize',18)
    legend('show','Location','northwest')
    hold off

    subplot(1,2,2)
    plot(brst_db.Var1,as_brst.*(brst_db.Var1.^(-bs_brst)),'DisplayName','breast: model','LineStyle','-','Color','b','LineWidth',1.5), hold on
    plot([brst_strct(:).lambda],[brst_strct(:).mus],'Marker','s','MarkerFaceColor','b','MarkerEdgeColor','k','DisplayName','breast: data','MarkerSize',16,'LineStyle','none'), hold on
    plot(lambda_in,mus_brst,'Marker','h','MarkerFaceColor','b','MarkerEdgeColor','k','DisplayName','breast: output','MarkerSize',16,'LineStyle','none'), hold on
    plot(tumr_db.Var1,as_tumr.*(tumr_db.Var1.^(-bs_tumr)),'DisplayName','tumor: model', 'LineStyle','-','Color','r','LineWidth',1.5), hold on
    plot([tumr_strct(:).lambda],[tumr_strct(:).mus],'Marker','s','MarkerFaceColor','r','MarkerEdgeColor','k','DisplayName','tumor: data','MarkerSize',16,'LineStyle','none'), hold on
    plot(lambda_in,mus_tumr,'Marker','h','MarkerFaceColor','r','MarkerEdgeColor','k','DisplayName','tumor: output','MarkerSize',16,'LineStyle','none'), hold on
    xlabel('wavelength (\lambda, nm)'), ylabel('\mu''_s (cm^{-1})'), title('\mu''_s vs. wavelength')
    axis tight, axis square, set(gca,'fontsize',18)
    legend('show','Location','northeast')
    hold off
end
end
