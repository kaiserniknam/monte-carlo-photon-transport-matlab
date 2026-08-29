function [] = Photon_45_5 ()
% Repository group: 05_digital_breast_phantoms
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: simulating phantoms prepared by Diego
% Similar to Photon_45_3: but in better intervals
% analyzing simulated data: Testing an idea

clc
close all
load('DB/Photon_45_4_151.mat');
code_num = 45;
p_threshold = 0.05;
wvlnt = 800; % sample wavelength (nm)
% Optical properties of breast and tumor [1-4]
mua_glnd = 0.55; mus_glnd = 332.7; g_glnd = 0.965; n_glnd = 1.4;
mua_adps = 0.82; mus_adps = 313.8; g_adps = 0.976; n_adps = 1.4;
mua_aire = 1e-4; mus_aire = 1.e-5; g_aire = 1.000; n_aire = 1.0;
% calculate optical properties of tumor
set_of_mua = nan(size(set_of_HCT_Cn));
set_of_mus = nan(size(set_of_HCT_Cn));
set_of___g = nan(size(set_of_HCT_Cn));
set_of___n = nan(size(set_of_HCT_Cn));
for counter = 1:length(set_of_HCT_Cn)
    [set_of_mua(counter),set_of_mus(counter),set_of___g(counter),~] = calc_muas_based_HCT(set_of_HCT_Cn(counter),wvlnt);
    set_of___n(counter) = refractive_index_water(wvlnt).*((Specific_Refractive_Increment_beta(wvlnt).*set_of_HCT_Cn(counter)./3)+1);
end
d_centers = 1/2*(d_edges(1:end-1)+d_edges(2:end-0));
clearvars counter d_edges

% 1-3: D. A. Boas, C. Pitris, and N. Ramanujam, Handbook of biomedical optics. CRC press, 2016.
lit_data = [
    0.28 548.8 0.981 1.4;
    0.34 384.0 0.983 1.4;
    0.34 233.0 0.962 1.4
    ];

figure(1), subplot(2,2,1) % mua
axis([0 max(set_of_HCT_Cn)+2.5 -0.5 3.5])
p_start = 6.5; p_finish = 11.5; Y_lim = ylim;
rectangle('Position',[p_start min(Y_lim) p_finish-p_start max(Y_lim)-min(Y_lim)],'FaceColor','g','EdgeColor','none','FaceAlpha',0.25), hold on
p_start = 15; p_finish = 20; Y_lim = ylim;
rectangle('Position',[p_start min(Y_lim) p_finish-p_start max(Y_lim)-min(Y_lim)],'FaceColor','y','EdgeColor','none','FaceAlpha',0.25), hold on
plot(set_of_HCT_Cn,set_of_mua,'DisplayName','tumor','LineWidth',2,'Color',get_color(4)), hold on
plot(set_of_HCT_Cn,sign(set_of_mua).*mua_aire,'DisplayName','air','LineWidth',2,'Color',get_color(1)), hold on
plot(set_of_HCT_Cn,sign(set_of_mua).*mua_adps,'DisplayName','adipose/fat','LineWidth',2,'Color',get_color(3)), hold on
plot(set_of_HCT_Cn,sign(set_of_mua).*mua_glnd,'DisplayName','fibroglandular','LineWidth',2,'Color',get_color(2)), hold on
plot(set_of_HCT_Cn,ones(size(set_of_mua)).*lit_data(1,1),'DisplayName','data #1','LineStyle','none','Marker','s','MarkerFaceColor',get_color(5),'MarkerSize',12), hold on
plot(set_of_HCT_Cn,ones(size(set_of_mua)).*lit_data(2,1),'DisplayName','data #2','LineStyle','none','Marker','s','MarkerFaceColor',get_color(6),'MarkerSize',12), hold on
plot(set_of_HCT_Cn,ones(size(set_of_mua)).*lit_data(3,1),'DisplayName','data #3','LineStyle','none','Marker','s','MarkerFaceColor',get_color(7),'MarkerSize',12), hold on
xlabel('HCT (%)'), ylabel('\mu_a (cm^{-1})'), title('\mu_a of tissues'), set(gca,'fontsize',16), hold on
legend('show','Location','northwest')
figure(1), subplot(2,2,2) % mus
axis([0 max(set_of_HCT_Cn)+2.5 -50 850])
p_start = 6.5; p_finish = 11.5; Y_lim = ylim;
rectangle('Position',[p_start min(Y_lim) p_finish-p_start max(Y_lim)-min(Y_lim)],'FaceColor','g','EdgeColor','none','FaceAlpha',0.25), hold on
p_start = 15; p_finish = 20; Y_lim = ylim;
rectangle('Position',[p_start min(Y_lim) p_finish-p_start max(Y_lim)-min(Y_lim)],'FaceColor','y','EdgeColor','none','FaceAlpha',0.25), hold on
plot(set_of_HCT_Cn,set_of_mus,'DisplayName','tumor','LineWidth',2,'Color',get_color(4)), hold on
plot(set_of_HCT_Cn,sign(set_of_mus).*mus_aire,'DisplayName','air','LineWidth',2,'Color',get_color(1)), hold on
plot(set_of_HCT_Cn,sign(set_of_mus).*mus_adps,'DisplayName','adipose/fat','LineWidth',2,'Color',get_color(3)), hold on
plot(set_of_HCT_Cn,sign(set_of_mus).*mus_glnd,'DisplayName','fibroglandular','LineWidth',2,'Color',get_color(2)), hold on
plot(set_of_HCT_Cn,ones(size(set_of_mua)).*lit_data(1,2),'DisplayName','data #1','LineStyle','none','Marker','s','MarkerFaceColor',get_color(5),'MarkerSize',12), hold on
plot(set_of_HCT_Cn,ones(size(set_of_mua)).*lit_data(2,2),'DisplayName','data #2','LineStyle','none','Marker','s','MarkerFaceColor',get_color(6),'MarkerSize',12), hold on
plot(set_of_HCT_Cn,ones(size(set_of_mua)).*lit_data(3,2),'DisplayName','data #3','LineStyle','none','Marker','s','MarkerFaceColor',get_color(7),'MarkerSize',12), hold on
xlabel('HCT (%)'), ylabel('\mu_s (cm^{-1})'), title('\mu_s of tissues'), set(gca,'fontsize',16), hold on
legend('show','Location','northwest')
figure(1), subplot(2,2,3) % g
axis([0 max(set_of_HCT_Cn)+2.5 0.95 1.01])
p_start = 6.5; p_finish = 11.5; Y_lim = ylim;
rectangle('Position',[p_start min(Y_lim) p_finish-p_start max(Y_lim)-min(Y_lim)],'FaceColor','g','EdgeColor','none','FaceAlpha',0.25), hold on
p_start = 15; p_finish = 20; Y_lim = ylim;
rectangle('Position',[p_start min(Y_lim) p_finish-p_start max(Y_lim)-min(Y_lim)],'FaceColor','y','EdgeColor','none','FaceAlpha',0.25), hold on
plot(set_of_HCT_Cn,set_of___g,'DisplayName','tumor','LineWidth',2,'Color',get_color(4)), hold on
plot(set_of_HCT_Cn,sign(set_of___g).*g_aire,'DisplayName','air','LineWidth',2,'Color',get_color(1)), hold on
plot(set_of_HCT_Cn,sign(set_of___g).*g_adps,'DisplayName','adipose/fat','LineWidth',2,'Color',get_color(3)), hold on
plot(set_of_HCT_Cn,sign(set_of___g).*g_glnd,'DisplayName','fibroglandular','LineWidth',2,'Color',get_color(2)), hold on
plot(set_of_HCT_Cn,ones(size(set_of_mua)).*lit_data(1,3),'DisplayName','data #1','LineStyle','none','Marker','s','MarkerFaceColor',get_color(5),'MarkerSize',12), hold on
plot(set_of_HCT_Cn,ones(size(set_of_mua)).*lit_data(2,3),'DisplayName','data #2','LineStyle','none','Marker','s','MarkerFaceColor',get_color(6),'MarkerSize',12), hold on
plot(set_of_HCT_Cn,ones(size(set_of_mua)).*lit_data(3,3),'DisplayName','data #3','LineStyle','none','Marker','s','MarkerFaceColor',get_color(7),'MarkerSize',12), hold on
xlabel('HCT (%)'), ylabel('g (a.u.)'), title('g of tissues'), set(gca,'fontsize',16), hold on
legend('show','Location','southwest','NumColumns',2)
figure(1), subplot(2,2,4) % n
axis([0 max(set_of_HCT_Cn)+2.5 0.95 1.45])
p_start = 6.5; p_finish = 11.5; Y_lim = ylim;
rectangle('Position',[p_start min(Y_lim) p_finish-p_start max(Y_lim)-min(Y_lim)],'FaceColor','g','EdgeColor','none','FaceAlpha',0.25), hold on
p_start = 15; p_finish = 20; Y_lim = ylim;
rectangle('Position',[p_start min(Y_lim) p_finish-p_start max(Y_lim)-min(Y_lim)],'FaceColor','y','EdgeColor','none','FaceAlpha',0.25), hold on
plot(set_of_HCT_Cn,set_of___n.*sign(set_of___g),'DisplayName','tumor','LineWidth',2,'Color',get_color(4)), hold on
plot(set_of_HCT_Cn,sign(set_of___g).*n_aire,'DisplayName','air','LineWidth',2,'Color',get_color(1)), hold on
plot(set_of_HCT_Cn,sign(set_of___g).*n_adps,'DisplayName','adipose/fat','LineWidth',2,'Color',get_color(3)), hold on
plot(set_of_HCT_Cn,sign(set_of___g).*n_glnd,'DisplayName','fibroglandular','LineWidth',2,'Color',get_color(2)), hold on
plot(set_of_HCT_Cn,ones(size(set_of_mua)).*lit_data(1,4),'DisplayName','data #1','LineStyle','none','Marker','s','MarkerFaceColor',get_color(5),'MarkerSize',12), hold on
plot(set_of_HCT_Cn,ones(size(set_of_mua)).*lit_data(2,4),'DisplayName','data #2','LineStyle','none','Marker','s','MarkerFaceColor',get_color(6),'MarkerSize',12), hold on
plot(set_of_HCT_Cn,ones(size(set_of_mua)).*lit_data(3,4),'DisplayName','data #3','LineStyle','none','Marker','s','MarkerFaceColor',get_color(7),'MarkerSize',12), hold on
xlabel('HCT (%)'), ylabel('n (a.u.)'), title('n of tissues'), set(gca,'fontsize',16), hold on
legend('show','Location','east','NumColumns',2)



% OD vs mua (several d)
idx_d = 1:4:39; Y_Limit = [2 17]; i_beam_X = 2; idx_HCT_Cn = 2:length(set_of_HCT_Cn);
for i_Pd = 1:3
    for i_d = 1:length(idx_d)% OD/mu/d vs mua (several d)
        TheColor = [(i_d-1)/(length(idx_d)-1),0,1-(i_d-1)/(length(idx_d)-1)];
        for i_Tumr_Z = 1:length(set_of_Tumr_Z)
            figure(i_Pd+1), subplot(2,2,i_Tumr_Z)
            p_comp = set_of_I_s(i_Pd,1:3,i_beam_X,i_Tumr_Z,idx_HCT_Cn,5,idx_d(i_d));
            y_refc = -log(set_of_I_s(i_Pd,1:3,i_beam_X,i_Tumr_Z,idx_HCT_Cn,4,idx_d(i_d)));
            y_temp = -log(set_of_I_s(i_Pd,1:3,i_beam_X,i_Tumr_Z,idx_HCT_Cn,2,idx_d(i_d)));
            y_temp(p_comp>p_threshold) = nan;

            plot(set_of_mua(idx_HCT_Cn),squeeze(nanmean(nanmean(y_refc,1),2)),'Color',TheColor,'HandleVisibility','off','LineWidth',0.5,'LineStyle','--'), hold on
            plot(set_of_mua(idx_HCT_Cn),squeeze(nanmean(nanmean(y_temp,1),2)),'Color',TheColor,'DisplayName',['d = ',num2str(d_centers(idx_d(i_d)))],'LineWidth',2,'Marker','o','MarkerFaceColor',TheColor,'MarkerEdgeColor','none','LineStyle','-','MarkerSize',6), hold on
            xlabel('\mu_a (cm^{-1})'), xlim([min(set_of_mua) max(set_of_mua)])
            ylabel('OD (a.u.)'), ylim([min(Y_Limit) max(Y_Limit)])
            set(gca,'fontsize',12), % legend('show','numColumns',1,'Orientation','horizontal'),
            title(['beamX = ',num2str(set_of_beam_X(i_beam_X)),' cm, Tumor depth = ',num2str(set_of_Tumr_Z(i_Tumr_Z)),' cm, Pd = ',num2str(set_of_percnt(i_Pd)),'%'])
            clearvars p_comp y_refc y_temp
        end
        clearvars TheColor
    end
end
clearvars i_Pd i_fig idx_d Y_Limit i_beam_X idx_HCT_Cn i_d i_Tumr_Z

i_fig = 5; idx_d = 1:4:39; Y_Limit = [2 17]; i_beam_X = 2; idx_HCT_Cn = 2:length(set_of_HCT_Cn);
for i_d = 1:length(idx_d)
    TheColor = [(i_d-1)/(length(idx_d)-1),0,1-(i_d-1)/(length(idx_d)-1)];
    for i_Tumr_Z = 1:length(set_of_Tumr_Z)
        figure(i_fig), subplot(2,2,i_Tumr_Z)
        p_comp = set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn,5,idx_d(i_d));
        y_refc = -log(set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn,4,idx_d(i_d)));
        y_temp = -log(set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn,2,idx_d(i_d)));
        y_temp(p_comp>p_threshold) = nan;

        plot(set_of_mua(idx_HCT_Cn),squeeze(nanmean(nanmean(y_refc,1),2)),'Color',TheColor,'HandleVisibility','off','LineWidth',0.5,'LineStyle','--'), hold on
        plot(set_of_mua(idx_HCT_Cn),squeeze(nanmean(nanmean(y_temp,1),2)),'Color',TheColor,'DisplayName',['d = ',num2str(d_centers(idx_d(i_d)))],'LineWidth',2,'Marker','o','MarkerFaceColor',TheColor,'MarkerEdgeColor','none','LineStyle','-','MarkerSize',6), hold on
        xlabel('\mu_a (cm^{-1})'), xlim([min(set_of_mua) max(set_of_mua)])
        ylabel('OD (a.u.)'), ylim([min(Y_Limit) max(Y_Limit)])
        set(gca,'fontsize',12), % legend('show','numColumns',1,'Orientation','horizontal'),
        title(['beamX = ',num2str(set_of_beam_X(i_beam_X)),' cm, Tumor depth = ',num2str(set_of_Tumr_Z(i_Tumr_Z)),' cm'])

        clearvars p_comp y_refc y_temp
    end
    clearvars TheColor
end
clearvars i_fig idx_d Y_Limit i_beam_X idx_HCT_Cn i_d i_Tumr_Z



% OD vs d (several mua)
idx_d = 1:1:39; Y_Limit = [2 17]; i_beam_X = 2; idx_HCT_Cn = 2:1:length(set_of_HCT_Cn);
for i_Pd = 1:3
    for i_HCT = 1:length(idx_HCT_Cn)
        TheColor = [(i_HCT-1)/(length(idx_HCT_Cn)-1),0,1-(i_HCT-1)/(length(idx_HCT_Cn)-1)];
        for i_Tumr_Z = 1:length(set_of_Tumr_Z)
            figure(i_Pd+5), subplot(2,2,i_Tumr_Z)
            p_comp =     (set_of_I_s(i_Pd,1:3,i_beam_X,i_Tumr_Z,idx_HCT_Cn(i_HCT),5,idx_d));
            y_refc = -log(set_of_I_s(i_Pd,1:3,i_beam_X,i_Tumr_Z,idx_HCT_Cn(i_HCT),4,idx_d));
            y_temp = -log(set_of_I_s(i_Pd,1:3,i_beam_X,i_Tumr_Z,idx_HCT_Cn(i_HCT),2,idx_d));
            y_temp(p_comp>p_threshold) = nan;

            plot(d_centers(idx_d),squeeze(nanmean(nanmean(y_refc,1),2)),'Color',TheColor,'HandleVisibility','off','LineWidth',0.5,'LineStyle','--'), hold on
            plot(d_centers(idx_d),squeeze(nanmean(nanmean(y_temp,1),2)),'Color',TheColor,'DisplayName',['\mu_a = ',num2str(set_of_mua(idx_HCT_Cn(i_HCT))),''],'LineWidth',2,'Marker','o','MarkerFaceColor',TheColor,'MarkerEdgeColor','none','LineStyle','-','MarkerSize',6), hold on
            xlabel('d (cm)'), xlim([min(d_centers) 3])
            ylabel('OD (a.u.)'), ylim([min(Y_Limit) max(Y_Limit)])
            set(gca,'fontsize',12), % legend('show','numColumns',1,'Orientation','horizontal'),
            title(['beamX = ',num2str(set_of_beam_X(i_beam_X)),' cm, Tumor depth = ',num2str(set_of_Tumr_Z(i_Tumr_Z)),' cm, Pd = ',num2str(set_of_percnt(i_Pd)),'%'])

            clearvars p_comp y_refc y_tempset_of_mua
        end
        clearvars TheColor
    end
end
clearvars i_Pd i_fig idx_d Y_Limit i_beam_X idx_HCT_Cn i_d i_Tumr_Z

i_fig = 9; idx_d = 1:1:39; Y_Limit = [2 17]; i_beam_X = 2; idx_HCT_Cn = 2:1:length(set_of_HCT_Cn);
for i_HCT = 1:length(idx_HCT_Cn)
    TheColor = [(i_HCT-1)/(length(idx_HCT_Cn)-1),0,1-(i_HCT-1)/(length(idx_HCT_Cn)-1)];
    for i_Tumr_Z = 1:length(set_of_Tumr_Z)
        figure(i_fig), subplot(2,2,i_Tumr_Z)
        p_comp =     (set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn(i_HCT),5,idx_d));
        y_refc = -log(set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn(i_HCT),4,idx_d));
        y_temp = -log(set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn(i_HCT),2,idx_d));
        y_temp(p_comp>p_threshold) = nan;

        plot(d_centers(idx_d),squeeze(nanmean(nanmean(y_refc,1),2)),'Color',TheColor,'HandleVisibility','off','LineWidth',0.5,'LineStyle','--'), hold on
        plot(d_centers(idx_d),squeeze(nanmean(nanmean(y_temp,1),2)),'Color',TheColor,'DisplayName',['\mu_a = ',num2str(set_of_mua(idx_HCT_Cn(i_HCT))),''],'LineWidth',2,'Marker','o','MarkerFaceColor',TheColor,'MarkerEdgeColor','none','LineStyle','-','MarkerSize',6), hold on
        xlabel('d (cm)'), xlim([min(d_centers) 3])
        ylabel('OD (a.u.)'), ylim([min(Y_Limit) max(Y_Limit)])
        set(gca,'fontsize',12), % legend('show','numColumns',1,'Orientation','horizontal'),
        title(['beamX = ',num2str(set_of_beam_X(i_beam_X)),' cm, Tumor depth = ',num2str(set_of_Tumr_Z(i_Tumr_Z)),' cm'])

        clearvars p_comp y_refc y_tempset_of_mua
    end
    clearvars TheColor
end
clearvars i_fig idx_d Y_Limit i_beam_X idx_HCT_Cn i_d i_Tumr_Z i_HCT






% 3-D map: OD, mua, d
Z_Limit = [2 17]; i_beam_X = 2;
for i_Pd = 1:3
    for i_Tumr_Z = 1:length(set_of_Tumr_Z)
        figure(i_Pd+9), subplot(2,2,i_Tumr_Z)
        p_comp =     (set_of_I_s(i_Pd,1:3,i_beam_X,i_Tumr_Z,:,5,:));
        z_refc = -log(set_of_I_s(i_Pd,1:3,i_beam_X,i_Tumr_Z,:,4,:));
        z_temp = -log(set_of_I_s(i_Pd,1:3,i_beam_X,i_Tumr_Z,:,2,:));
        z_temp(p_comp>p_threshold) = nan;

        mesh(set_of_mua,d_centers,squeeze(nanmean(z_refc,2)).','EdgeColor','k','DisplayName','No Tumor','LineStyle','-'), hold on
        mesh(set_of_mua,d_centers,squeeze(nanmean(z_temp,2)).','EdgeColor','r','DisplayName','with Tumor','LineStyle','-'), hold on
        xlabel('\mu_a (cm^{-1})'), xlim([min(set_of_mua) max(set_of_mua)])
        ylabel('d (cm)'), ylim([min(d_centers) 3.5])
        zlabel('OD (a.u.)'), zlim([min(Z_Limit) max(Z_Limit)])
        set(gca,'fontsize',12), legend('show','numColumns',1), view([-300 35])
        title(['beamX = ',num2str(set_of_beam_X(i_beam_X)),' cm, Tumor depth = ',num2str(set_of_Tumr_Z(i_Tumr_Z)),' cm, Pd = ',num2str(set_of_percnt(i_Pd)),'%'])
        clearvars p_comp z_refc z_temp
    end
end
clearvars i_Pd i_fig Z_Limit i_beam_X  i_Tumr_Z

Z_Limit = [2 17]; i_beam_X = 2;
for i_Tumr_Z = 1:length(set_of_Tumr_Z)
    figure(4+9), subplot(2,2,i_Tumr_Z)
    p_comp =     (set_of_I_s(1:3,1:3,i_beam_X,i_Tumr_Z,:,5,:));
    z_refc = -log(set_of_I_s(1:3,1:3,i_beam_X,i_Tumr_Z,:,4,:));
    z_temp = -log(set_of_I_s(1:3,1:3,i_beam_X,i_Tumr_Z,:,2,:));
    z_temp(p_comp>p_threshold) = nan;

    mesh(set_of_mua,d_centers,squeeze(nanmean(nanmean(z_refc,1),2)).','EdgeColor','k','DisplayName','No Tumor','LineStyle','-'), hold on
    mesh(set_of_mua,d_centers,squeeze(nanmean(nanmean(z_temp,1),2)).','EdgeColor','r','DisplayName','with Tumor','LineStyle','-'), hold on
        xlabel('\mu_a (cm^{-1})'), xlim([min(set_of_mua) max(set_of_mua)])
    ylabel('d (cm)'), ylim([min(d_centers) 3.5])
    zlabel('OD (a.u.)'), zlim([min(Z_Limit) max(Z_Limit)])
    set(gca,'fontsize',12), legend('show','numColumns',1), view([-300 35])
    title(['beamX = ',num2str(set_of_beam_X(i_beam_X)),' cm, Tumor depth = ',num2str(set_of_Tumr_Z(i_Tumr_Z)),' cm'])
    clearvars p_comp z_refc z_temp
end
clearvars i_Pd i_fig Z_Limit i_beam_X  i_Tumr_Z






% OD_relative vs mua (several d)
idx_d = 1:4:39; Y_Limit = [-4 +4]; i_beam_X = 2; idx_HCT_Cn = 2:length(set_of_HCT_Cn);
for i_Pd = 1:3
    for i_d = 1:length(idx_d)
        TheColor = [(i_d-1)/(length(idx_d)-1),0,1-(i_d-1)/(length(idx_d)-1)];
        for i_Tumr_Z = 1:length(set_of_Tumr_Z)
            figure(i_Pd+13), subplot(2,2,i_Tumr_Z)
            p_comp = set_of_I_s(i_Pd,1:3,i_beam_X,i_Tumr_Z,idx_HCT_Cn,5,idx_d(i_d));
            y_refc = (set_of_I_s(i_Pd,1:3,i_beam_X,i_Tumr_Z,idx_HCT_Cn,4,idx_d(i_d)));
            y_temp = (set_of_I_s(i_Pd,1:3,i_beam_X,i_Tumr_Z,idx_HCT_Cn,2,idx_d(i_d)));
            y_temp = -log(y_temp./y_refc);
            y_temp(p_comp>p_threshold) = nan;

            plot(set_of_mua(idx_HCT_Cn),zeros(size(set_of_mua(idx_HCT_Cn))),  'Color','k','HandleVisibility','off','LineWidth',0.5,'LineStyle','--'), hold on
            plot(set_of_mua(idx_HCT_Cn),squeeze(nanmean(nanmean(y_temp,1),2)),'Color',TheColor,'DisplayName',['d = ',num2str(d_centers(idx_d(i_d)))],'LineWidth',2,'Marker','o','MarkerFaceColor',TheColor,'MarkerEdgeColor','none','LineStyle','-','MarkerSize',6), hold on
            xlabel('\mu_a (cm^{-1})'), xlim([min(set_of_mua) max(set_of_mua)])
            ylabel('OD (a.u.)'), ylim([min(Y_Limit) max(Y_Limit)])
            set(gca,'fontsize',12), % legend('show','numColumns',1,'Orientation','horizontal'),
            title(['beamX = ',num2str(set_of_beam_X(i_beam_X)),' cm, Tumor depth = ',num2str(set_of_Tumr_Z(i_Tumr_Z)),' cm, Pd = ',num2str(set_of_percnt(i_Pd)),'%'])
            clearvars p_comp y_refc y_temp
        end
        clearvars TheColor
    end
end
clearvars i_Pd i_fig idx_d Y_Limit i_beam_X idx_HCT_Cn i_d i_Tumr_Z

i_fig = 17; idx_d = 1:4:39; Y_Limit = [-4 +4]; i_beam_X = 2; idx_HCT_Cn = 2:length(set_of_HCT_Cn);
for i_d = 1:length(idx_d)
    TheColor = [(i_d-1)/(length(idx_d)-1),0,1-(i_d-1)/(length(idx_d)-1)];
    for i_Tumr_Z = 1:length(set_of_Tumr_Z)
        figure(i_fig), subplot(2,2,i_Tumr_Z)
        p_comp = set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn,5,idx_d(i_d));
        y_refc = (set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn,4,idx_d(i_d)));
        y_temp = (set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn,2,idx_d(i_d)));
        y_temp = -log(y_temp./y_refc);
        y_temp(p_comp>p_threshold) = nan;

        plot(set_of_mua(idx_HCT_Cn),zeros(size(set_of_mua(idx_HCT_Cn))),  'Color','k','HandleVisibility','off','LineWidth',0.5,'LineStyle','--'), hold on
        plot(set_of_mua(idx_HCT_Cn),squeeze(nanmean(nanmean(y_temp,1),2)),'Color',TheColor,'DisplayName',['d = ',num2str(d_centers(idx_d(i_d)))],'LineWidth',2,'Marker','o','MarkerFaceColor',TheColor,'MarkerEdgeColor','none','LineStyle','-','MarkerSize',6), hold on
        xlabel('\mu_a (cm^{-1})'), xlim([min(set_of_mua) max(set_of_mua)])
        ylabel('OD (a.u.)'), ylim([min(Y_Limit) max(Y_Limit)])
        set(gca,'fontsize',12), % legend('show','numColumns',1,'Orientation','horizontal'),
        title(['beamX = ',num2str(set_of_beam_X(i_beam_X)),' cm, Tumor depth = ',num2str(set_of_Tumr_Z(i_Tumr_Z)),' cm'])

        clearvars p_comp y_refc y_temp
    end
    clearvars TheColor
end
clearvars i_fig idx_d Y_Limit i_beam_X idx_HCT_Cn i_d i_Tumr_Z



% OD_relative vs d (several mua)
idx_d = 1:1:39; Y_Limit = [-4 +4]; i_beam_X = 2; idx_HCT_Cn = 2:1:length(set_of_HCT_Cn);
for i_Pd = 1:3
    for i_HCT = 1:length(idx_HCT_Cn)
        TheColor = [(i_HCT-1)/(length(idx_HCT_Cn)-1),0,1-(i_HCT-1)/(length(idx_HCT_Cn)-1)];
        for i_Tumr_Z = 1:length(set_of_Tumr_Z)
            figure(i_Pd+17), subplot(2,2,i_Tumr_Z)
            p_comp =     (set_of_I_s(i_Pd,1:3,i_beam_X,i_Tumr_Z,idx_HCT_Cn(i_HCT),5,idx_d));
            y_refc = (set_of_I_s(i_Pd,1:3,i_beam_X,i_Tumr_Z,idx_HCT_Cn(i_HCT),4,idx_d));
            y_temp = (set_of_I_s(i_Pd,1:3,i_beam_X,i_Tumr_Z,idx_HCT_Cn(i_HCT),2,idx_d));
            y_temp = -log(y_temp./y_refc);
            y_temp(p_comp>p_threshold) = nan;

            plot(set_of_mua(idx_HCT_Cn),zeros(size(set_of_mua(idx_HCT_Cn))),'Color','k','HandleVisibility','off','LineWidth',0.5,'LineStyle','--'), hold on
            plot(d_centers(idx_d),squeeze(nanmean(nanmean(y_temp,1),2)),    'Color',TheColor,'DisplayName',['\mu_a = ',num2str(set_of_mua(idx_HCT_Cn(i_HCT))),''],'LineWidth',2,'Marker','o','MarkerFaceColor',TheColor,'MarkerEdgeColor','none','LineStyle','-','MarkerSize',6), hold on
            xlabel('d (cm)'), xlim([min(d_centers) 3])
            ylabel('OD (a.u.)'), ylim([min(Y_Limit) max(Y_Limit)])
            set(gca,'fontsize',12), % legend('show','numColumns',1,'Orientation','horizontal'),
            title(['beamX = ',num2str(set_of_beam_X(i_beam_X)),' cm, Tumor depth = ',num2str(set_of_Tumr_Z(i_Tumr_Z)),' cm, Pd = ',num2str(set_of_percnt(i_Pd)),'%'])

            clearvars p_comp y_refc y_tempset_of_mua
        end
        clearvars TheColor
    end
end
clearvars i_Pd i_fig idx_d Y_Limit i_beam_X idx_HCT_Cn i_d i_Tumr_Z

i_fig = 21; idx_d = 1:1:39; Y_Limit = [-4 +4]; i_beam_X = 2; idx_HCT_Cn = 2:1:length(set_of_HCT_Cn);
for i_HCT = 1:length(idx_HCT_Cn)
    TheColor = [(i_HCT-1)/(length(idx_HCT_Cn)-1),0,1-(i_HCT-1)/(length(idx_HCT_Cn)-1)];
    for i_Tumr_Z = 1:length(set_of_Tumr_Z)
        figure(i_fig), subplot(2,2,i_Tumr_Z)
        p_comp =     (set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn(i_HCT),5,idx_d));
        y_refc = (set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn(i_HCT),4,idx_d));
        y_temp = (set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn(i_HCT),2,idx_d));
        y_temp = -log(y_temp./y_refc);
        y_temp(p_comp>p_threshold) = nan;

        plot(set_of_mua(idx_HCT_Cn),zeros(size(set_of_mua(idx_HCT_Cn))),'Color','k','HandleVisibility','off','LineWidth',0.5,'LineStyle','--'), hold on
        plot(d_centers(idx_d),squeeze(nanmean(nanmean(y_temp,1),2)),'Color',TheColor,'DisplayName',['\mu_a = ',num2str(set_of_mua(idx_HCT_Cn(i_HCT))),''],'LineWidth',2,'Marker','o','MarkerFaceColor',TheColor,'MarkerEdgeColor','none','LineStyle','-','MarkerSize',6), hold on
        xlabel('d (cm)'), xlim([min(d_centers) 3])
        ylabel('OD (a.u.)'), ylim([min(Y_Limit) max(Y_Limit)])
        set(gca,'fontsize',12), % legend('show','numColumns',1,'Orientation','horizontal'),
        title(['beamX = ',num2str(set_of_beam_X(i_beam_X)),' cm, Tumor depth = ',num2str(set_of_Tumr_Z(i_Tumr_Z)),' cm'])

        clearvars p_comp y_refc y_tempset_of_mua
    end
    clearvars TheColor
end
clearvars i_fig idx_d Y_Limit i_beam_X idx_HCT_Cn i_d i_Tumr_Z i_HCT









% OD/mu/d vs mua (several d)
idx_d = 1:4:39; Y_Limit = [0 50]; i_beam_X = 2; idx_HCT_Cn = 2:length(set_of_HCT_Cn);
for i_Pd = 1:3
    for i_d = 1:length(idx_d)
        TheColor = [(i_d-1)/(length(idx_d)-1),0,1-(i_d-1)/(length(idx_d)-1)];
        for i_Tumr_Z = 1:length(set_of_Tumr_Z)
            figure(i_Pd+21), subplot(2,2,i_Tumr_Z)
            p_comp = set_of_I_s(i_Pd,1:3,i_beam_X,i_Tumr_Z,idx_HCT_Cn,5,idx_d(i_d));
            y_refc = -log(set_of_I_s(i_Pd,1:3,i_beam_X,i_Tumr_Z,idx_HCT_Cn,4,idx_d(i_d)));
            y_temp = -log(set_of_I_s(i_Pd,1:3,i_beam_X,i_Tumr_Z,idx_HCT_Cn,2,idx_d(i_d)));
            y_temp(p_comp>p_threshold) = nan;

            plot(set_of_mua(idx_HCT_Cn),squeeze(nanmean(nanmean(y_refc,1),2))./set_of_mua(idx_HCT_Cn).'./d_centers(idx_d(i_d)),'Color',TheColor,'HandleVisibility','off','LineWidth',0.5,'LineStyle','--'), hold on
            plot(set_of_mua(idx_HCT_Cn),squeeze(nanmean(nanmean(y_temp,1),2))./set_of_mua(idx_HCT_Cn).'./d_centers(idx_d(i_d)),'Color',TheColor,'DisplayName',['d = ',num2str(d_centers(idx_d(i_d)))],'LineWidth',2,'Marker','o','MarkerFaceColor',TheColor,'MarkerEdgeColor','none','LineStyle','-','MarkerSize',6), hold on
            xlabel('\mu_a (cm^{-1})'), xlim([min(set_of_mua) max(set_of_mua)])
            ylabel('OD/\mu_a/d (cm)'), ylim([min(Y_Limit) max(Y_Limit)]), set(gca,'ytick',0:5:50), % yscale("log")
            set(gca,'fontsize',12), % legend('show','numColumns',1,'Orientation','horizontal'),
            title(['beamX = ',num2str(set_of_beam_X(i_beam_X)),' cm, Tumor depth = ',num2str(set_of_Tumr_Z(i_Tumr_Z)),' cm, Pd = ',num2str(set_of_percnt(i_Pd)),'%'])
            clearvars p_comp y_refc y_temp
        end
        clearvars TheColor
    end
end
clearvars i_Pd i_fig idx_d Y_Limit i_beam_X idx_HCT_Cn i_d i_Tumr_Z

i_fig = 25; idx_d = 1:4:39; Y_Limit = [0 50]; i_beam_X = 2; idx_HCT_Cn = 2:length(set_of_HCT_Cn);
for i_d = 1:length(idx_d)
    TheColor = [(i_d-1)/(length(idx_d)-1),0,1-(i_d-1)/(length(idx_d)-1)];
    for i_Tumr_Z = 1:length(set_of_Tumr_Z)
        figure(i_fig), subplot(2,2,i_Tumr_Z)
        p_comp = set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn,5,idx_d(i_d));
        y_refc = -log(set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn,4,idx_d(i_d)));
        y_temp = -log(set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn,2,idx_d(i_d)));
        y_temp(p_comp>p_threshold) = nan;

        plot(set_of_mua(idx_HCT_Cn),squeeze(nanmean(nanmean(y_refc,1),2))./set_of_mua(idx_HCT_Cn).'./d_centers(idx_d(i_d)),'Color',TheColor,'HandleVisibility','off','LineWidth',0.5,'LineStyle','--'), hold on
        plot(set_of_mua(idx_HCT_Cn),squeeze(nanmean(nanmean(y_temp,1),2))./set_of_mua(idx_HCT_Cn).'./d_centers(idx_d(i_d)),'Color',TheColor,'DisplayName',['d = ',num2str(d_centers(idx_d(i_d)))],'LineWidth',2,'Marker','o','MarkerFaceColor',TheColor,'MarkerEdgeColor','none','LineStyle','-','MarkerSize',6), hold on
        xlabel('\mu_a (cm^{-1})'), xlim([min(set_of_mua) max(set_of_mua)])
        ylabel('OD/\mu_a/d (cm)'), ylim([min(Y_Limit) max(Y_Limit)]), set(gca,'ytick',0:5:50), % yscale("log")
        set(gca,'fontsize',12), % legend('show','numColumns',1,'Orientation','horizontal'),
        title(['beamX = ',num2str(set_of_beam_X(i_beam_X)),' cm, Tumor depth = ',num2str(set_of_Tumr_Z(i_Tumr_Z)),' cm'])

        clearvars p_comp y_refc y_temp
    end
    clearvars TheColor
end
clearvars i_fig idx_d Y_Limit i_beam_X idx_HCT_Cn i_d i_Tumr_Z



% OD/mu/d vs d (several mua)
idx_d = 1:1:39; Y_Limit = [0 50]; i_beam_X = 2; idx_HCT_Cn = 2:1:length(set_of_HCT_Cn);
for i_Pd = 1:3
    for i_HCT = 1:length(idx_HCT_Cn)
        TheColor = [(i_HCT-1)/(length(idx_HCT_Cn)-1),0,1-(i_HCT-1)/(length(idx_HCT_Cn)-1)];
        for i_Tumr_Z = 1:length(set_of_Tumr_Z)
            figure(i_Pd+25), subplot(2,2,i_Tumr_Z)
            p_comp =     (set_of_I_s(i_Pd,1:3,i_beam_X,i_Tumr_Z,idx_HCT_Cn(i_HCT),5,idx_d));
            y_refc = -log(set_of_I_s(i_Pd,1:3,i_beam_X,i_Tumr_Z,idx_HCT_Cn(i_HCT),4,idx_d));
            y_temp = -log(set_of_I_s(i_Pd,1:3,i_beam_X,i_Tumr_Z,idx_HCT_Cn(i_HCT),2,idx_d));
            y_temp(p_comp>p_threshold) = nan;

            plot(d_centers(idx_d),squeeze(nanmean(nanmean(y_refc,1),2))./set_of_mua(idx_HCT_Cn(i_HCT))./d_centers(idx_d).','Color',TheColor,'HandleVisibility','off','LineWidth',0.5,'LineStyle','--'), hold on
            plot(d_centers(idx_d),squeeze(nanmean(nanmean(y_temp,1),2))./set_of_mua(idx_HCT_Cn(i_HCT))./d_centers(idx_d).','Color',TheColor,'DisplayName',['\mu_a = ',num2str(set_of_mua(idx_HCT_Cn(i_HCT))),''],'LineWidth',2,'Marker','o','MarkerFaceColor',TheColor,'MarkerEdgeColor','none','LineStyle','-','MarkerSize',6), hold on
            xlabel('d (cm)'), xlim([min(d_centers) 3])
            ylabel('OD/\mu_a/d (cm)'), ylim([min(Y_Limit) max(Y_Limit)])
            set(gca,'fontsize',12), % legend('show','numColumns',1,'Orientation','horizontal'),
            title(['beamX = ',num2str(set_of_beam_X(i_beam_X)),' cm, Tumor depth = ',num2str(set_of_Tumr_Z(i_Tumr_Z)),' cm, Pd = ',num2str(set_of_percnt(i_Pd)),'%'])

            clearvars p_comp y_refc y_tempset_of_mua
        end
        clearvars TheColor
    end
end
clearvars i_Pd i_fig idx_d Y_Limit i_beam_X idx_HCT_Cn i_d i_Tumr_Z

i_fig = 29; idx_d = 1:1:39; Y_Limit = [0 50]; i_beam_X = 2; idx_HCT_Cn = 2:1:length(set_of_HCT_Cn);
for i_HCT = 1:length(idx_HCT_Cn)
    TheColor = [(i_HCT-1)/(length(idx_HCT_Cn)-1),0,1-(i_HCT-1)/(length(idx_HCT_Cn)-1)];
    for i_Tumr_Z = 1:length(set_of_Tumr_Z)
        figure(i_fig), subplot(2,2,i_Tumr_Z)
        p_comp =     (set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn(i_HCT),5,idx_d));
        y_refc = -log(set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn(i_HCT),4,idx_d));
        y_temp = -log(set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn(i_HCT),2,idx_d));
        y_temp(p_comp>p_threshold) = nan;

        plot(d_centers(idx_d),squeeze(nanmean(nanmean(y_refc,1),2))./set_of_mua(idx_HCT_Cn(i_HCT))./d_centers(idx_d).','Color',TheColor,'HandleVisibility','off','LineWidth',0.5,'LineStyle','--'), hold on
        plot(d_centers(idx_d),squeeze(nanmean(nanmean(y_temp,1),2))./set_of_mua(idx_HCT_Cn(i_HCT))./d_centers(idx_d).','Color',TheColor,'DisplayName',['\mu_a = ',num2str(set_of_mua(idx_HCT_Cn(i_HCT))),''],'LineWidth',2,'Marker','o','MarkerFaceColor',TheColor,'MarkerEdgeColor','none','LineStyle','-','MarkerSize',6), hold on
        xlabel('d (cm)'), xlim([min(d_centers) 3])
        ylabel('OD/\mu_a/d (cm)'), ylim([min(Y_Limit) max(Y_Limit)])
        set(gca,'fontsize',12), % legend('show','numColumns',1,'Orientation','horizontal'),
        title(['beamX = ',num2str(set_of_beam_X(i_beam_X)),' cm, Tumor depth = ',num2str(set_of_Tumr_Z(i_Tumr_Z)),' cm'])

        clearvars p_comp y_refc y_tempset_of_mua
    end
    clearvars TheColor
end
clearvars i_fig idx_d Y_Limit i_beam_X idx_HCT_Cn i_d i_Tumr_Z i_HCT



% OD_relative vs mua (several d)
idx_d = 1:4:39; Y_Limit = [-0.5 +0.5]; i_beam_X = 2; idx_HCT_Cn = 2:length(set_of_HCT_Cn);
for i_Pd = 1:3
    for i_d = 1:length(idx_d)
        TheColor = [(i_d-1)/(length(idx_d)-1),0,1-(i_d-1)/(length(idx_d)-1)];
        for i_Tumr_Z = 1:length(set_of_Tumr_Z)
            figure(i_Pd+29), subplot(2,2,i_Tumr_Z)
            p_comp = set_of_I_s(i_Pd,1:3,i_beam_X,i_Tumr_Z,idx_HCT_Cn,5,idx_d(i_d));
            y_refc = (set_of_I_s(i_Pd,1:3,i_beam_X,i_Tumr_Z,idx_HCT_Cn,4,idx_d(i_d)));
            y_temp = (set_of_I_s(i_Pd,1:3,i_beam_X,i_Tumr_Z,idx_HCT_Cn,2,idx_d(i_d)));
            y_temp = -log(y_temp./y_refc);
            y_temp(p_comp>p_threshold) = nan;

            plot(set_of_mua(idx_HCT_Cn),zeros(size(set_of_mua(idx_HCT_Cn))),  'Color','k','HandleVisibility','off','LineWidth',0.5,'LineStyle','--'), hold on
            plot(set_of_mua(idx_HCT_Cn),squeeze(nanmean(nanmean(y_temp,1),2))./d_centers(idx_d(i_d))./set_of_mua(idx_HCT_Cn).','Color',TheColor,'DisplayName',['d = ',num2str(d_centers(idx_d(i_d)))],'LineWidth',2,'Marker','o','MarkerFaceColor',TheColor,'MarkerEdgeColor','none','LineStyle','-','MarkerSize',6), hold on
            xlabel('\mu_a (cm^{-1})'), xlim([min(set_of_mua) max(set_of_mua)])
            ylabel('OD/\mu_a/d (cm)'), ylim([min(Y_Limit) max(Y_Limit)])
            set(gca,'fontsize',12), % legend('show','numColumns',1,'Orientation','horizontal'),
            title(['beamX = ',num2str(set_of_beam_X(i_beam_X)),' cm, Tumor depth = ',num2str(set_of_Tumr_Z(i_Tumr_Z)),' cm, Pd = ',num2str(set_of_percnt(i_Pd)),'%'])
            clearvars p_comp y_refc y_temp
        end
        clearvars TheColor
    end
end
clearvars i_Pd i_fig idx_d Y_Limit i_beam_X idx_HCT_Cn i_d i_Tumr_Z

i_fig = 33; idx_d = 1:4:39; Y_Limit = [-0.5 +0.5]; i_beam_X = 2; idx_HCT_Cn = 2:length(set_of_HCT_Cn);
for i_d = 1:length(idx_d)
    TheColor = [(i_d-1)/(length(idx_d)-1),0,1-(i_d-1)/(length(idx_d)-1)];
    for i_Tumr_Z = 1:length(set_of_Tumr_Z)
        figure(i_fig), subplot(2,2,i_Tumr_Z)
        p_comp = set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn,5,idx_d(i_d));
        y_refc = (set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn,4,idx_d(i_d)));
        y_temp = (set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn,2,idx_d(i_d)));
        y_temp = -log(y_temp./y_refc);
        y_temp(p_comp>p_threshold) = nan;

        plot(set_of_mua(idx_HCT_Cn),zeros(size(set_of_mua(idx_HCT_Cn))),  'Color','k','HandleVisibility','off','LineWidth',0.5,'LineStyle','--'), hold on
        plot(set_of_mua(idx_HCT_Cn),squeeze(nanmean(nanmean(y_temp,1),2))./d_centers(idx_d(i_d))./set_of_mua(idx_HCT_Cn).','Color',TheColor,'DisplayName',['d = ',num2str(d_centers(idx_d(i_d)))],'LineWidth',2,'Marker','o','MarkerFaceColor',TheColor,'MarkerEdgeColor','none','LineStyle','-','MarkerSize',6), hold on
        xlabel('\mu_a (cm^{-1})'), xlim([min(set_of_mua) max(set_of_mua)])
        ylabel('OD/\mu_a/d (cm)'), ylim([min(Y_Limit) max(Y_Limit)])
        set(gca,'fontsize',12), % legend('show','numColumns',1,'Orientation','horizontal'),
        title(['beamX = ',num2str(set_of_beam_X(i_beam_X)),' cm, Tumor depth = ',num2str(set_of_Tumr_Z(i_Tumr_Z)),' cm'])

        clearvars p_comp y_refc y_temp
    end
    clearvars TheColor
end
clearvars i_fig idx_d Y_Limit i_beam_X idx_HCT_Cn i_d i_Tumr_Z



% OD_relative vs d (several mua)
idx_d = 1:1:39; Y_Limit = [-0.5 +0.5]; i_beam_X = 2; idx_HCT_Cn = 2:1:length(set_of_HCT_Cn);
for i_Pd = 1:3
    for i_HCT = 1:length(idx_HCT_Cn)
        TheColor = [(i_HCT-1)/(length(idx_HCT_Cn)-1),0,1-(i_HCT-1)/(length(idx_HCT_Cn)-1)];
        for i_Tumr_Z = 1:length(set_of_Tumr_Z)
            figure(i_Pd+33), subplot(2,2,i_Tumr_Z)
            p_comp =     (set_of_I_s(i_Pd,1:3,i_beam_X,i_Tumr_Z,idx_HCT_Cn(i_HCT),5,idx_d));
            y_refc = (set_of_I_s(i_Pd,1:3,i_beam_X,i_Tumr_Z,idx_HCT_Cn(i_HCT),4,idx_d));
            y_temp = (set_of_I_s(i_Pd,1:3,i_beam_X,i_Tumr_Z,idx_HCT_Cn(i_HCT),2,idx_d));
            y_temp = -log(y_temp./y_refc);
            y_temp(p_comp>p_threshold) = nan;

            plot(set_of_mua(idx_HCT_Cn),zeros(size(set_of_mua(idx_HCT_Cn))),'Color','k','HandleVisibility','off','LineWidth',0.5,'LineStyle','--'), hold on
            plot(d_centers(idx_d),squeeze(nanmean(nanmean(y_temp,1),2))/set_of_mua(idx_HCT_Cn(i_HCT))./d_centers(idx_d).',    'Color',TheColor,'DisplayName',['\mu_a = ',num2str(set_of_mua(idx_HCT_Cn(i_HCT))),''],'LineWidth',2,'Marker','o','MarkerFaceColor',TheColor,'MarkerEdgeColor','none','LineStyle','-','MarkerSize',6), hold on
            xlabel('d (cm)'), xlim([min(d_centers) 3])
            ylabel('OD/\mu_a/d (cm)'), ylim([min(Y_Limit) max(Y_Limit)])
            set(gca,'fontsize',12), % legend('show','numColumns',1,'Orientation','horizontal'),
            title(['beamX = ',num2str(set_of_beam_X(i_beam_X)),' cm, Tumor depth = ',num2str(set_of_Tumr_Z(i_Tumr_Z)),' cm, Pd = ',num2str(set_of_percnt(i_Pd)),'%'])

            clearvars p_comp y_refc y_tempset_of_mua
        end
        clearvars TheColor
    end
end
clearvars i_Pd i_fig idx_d Y_Limit i_beam_X idx_HCT_Cn i_d i_Tumr_Z

i_fig = 37; idx_d = 1:1:39; Y_Limit = [-0.5 +0.5]; i_beam_X = 2; idx_HCT_Cn = 2:1:length(set_of_HCT_Cn);
for i_HCT = 1:length(idx_HCT_Cn)
    TheColor = [(i_HCT-1)/(length(idx_HCT_Cn)-1),0,1-(i_HCT-1)/(length(idx_HCT_Cn)-1)];
    for i_Tumr_Z = 1:length(set_of_Tumr_Z)
        figure(i_fig), subplot(2,2,i_Tumr_Z)
        p_comp =     (set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn(i_HCT),5,idx_d));
        y_refc = (set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn(i_HCT),4,idx_d));
        y_temp = (set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn(i_HCT),2,idx_d));
        y_temp = -log(y_temp./y_refc);
        y_temp(p_comp>p_threshold) = nan;

        plot(set_of_mua(idx_HCT_Cn),zeros(size(set_of_mua(idx_HCT_Cn))),'Color','k','HandleVisibility','off','LineWidth',0.5,'LineStyle','--'), hold on
        plot(d_centers(idx_d),squeeze(nanmean(nanmean(y_temp,1),2))/set_of_mua(idx_HCT_Cn(i_HCT))./d_centers(idx_d).','Color',TheColor,'DisplayName',['\mu_a = ',num2str(set_of_mua(idx_HCT_Cn(i_HCT))),''],'LineWidth',2,'Marker','o','MarkerFaceColor',TheColor,'MarkerEdgeColor','none','LineStyle','-','MarkerSize',6), hold on
        xlabel('d (cm)'), xlim([min(d_centers) 3])
        ylabel('OD/\mu_a/d (cm)'), ylim([min(Y_Limit) max(Y_Limit)])
        set(gca,'fontsize',12), % legend('show','numColumns',1,'Orientation','horizontal'),
        title(['beamX = ',num2str(set_of_beam_X(i_beam_X)),' cm, Tumor depth = ',num2str(set_of_Tumr_Z(i_Tumr_Z)),' cm'])

        clearvars p_comp y_refc y_tempset_of_mua
    end
    clearvars TheColor
end
clearvars i_fig idx_d Y_Limit i_beam_X idx_HCT_Cn i_d i_Tumr_Z i_HCT













figure(100)
for i_percent = 1:length(set_of_percnt)
    for i_versns = 1:length(set_of_versns)
        for i_beam_X = 2
            for i_Tumr_Z = 1:length(set_of_Tumr_Z)
                for i_HCT_Cn = 1:length(set_of_HCT_Cn)
                    x_temp = squeeze(set_of_I_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,1,:));
                    y_temp = squeeze(set_of_I_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,2,:));
                    x_refc = squeeze(set_of_I_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,3,:));
                    y_refc = squeeze(set_of_I_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,4,:));
                    p_valu = squeeze(set_of_I_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,5,:));
                    y_temp(p_valu>p_threshold) = nan;
                    y_refc(p_valu>p_threshold) = nan;

                    subplot(2,2,1), plot(x_temp,-log(y_temp),                              'Color',[i_HCT_Cn/length(set_of_HCT_Cn) 0 i_Tumr_Z/length(set_of_Tumr_Z)]), hold on, axis([0 4 0 020]), xlabel('d (cm)'), ylabel('OD (a.u.)'), title('inclusion')
                    subplot(2,2,2), plot(x_refc,-log(y_refc),                              'Color',[i_HCT_Cn/length(set_of_HCT_Cn) 0 i_Tumr_Z/length(set_of_Tumr_Z)]), hold on, axis([0 4 0 020]), xlabel('d (cm)'), ylabel('OD (a.u.)'), title('no tumor')
                    subplot(2,2,3), plot(x_temp,-log(y_temp)./set_of_mua(i_HCT_Cn)./x_temp,'Color',[i_HCT_Cn/length(set_of_HCT_Cn) 0 i_Tumr_Z/length(set_of_Tumr_Z)]), hold on, axis([0 4 0 140]), xlabel('d (cm)'), ylabel('OD/\mu_a/d (a.u.)'), title('inclusion')
                    subplot(2,2,4), plot(x_refc,-log(y_refc)./set_of_mua(i_HCT_Cn)./x_refc,'Color',[i_HCT_Cn/length(set_of_HCT_Cn) 0 i_Tumr_Z/length(set_of_Tumr_Z)]), hold on, axis([0 4 0 140]), xlabel('d (cm)'), ylabel('OD/\mu_a/d (a.u.)'), title('no tumor')
                    clearvars x_temp y_temp x_refc y_refc
                end
            end
        end
    end
end
clearvars i_percent i_versns i_beam_X i_Tumr_Z i_HCT_Cn

figure(101)
for i_percent = 1:length(set_of_percnt)
    for i_versns = 1:length(set_of_versns)
        for i_beam_X = 2
            for i_Tumr_Z = 1:length(set_of_Tumr_Z)
                for i_HCT_Cn = 1:length(set_of_HCT_Cn)
                    x_temp = squeeze(set_of_I_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,1,:));
                    y_temp = squeeze(set_of_I_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,2,:));
                    x_refc = squeeze(set_of_I_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,3,:));
                    y_refc = squeeze(set_of_I_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,4,:));
                    p_valu = squeeze(set_of_I_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,5,:));
                    y_temp(p_valu>p_threshold) = nan;
                    y_refc(p_valu>p_threshold) = nan;

                    subplot(2,1,1), plot(x_temp,-log(y_temp./y_refc),                              'Color',[i_HCT_Cn/length(set_of_HCT_Cn) 0 i_Tumr_Z/length(set_of_Tumr_Z)]), hold on, axis([0 4 -06 +06]), xlabel('d (cm)'), ylabel('OD (a.u.)'), title('inclusion wrt no-tumor')
                    subplot(2,1,2), plot(x_temp,-log(y_temp./y_refc)./set_of_mua(i_HCT_Cn)./x_refc,'Color',[i_HCT_Cn/length(set_of_HCT_Cn) 0 i_Tumr_Z/length(set_of_Tumr_Z)]), hold on, axis([0 4 -10 +60]), xlabel('d (cm)'), ylabel('OD/\mu_a/d (a.u.)'), title('inclusion wrt no-tumor')
                    clearvars x_temp y_temp x_refc y_refc
                end
            end
        end
    end
end
clearvars i_percent i_versns i_beam_X i_Tumr_Z i_HCT_Cn

figure(102)
set_of_DPF_s = nan(size(set_of_I_s,1),size(set_of_I_s,2),size(set_of_I_s,3),size(set_of_I_s,4),size(set_of_I_s,5),size(set_of_I_s,7));
for i_percent = 1:length(set_of_percnt)
    for i_versns = 1:length(set_of_versns)
        for i_beam_X = 1:length(set_of_beam_X)
            for i_Tumr_Z = 1:length(set_of_Tumr_Z)
                for i_HCT_Cn = 1:length(set_of_HCT_Cn)
                    x_temp = squeeze(set_of_I_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,1,:));
                    y_temp = squeeze(set_of_I_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,2,:));
                    x_refc = squeeze(set_of_I_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,3,:));
                    y_refc = squeeze(set_of_I_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,4,:));
                    p_valu = squeeze(set_of_I_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,5,:));
                    y_temp(p_valu>p_threshold) = nan;
                    y_refc(p_valu>p_threshold) = nan;

                                        set_of_DPF_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,:) = -log(y_temp./y_refc)./x_temp./set_of_mua(i_HCT_Cn);
                    plot(x_temp,squeeze(set_of_DPF_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,:)),'Color',[i_HCT_Cn/length(set_of_HCT_Cn) 0 i_Tumr_Z/length(set_of_Tumr_Z)]), hold on, axis([0 4 -10 +60]), xlabel('d (cm)'), ylabel('DPF (a.u.)'), title('inclusion wrt no-tumor')
                    clearvars x_temp y_temp x_refc y_refc
                end
            end
        end
    end
end

figure(103)
set_of_err_s = nan(size(set_of_I_s,1),size(set_of_I_s,2),size(set_of_I_s,3),size(set_of_I_s,4),size(set_of_I_s,5),size(set_of_I_s,7));
TheDPF = squeeze(nanmean(nanmean(nanmean(nanmean(set_of_DPF_s(:,:,2,:,:,:),1),2),4),5));
for i_percent = 1:length(set_of_percnt)
    for i_versns = 1:length(set_of_versns)
        for i_beam_X = 2
            for i_Tumr_Z = 1:length(set_of_Tumr_Z)
                for i_HCT_Cn = 1:length(set_of_HCT_Cn)
                    x_temp = squeeze(set_of_I_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,1,:));
                    y_temp = squeeze(set_of_I_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,2,:));
                    x_refc = squeeze(set_of_I_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,3,:));
                    y_refc = squeeze(set_of_I_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,4,:));
                    p_valu = squeeze(set_of_I_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,5,:));
                    y_temp(p_valu>p_threshold) = nan;
                    y_refc(p_valu>p_threshold) = nan;

                    set_of_err_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,:) = (-log(y_temp./y_refc)./x_temp./squeeze(set_of_DPF_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,:))-set_of_mua(i_HCT_Cn))/set_of_mua(i_HCT_Cn)*100;
                    plot(x_temp,squeeze(set_of_err_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,:)),'Color',[i_HCT_Cn/length(set_of_HCT_Cn) 0 1-i_HCT_Cn/length(set_of_HCT_Cn)]), hold on, xlim([0 4]), xlabel('d (cm)'), ylabel('error (%)'), title('inclusion wrt no-tumor')
                    clearvars x_temp y_temp x_refc y_refc
                end
            end
        end
    end
end

figure(104)
set_of_err_s = nan(size(set_of_I_s,1),size(set_of_I_s,2),size(set_of_I_s,3),size(set_of_I_s,4),size(set_of_I_s,5),size(set_of_I_s,7));
TheDPF = squeeze(nanmean(nanmean(nanmean(nanmean(set_of_DPF_s(:,:,2,:,:,:),1),2),4),5));
for i_percent = 1:length(set_of_percnt)
    for i_versns = 1:length(set_of_versns)
        for i_beam_X = 2
            for i_Tumr_Z = 1:length(set_of_Tumr_Z)
                for i_HCT_Cn = 1:length(set_of_HCT_Cn)
                    x_temp = squeeze(set_of_I_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,1,:));
                    y_temp = squeeze(set_of_I_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,2,:));
                    x_refc = squeeze(set_of_I_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,3,:));
                    y_refc = squeeze(set_of_I_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,4,:));
                    p_valu = squeeze(set_of_I_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,5,:));
                    y_temp(p_valu>p_threshold) = nan;
                    y_refc(p_valu>p_threshold) = nan;

                    set_of_err_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,:) = (-log(y_temp./y_refc)./x_temp./squeeze(TheDPF)-set_of_mua(i_HCT_Cn))/set_of_mua(i_HCT_Cn)*100;
                    plot(x_temp,squeeze(set_of_err_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,:)),'Color',[i_HCT_Cn/length(set_of_HCT_Cn) 0 1-i_HCT_Cn/length(set_of_HCT_Cn)]), hold on, xlim([0 4]), xlabel('d (cm)'), ylabel('error (%)'), title('inclusion wrt no-tumor')
                    clearvars x_temp y_temp x_refc y_refc
                end
            end
        end
    end
end











for i_fig = 17:-1:1, figure(i_fig), saveas(gcf,['Photon_',num2str(code_num),'_5_fig_',sprintf('%1.0f',i_fig),'.jpg' ],'jpeg'); end
end

function [mua,mus,g,musp] = calc_muas_based_HCT(HCT,lambda)
% Martina Meinke, Gerhard Müller, Jürgen Helfmann, and Moritz Friebel, "Empirical model functions to calculate hematocrit-dependent optical properties of human blood," Appl. Opt. 46, 1742-1753 (2007)
[mua_st,mus_st,g_st,musp_st] = Data_of_Standard_Optical_Parameters (lambda);
mua = nan(size(lambda));
mus = nan(size(lambda));
g   = nan(size(lambda));
musp= nan(size(lambda));
for i_lambda = 1:length(lambda)
    % mu_a
    if     ((250<=lambda(i_lambda)&&lambda(i_lambda)<=400)||(430<=lambda(i_lambda)&&lambda(i_lambda)<=600))&&(0.84<=HCT&&HCT<=42.1)
        mua(i_lambda) = 0.1233.*mua_st(i_lambda).*HCT;
    elseif  (600<=lambda(i_lambda)&&lambda(i_lambda)<=1100)&&(0.84<=HCT&&HCT<=42.1)
        mua(i_lambda) = 0.1206.*mua_st(i_lambda).*HCT;
    elseif  (400< lambda(i_lambda)&&lambda(i_lambda)< 430 )&&(0.84<=HCT&&HCT<=42.1)
        mua(i_lambda) = mean([0.1233 0.1206]).*mua_st(i_lambda).*HCT;
    else
        mua(i_lambda) = nan;
    end
    % mu_s
    if      (250<=lambda(i_lambda)&&lambda(i_lambda)<=1100)&&(0.84<=HCT&&HCT<=17.1)
        mus(i_lambda) = (-0.0015.*HCT.^2+0.1268.*HCT).*mus_st(i_lambda);
    elseif  (600<=lambda(i_lambda)&&lambda(i_lambda)<=1100)&&(17.1<=HCT&&HCT<=42.1)
        mus(i_lambda) = (-0.0015.*HCT.^2+0.1268.*HCT).*mus_st(i_lambda);
    elseif  (250<=lambda(i_lambda)&&lambda(i_lambda)< 600 )&&(17.1<=HCT&&HCT<=42.1)
        mus(i_lambda) = (-0.0015.*HCT.^2+0.1268.*HCT).*mus_st(i_lambda);
    else
        mus(i_lambda) = nan;
    end
    % mu_sp
    if      (250<=lambda(i_lambda)&&lambda(i_lambda)<=1100)&&(0.84<=HCT&&HCT<=17.1)
        musp(i_lambda) = 0.1167.*HCT.*musp_st(i_lambda);
    elseif  (600<=lambda(i_lambda)&&lambda(i_lambda)<=1100)&&(17.1<=HCT&&HCT<=42.1)
        musp(i_lambda) = 0.1167.*HCT.*musp_st(i_lambda);
    elseif  (250<=lambda(i_lambda)&&lambda(i_lambda)< 600 )&&(17.1<=HCT&&HCT<=42.1)
        musp(i_lambda) = 0.1167.*HCT.*musp_st(i_lambda);
    else
        musp(i_lambda) = nan;
    end
    % g
    g = (1-musp./mus);
    % if      (600<=lambda(i_lambda)&&lambda(i_lambda)<=1100)&&(0.84<=HCT&&HCT<=42.1)
    %     g(i_lambda) = (((-2.684e-6).*HCT.^2)+((-2.373e-4).*HCT)+1.003).*g_st(i_lambda);
    % elseif  (250<=lambda(i_lambda)&&lambda(i_lambda)< 600 )&&(0.84<=HCT&&HCT<=42.1)
    %     g(i_lambda) = (((-2.684e-6).*HCT.^2)+((-2.373e-4).*HCT)+1.003).*g_st(i_lambda);
    % else
    %     g(i_lambda) = nan;
    % end
end
end
function [mua_st,mus_st,g_st,musp_st] = Data_of_Standard_Optical_Parameters (lmbda)
% Martina Meinke, Gerhard Müller, Jürgen Helfmann, and Moritz Friebel, "Empirical model functions to calculate hematocrit-dependent optical properties of human blood," Appl. Opt. 46, 1742-1753 (2007)
db = [
250	7.52	27.9	0.877	3.43
255	7.93	27.2	0.88	3.27
260	8.6	26.6	0.883	3.1
265	9.41	25.7	0.883	3.01
270	9.82	25.4	0.882	2.99
275	9.93	25.2	0.884	2.93
280	9.74	25.3	0.889	2.8
285	9.19	25.5	0.899	2.59
290	8.24	25.9	0.912	2.28
295	7.04	26.7	0.9268	1.96
300	6	27.4	0.9389	1.68
305	5.52	27.7	0.9449	1.53
310	5.56	27.5	0.9442	1.54
315	5.95	26.9	0.9398	1.62
320	6.5	26	0.9335	1.73
325	7.13	25.4	0.928	1.83
330	7.69	25.2	0.9239	1.92
335	8.1	24.9	0.9208	1.97
340	8.38	24.7	0.9192	1.99
345	8.48	24.6	0.9191	1.99
350	8.37	24.6	0.9215	1.93
355	8.07	24.7	0.9257	1.84
360	7.69	24.7	0.9308	1.71
365	7.4	24.7	0.9345	1.61
370	7.37	24.5	0.9356	1.58
375	7.72	24	0.9323	1.62
380	8.58	23.2	0.9242	1.76
385	9.98	22.2	0.9105	1.99
390	11.86	21.1	0.892	2.29
395	14.21	20.1	0.869	2.64
400	16.86	19.1	0.845	2.96
405	19.66	18.3	0.824	3.23
410	21.8	17.9	0.81	3.39
415	22.66	18.6	0.812	3.48
420	21.54	19.8	0.831	3.35
425	18.62	21.8	0.861	3.03
430	14.93	23.8	0.893	2.55
435	11.52	25.3	0.92	2.02
440	8.95	26.8	0.9405	1.6
445	7.11	28	0.9543	1.28
450	5.79	29	0.9629	1.08
455	4.83	30	0.9685	0.946
460	4.11	30.5	0.9724	0.841
465	3.56	31	0.9756	0.757
470	3.13	31.3	0.9777	0.699
475	2.79	31.5	0.9793	0.652
480	2.52	31.8	0.9806	0.619
490	2.2	32	0.982	0.575
500	2.04	32.2	0.9831	0.545
510	1.98	32.2	0.9835	0.531
520	2.51	31.5	0.9838	0.51
530	3.97	29.5	0.9794	0.609
540	5.05	28.2	0.9755	0.691
550	4.4	29	0.9779	0.642
560	3.6	30.1	0.9804	0.59
570	4.56	28.8	0.9777	0.641
580	4.56	28.9	0.9771	0.662
590	1.9	32.1	0.9827	0.556
600	0.478	33.9	0.9854	0.496
610	0.17	34.2	0.9858	0.487
620	0.0812	34.3	0.9861	0.477
630	0.0496	34.2	0.9863	0.469
640	0.0348	34.3	0.9865	0.462
660	0.0251	34.1	0.9868	0.45
670	0.0239	33.9	0.9871	0.439
690	0.0243	33.6	0.9872	0.43
700	0.0246	33.4	0.9872	0.427
720	0.0284	32.9	0.9871	0.426
740	0.036	32.4	0.987	0.422
760	0.0461	31.8	0.9868	0.42
780	0.0558	31.2	0.9867	0.415
800	0.0641	30.8	0.9868	0.407
820	0.0762	30.5	0.9867	0.407
840	0.0853	30.5	0.9867	0.405
860	0.0953	29.7	0.9864	0.403
880	0.104	29.8	0.9861	0.413
900	0.106	28.4	0.9857	0.405
920	0.111	28.1	0.9854	0.411
940	0.117	27.7	0.9847	0.423
960	0.125	26.8	0.984	0.429
980	0.133	26.1	0.9836	0.428
1000	0.128	25.8	0.9837	0.421
1020	0.118	25.7	0.9839	0.412
1040	0.104	25.3	0.9837	0.412
1060	0.0879	25	0.9838	0.406
1080	0.0795	24.6	0.9841	0.392
1100	0.074	24.5	0.9842	0.387
];
mua_st  = interp1(db(:,1),db(:,2),lmbda,'linear','extrap'); mua_st  =  mua_st*10;
mus_st  = interp1(db(:,1),db(:,3),lmbda,'linear','extrap'); mus_st  =  mus_st*10;
g_st    = interp1(db(:,1),db(:,4),lmbda,'linear','extrap'); g_st    =       g_st;
musp_st = interp1(db(:,1),db(:,5),lmbda,'linear','extrap'); musp_st = musp_st*10;
end
function [n_water] = refractive_index_water(lambda_nm)
% Computes the refractive index of water as a function of wavelength in nm
% George M. Hale and Marvin R. Querry, "Optical Constants of Water in the 200-nm to 200-μm Wavelength Region," Appl. Opt. 12, 555-563 (1973)
% from https://refractiveindex.info/?book=H2O&page=Hale&shelf=main&utm_source=chatgpt.com
db = [...
    0.200	1.396
    0.225	1.373
    0.250	1.362
    0.275	1.354
    0.300	1.349
    0.325	1.346
    0.350	1.343
    0.375	1.341
    0.400	1.339
    0.425	1.338
    0.450	1.337
    0.475	1.336
    0.500	1.335
    0.525	1.334
    0.550	1.333
    0.575	1.333
    0.600	1.332
    0.625	1.332
    0.650	1.331
    0.675	1.331
    0.700	1.331
    0.725	1.330
    0.750	1.330
    0.775	1.330
    0.800	1.329
    0.825	1.329
    0.850	1.329
    0.875	1.328
    0.900	1.328
    0.925	1.328
    0.950	1.327
    0.975	1.327
    1.0	1.327
    1.2	1.324
    1.4	1.321
    1.6	1.317
    1.8	1.312
    2.0	1.306
    2.2	1.296
    2.4	1.279
    2.6	1.242
    2.65	1.219
    2.70	1.188
    2.75	1.157
    2.80	1.142
    2.85	1.149
    2.90	1.201
    2.95	1.292
    3.00	1.371
    3.05	1.426
    3.10	1.467
    3.15	1.483
    3.20	1.478
    3.25	1.467
    3.30	1.450
    3.35	1.432
    3.40	1.420
    3.45	1.410
    3.50	1.400
    3.6	1.385
    3.7	1.374
    3.8	1.364
    3.9	1.357
    4.0	1.351
    4.1	1.346
    4.2	1.342
    4.3	1.338
    4.4	1.334
    4.5	1.332
    4.6	1.330
    4.7	1.330
    4.8	1.330
    4.9	1.328
    5.0	1.325
    5.1	1.322
    5.2	1.317
    5.3	1.312
    5.4	1.305
    5.5	1.298
    5.6	1.289
    5.7	1.277
    5.8	1.262
    5.9	1.248
    6.0	1.265
    6.1	1.319
    6.2	1.363
    6.3	1.357
    6.4	1.347
    6.5	1.339
    6.6	1.334
    6.7	1.329
    6.8	1.324
    6.9	1.321
    7.0	1.317
    7.1	1.314
    7.2	1.312
    7.3	1.309
    7.4	1.307
    7.5	1.304
    7.6	1.302
    7.7	1.299
    7.8	1.297
    7.9	1.294
    8.0	1.291
    8.2	1.286
    8.4	1.281
    8.6	1.275
    8.8	1.269
    9.0	1.262
    9.2	1.255
    9.4	1.247
    9.6	1.239
    9.8	1.229
    10.0	1.218
    10.5	1.185
    11.0	1.153
    11.5	1.126
    12.0	1.111
    12.5	1.123
    13.0	1.146
    13.5	1.177
    14.0	1.210
    14.5	1.241
    15.0	1.270
    15.5	1.297
    16.0	1.325
    16.5	1.351
    17.0	1.376
    17.5	1.401
    18.0	1.423
    18.5	1.443
    19.0	1.461
    19.5	1.476
    20.0	1.480
    21.0	1.487
    22	1.500
    23	1.511
    24	1.521
    25	1.531
    26	1.539
    27	1.545
    28	1.549
    29	1.551
    30	1.551
    32	1.546
    34	1.536
    36	1.527
    38	1.522
    40	1.519
    42	1.522
    44	1.530
    46	1.541
    48	1.555
    50	1.587
    60	1.703
    70	1.821
    80	1.886
    90	1.924
    100	1.957
    110	1.966
    120	2.004
    130	2.036
    140	2.056
    150	2.069
    160	2.081
    170	2.094
    180	2.107
    190	2.119
    200	2.130
    ];
n_water = interp1(db(:,1).*1000,db(:,2),lambda_nm,'linear','extrap');
end
function [beta_st] = Specific_Refractive_Increment_beta(lmbda)
% Moritz Friebel and Martina Meinke, "Model function to calculate the refractive index of native hemoglobin in the wavelength range of 250-1100 nm dependent on concentration," Appl. Opt. 45, 2838-2842 (2006)
db = [...
250	0.00221
255	0.002155
260	0.002105
265	0.002069
270	0.002048
275	0.002042
280	0.002044
285	0.002047
290	0.002047
295	0.002037
300	0.00202
305	0.001999
310	0.001998
320	0.002007
330	0.002021
340	0.00201
350	0.001989
355	0.001985
360	0.001983
365	0.001912
370	0.00186
375	0.001816
380	0.001774
385	0.001732
390	0.001694
395	0.001668
400	0.001664
405	0.001701
410	0.001799
415	0.001985
420	0.002117
425	0.002195
430	0.002273
435	0.002227
440	0.00221
445	0.002184
450	0.002156
455	0.002131
460	0.002109
465	0.002092
470	0.002078
475	0.002067
480	0.002056
485	0.002045
490	0.002033
495	0.002019
500	0.002005
510	0.002009
520	0.001983
530	0.001966
540	0.001981
550	0.001998
560	0.001992
570	0.001988
580	0.002004
590	0.002015
600	0.001988
610	0.001967
620	0.001964
630	0.00196
640	0.001954
660	0.001958
680	0.00197
700	0.001992
720	0.001979
740	0.001955
760	0.001958
780	0.00196
800	0.001939
820	0.00192
840	0.001935
860	0.001951
880	0.001982
900	0.001998
920	0.002011
940	0.002015
960	0.002021
980	0.002017
1000	0.002052
1020	0.002049
1040	0.002044
1060	0.00204
1080	0.002044
1100	0.002056
];
beta_st = interp1(db(:,1),db(:,2),lmbda,'linear','extrap');
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
function params = fit2DFunction(fun, x, y, z, init_guess, lb, ub)
% fit2DFunction Fits a parametric function fun(x, y, params) to data
%
% Inputs:
%   fun        - function handle, e.g., @(p,x,y) yourModel(x,y,p)
%   x, y       - input data (same size as z or vectorized)
%   z          - observed output data (same size as x and y)
%   init_guess - initial guess for parameters (vector)
%   lb, ub     - (optional) lower and upper bounds for parameters
%
% Output:
%   params     - fitted parameters

    % Flatten x, y, z if not already
    x = x(:); y = y(:); z = z(:);

    % Objective function: maps params to model output at (x,y)
    model = @(p, xy) fun(p, xy(:,1), xy(:,2));

    % Combine x and y into one matrix
    xy = [x, y];

    if nargin < 6
        % No bounds
        params = lsqcurvefit(model, init_guess, xy, z);
    else
        % With bounds
        params = lsqcurvefit(model, init_guess, xy, z, lb, ub);
    end
end
