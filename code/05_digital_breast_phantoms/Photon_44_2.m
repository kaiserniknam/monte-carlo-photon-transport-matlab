function [] = Photon_44_2 ()
% Repository group: 05_digital_breast_phantoms
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: simulating phantoms prepared by Diego
% Similar to Photon_42: comparisons between Tumor and No Tumor cases, with the tumor composed of blood. The source is placed at four different locations. Tumor radius is 1.5/2 cm. Simulations are run for different hematocrit (HCT) levels.
% analyzing simulated data: only for source on top of tumor at Z = 1.25 cm

clc
close all

% Optical properties of breast and tumor [1-4]
% mua_glnd = 0.55; mus_glnd = 332.7; g_glnd = 0.965; n_glnd = 1.4;
% mua_adps = 0.82; mus_adps = 313.8; g_adps = 0.976; n_adps = 1.4;
set_of_percnt = [25,50,75]; % percent of adipose
set_of_versns = [0,1,2]; % version of each density
set_of_beam_X = [-6,-3,0,+3]; % beam_X location (in cm)
set_of_Tumr_Z = [0.75,1.25,1.75]; % Tumor_Z location (in cm)
set_of_HCT_Cn = [0,10,20,30,40]; % HCT (in %)
N_bins = 151; Lx = 20.4; Ly = 6.48; Lz = 5.16;
d_edges = linspace(0,sqrt((Lx/2).^2+(Ly/2).^2+(Lz/2).^2),N_bins+1);
[C,D] = ndgrid(set_of_HCT_Cn,1/2.*(d_edges(1:end-1)+d_edges(2:end-0)));
set_of_Is = nan(length(set_of_percnt),length(set_of_versns),length(set_of_HCT_Cn),length(d_edges)-1);

for i_percent = 1:length(set_of_percnt)
    for i_versns = 1:length(set_of_versns)
        for i_beam_X = 3 % 1:length(set_of_beam_X)
            for i_Tumr_Z = 2 % 1:length(set_of_Tumr_Z)
                for i_HCT_Cn = 1:length(set_of_HCT_Cn)
                    the_percnt = set_of_percnt(i_percent);
                    the_verson = set_of_versns(i_versns);
                    the_beam_X = set_of_beam_X(i_beam_X);
                    the_Tumr_Z = set_of_Tumr_Z(i_Tumr_Z);
                    the_HCT_Cn = set_of_HCT_Cn(i_HCT_Cn);
                    data = load(['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_44_Pd',num2str(the_percnt),'_',num2str(the_verson),'_',num2str(the_HCT_Cn),'_',num2str(the_beam_X),'_',num2str(the_Tumr_Z),'_y_-z_x.mat']);


                    i_fig = 1; fun_x = @mean; fun_y = @sum; TheOutFun = @(x)(-log(x)); x_label = 'd (cm)'; y_label = 'OD (a.u.)';
                    the_Color = [(i_HCT_Cn-1)/(length(set_of_HCT_Cn)-1),0,1-(i_HCT_Cn-1)/(length(set_of_HCT_Cn)-1)];
                    % Lx = size(data.M_raw,1)*data.dl; Ly = size(data.M_raw,2)*data.dl; Lz = size(data.M_raw,3)*data.dl; disp(num2str([Lx Ly Lz]))
                    x_temp = sqrt(...
                        (data.p_in(:,1)-data.p_ot(:,1)).^2+ ...
                        (data.p_in(:,2)-data.p_ot(:,2)).^2+ ...
                        (data.p_in(:,3)-data.p_ot(:,3)).^2); y_temp = data.w_ot/data.no_of_photons;
                    [~,~,index_in] = histcounts(x_temp,d_edges); clearvars Lx Ly Lz
                    x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
                    y_bind = accumarray(index_in,y_temp,[],fun_y,nan);
                    % I vs. d
                    figure(i_fig), subplot(3,3,(i_percent-1)*3+i_versns)
                    set_of_Is(i_percent,i_versns,i_HCT_Cn,1:length(x_bind)) = TheOutFun(y_bind);
                    plot(x_bind,TheOutFun(y_bind),'Color',the_Color,'DisplayName',['HCT = ',num2str(the_HCT_Cn)],'Marker','none','LineStyle','-','MarkerEdgeColor','none','MarkerFaceColor','none','MarkerSize',8,'LineWidth',1), hold on
                    if i_percent==1, title(['version = ',num2str(the_verson)]); else, title("");end
                    if i_percent==length(set_of_percnt), set(gca,'xtick',0:5); xlabel(x_label); else, set(gca,'xtick',0:5); xlabel(''); end
                    if i_versns==1, ylabel([{['Pd = ',num2str(the_percnt)]},{y_label}]); else, ylabel(""); end
                    set(gca,'fontsize',12), axis tight, grid on; axis([0 5 0 20])
                    legend('show','Location','southeast','NumColumns',2)
                    clearvars i_fig index_in x_temp y_temp x_bind y_bind fun_x fun_y TheOutFun x_label y_label


                    i_fig = 2; fun_x = @mean; fun_y = @sum; TheOutFun = @(x)(-log(x)); x_label = 'd (cm)'; y_label = 'OD/d (cm^{-1})';
                    % Lx = size(data.M_raw,1)*data.dl; Ly = size(data.M_raw,2)*data.dl; Lz = size(data.M_raw,3)*data.dl;
                    x_temp = sqrt(...
                        (data.p_in(:,1)-data.p_ot(:,1)).^2+ ...
                        (data.p_in(:,2)-data.p_ot(:,2)).^2+ ...
                        (data.p_in(:,3)-data.p_ot(:,3)).^2); y_temp = data.w_ot/data.no_of_photons;
                    [~,~,index_in] = histcounts(x_temp,d_edges); clearvars Lx Ly Lz
                    x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
                    y_bind = accumarray(index_in,y_temp,[],fun_y,nan);
                    % I vs. d
                    figure(i_fig), subplot(3,3,(i_percent-1)*3+i_versns)
                    plot(x_bind,TheOutFun(y_bind)./x_bind,'Color',the_Color,'DisplayName',['HCT = ',num2str(the_HCT_Cn)],'Marker','none','LineStyle','-','MarkerEdgeColor','none','MarkerFaceColor','none','MarkerSize',8,'LineWidth',1), hold on
                    if i_percent==1, title(['version = ',num2str(the_verson)]); else, title("");end
                    if i_percent==length(set_of_percnt), set(gca,'xtick',0:5); xlabel(x_label); else, set(gca,'xtick',0:5); xlabel(''); end
                    if i_versns==1, ylabel([{['Pd = ',num2str(the_percnt)]},{y_label}]); else, ylabel(""); end
                    set(gca,'fontsize',12), axis tight, grid on; axis([0 5 0 20]), set(gca,'ytick',0:5:20)
                    legend('show','Location','northeast','NumColumns',2), yscale log
                    clearvars i_fig index_in x_temp y_temp x_bind y_bind fun_x fun_y TheOutFun x_label y_label


                    clearvars data the_percnt the_verson the_HCT_Cn the_beam_X the_Tumr_Z
                end
            end
        end
    end
end
clearvars i_tissue i_percent i_versns i_beam_X i_Tumr_Z i_HCT_Cn d_edges N_bins

i_fig = 3;
for i_percent = 1:length(set_of_percnt)
    for i_versns = 1:length(set_of_versns)
        figure(i_fig), subplot(3,3,(i_percent-1)*3+i_versns)
        pcolor(C,D,squeeze(set_of_Is(i_percent,i_versns,:,:))./C./D), colormap jet; shading interp; colorbar; clim([0 1.5])
        if i_percent==1, title(['version = ',num2str((i_versns))]); else, title("");end
        if i_percent==length(set_of_percnt), xlabel('HCT (%)'); else, xlabel(''); end; set(gca,'xtick',[min(set_of_HCT_Cn),mean((set_of_HCT_Cn)),max(set_of_HCT_Cn)]);
        if i_versns==1, ylabel([{['Pd = ',num2str(set_of_percnt(i_percent))]},{'d (cm)'}]); else, ylabel(""); end
        set(gca,'fontsize',12), axis tight, grid on; axis([min(set_of_HCT_Cn) max(set_of_HCT_Cn) 0 5]), set(gca,'ytick',0:1:5)
    end
end
end
