function [] = Photon_33_4 ()
% Repository group: 07_dpf_parameter_studies
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: w & s vs. d over mu_s & mu_a: more photons, more combinations
% the analysis: population stat - to plot p and q for mu_a's and mu_s's


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
set_of_mua = (0.0:.05:0.5);
set_of_mus = (000:050:500);
N_bins = 250;
set_of_conv_DPF = nan(length(set_of_mua),length(set_of_mus),2);


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
        t_db.no_of_photons = t_db.no_of_photons-sum(u~=u_unique(i_fx));
        clearvars u u_unique freq i_f i_fx

        % 1-D sorting
        d_diff_edges = linspace(0, sqrt((Lx/2).^2+(Ly/2).^2      ),N_bins+1); [~,~,ind_diff] = histcounts(sqrt(t_db.x(t_db.c==0).^2+t_db.y(t_db.c==0).^2+t_db.z(t_db.c==0).^2),d_diff_edges); % d_diffuse bins
        clearvars d_trns_edges d_diff_edges


        % I vs. d
        i_fig = 1; fun_x = @mean; fun_y = @sum; TheOutFun = @(x)(-log(x)); TheLineStyle = 'none'; TheMarker = 'o'; Therloess = 0.0; blnFit = false; blnShow = false;
        index_in = ind_diff; TheCode = 0; TheLegend = 'diffuse';
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.w(t_db.c==TheCode)./t_db.no_of_photons;
        [m,b,~,~,~] = make_a_subplot(i_fig,nan,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,0.5,500,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit,TheLegend,'k',blnShow);
        set_of_conv_DPF(i_a,i_s,1) = m/mua; % a
        set_of_conv_DPF(i_a,i_s,2) = b/mua; % b
        clearvars t_db the_filename mua mus ind_diff ind_trns m b blnFit blnShow fun_x fun_y i_fig index_in TheCode TheLineStyle TheMarker TheOutFun Therloess x_temp y_temp TheLegend
    end
end
clearvars i_a i_s


set_of_conv_DPF(isinf(set_of_conv_DPF))= nan;


figure(1)
subplot(1,2,1), contourf(set_of_mua,set_of_mus,squeeze(set_of_conv_DPF(:,:,1)).'), shading flat
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
title('a'),
set(gca,'fontsize',16), axis square, axis tight,
colormap jet; colorbar; clim([0 75])
subplot(1,2,2), contourf(set_of_mua,set_of_mus,squeeze(set_of_conv_DPF(:,:,2)).'), shading flat
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
title('b'),
set(gca,'fontsize',16), axis square, axis tight,
colormap jet; colorbar; clim([0 75])


figure(2)
subplot(1,2,1), surfc(set_of_mua,set_of_mus,squeeze(set_of_conv_DPF(:,:,1)).'), shading flat
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
zlabel('a (a.u.)'), set(gca,'ztick',0:25:75)
set(gca,'fontsize',16), axis square, axis tight, colormap jet; view([145 22.5])
axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus) 0 75])


subplot(1,2,2), surfc(set_of_mua,set_of_mus,squeeze(set_of_conv_DPF(:,:,2)).'), shading flat
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
zlabel('b (a.u.)'), set(gca,'ztick',0:25:75)
set(gca,'fontsize',16), axis square, axis tight, colormap jet; view([145 22.5])
axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus) 0 75])
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
