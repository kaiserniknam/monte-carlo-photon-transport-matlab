function [] = Photon_81_6 ()
% Repository group: 07_dpf_parameter_studies
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: w & s vs. d over mu_s & mu_a: more photons, more combinations
% the analysis: population stat
% this version: Farrell's approx.

clc
close all
set(0,'DefaultFigureWindowStyle','docked')

dlta = 0.17; % simulation radial bin width (cm)
set_of_mua = (0.00:0.01:0.25);
set_of_mus = ( 001: 001:0100);

set_of_dist = nan(length(set_of_mua),length(set_of_mus),150);
set_of_refc = nan(length(set_of_mua),length(set_of_mus),150);
set_of_frel = nan(length(set_of_mua),length(set_of_mus),150);
rsq_of_diff = nan(length(set_of_mua),length(set_of_mus));

for i_s = 1:1:length(set_of_mus)
    figure
    for i_a = 1:1:length(set_of_mua)
        mua = set_of_mua(i_a);
        mus = set_of_mus(i_s);
        disp(['mua = ',num2str(mua,'%.2f'),', mus = ',num2str(mus,'%.0f')])
        TheColor = get_color(mua/max(set_of_mua),mus/max(set_of_mus)); % make my own colormap

        % read dbase
        the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_81_mua_',sprintf('%.2f',mua),'_mus_',sprintf('%.2f',mus),'.mat'];
        [y_bind, x_bind] = get_OD(the_filename, dlta);
        set_of_dist(1:length(x_bind)) = x_bind;
        set_of_refc(1:length(y_bind)) = y_bind;

        % ===== Farrell Eq. (15)
        g = 0.93;
        musp = (1-g)*mus;
        mutp = mua + musp;
        D = 1/(3*mutp);
        mueff = sqrt(3*mua*mutp);
        A  = 1; % refractive-index matched
        zb = 2*A*D;
        z0 = 1/mutp;
        r1 = sqrt(x_bind.^2 + z0.^2);
        r2 = sqrt(x_bind.^2 + (z0 + 2*zb).^2);
        frel = (1/(4*pi))*...
            ( ...
              z0 .* (mueff + 1./r1) .* exp(-mueff*r1) ./ r1.^2 + ...
              (z0+2*zb) .* (mueff + 1./r2) .* exp(-mueff*r2) ./ r2.^2 ...
            );
        set_of_frel(i_a,i_s,1:length(frel)) = frel;
        rsq_of_diff(i_a,i_s) = rsquared(y_bind(1:min(48,length(x_bind))),frel(1:min(48,length(x_bind))));

        % semilogy(x_bind,y_bind,'Color',TheColor,'DisplayName',['\mu_a=',num2str(mua,'%.2f'),', \mu_s=',num2str(mus,'%.0f')],                                                'Marker','none','LineStyle','-','MarkerEdgeColor','w','MarkerFaceColor',TheColor,'MarkerSize',8,'LineWidth',1.5), hold on
        % semilogy(x_bind,frel,  'Color',TheColor,'DisplayName',['\mu_a=',num2str(mua,'%.2f'),', \mu_s=',num2str(mus,'%.0f'),': r^2 ~ ',num2str(rsq_of_diff(i_a,i_s),'%.2f')],'Marker','none','LineStyle',':','MarkerEdgeColor','w','MarkerFaceColor',TheColor,'MarkerSize',8,'LineWidth',2.0), hold on
        semilogy(x_bind,y_bind,'Color',TheColor,'DisplayName',['\mu_a=',num2str(mua,'%.2f')],                                                'Marker','none','LineStyle','-','MarkerEdgeColor','w','MarkerFaceColor',TheColor,'MarkerSize',8,'LineWidth',1.5), hold on
        semilogy(x_bind,frel,  'Color',TheColor,'DisplayName',['\mu_a=',num2str(mua,'%.2f'),': r^2 ~ ',num2str(rsq_of_diff(i_a,i_s),'%.2f')],'Marker','none','LineStyle',':','MarkerEdgeColor','w','MarkerFaceColor',TheColor,'MarkerSize',8,'LineWidth',2.0), hold on
        xlabel('separation (cm)'), ylabel('R (1/cm^2)')
        set(gca,'fontsize',20), grid on, axis square, hold on
        axis([0 12 10^(-12) 1])
        legend("show",'Location','northeastoutside','Orientation','horizontal',NumColumns=2)
        set(gca, 'YScale', 'log')

        clearvars g musp mutp D mueff A zb z0 r1 r2 frel
        clearvars mua mus TheColor
        clearvars x_bind y_bind the_filename
    end
    title(['\mu_s = ',num2str(set_of_mus(i_s),'%.0f')])
    pause(10)
    close
end

figure
pcolor(set_of_mua,set_of_mus,rsq_of_diff.'), shading flat
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
title('percent of diffused photons'),
set(gca,'fontsize',20), axis square, axis tight,
colormap jet; h = colorbar; % set(gca,'colorscale','log');
ylabel(h,'r^2','FontSize',16), clearvars h
end

function [y_bind, x_bind] = get_OD(the_filename, dlta_d)
% Bin exiting photon weights as a function of diffuse source-detector distance.
Lx = 29.1; Ly = 29.1; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
t_db = load(the_filename);

% 1-D radial binning of photons exiting from the top surface.
d_diff_edges = 0:dlta_d:sqrt((Lx/2).^2 + (Ly/2).^2)+dlta_d;
idx = t_db.z <= 0;
x_temp = sqrt(t_db.x(idx).^2 + ...
              t_db.y(idx).^2 + ...
              t_db.z(idx).^2);
y_temp = t_db.w(idx);
[~, ~, index_in] = histcounts(x_temp, d_diff_edges);

% Mean distance and summed detected weight in each radial bin.
x_bind = accumarray(index_in, x_temp, [], @mean, nan);
y_bind = accumarray(index_in, y_temp, [], @sum,  nan);

% Convert summed photon weight to normalized diffuse reflectance.
ring_area = pi.*dlta_d.*dlta_d + 2.*pi.*x_bind.*dlta_d;
% source_area = pi.*dlta_d.*dlta_d;
% y_bind = (y_bind ./ ring_area) ./ (t_db.no_of_photons ./ source_area);
y_bind = (y_bind ./ ring_area) ./ (t_db.no_of_photons               );
x_bind = (d_diff_edges(1:end-1)+d_diff_edges(2:end-0))./2; x_bind = x_bind(1:length(y_bind)).';

end
function [r2] = rsquared(y, yhat)
% Compute coefficient of determination after removing invalid values.
y = y(:);
yhat = yhat(:);
idx = ~isnan(y) & ~isnan(yhat) & ~isinf(y) & ~isinf(yhat);
y = -log10(y(idx));
yhat = -log10(yhat(idx));

SS_res = sum((y - yhat).^2);
SS_tot = sum((y - mean(y)).^2);
r2 = 1 - SS_res/SS_tot;
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
