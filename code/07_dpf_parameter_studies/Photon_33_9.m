function [] = Photon_33_9 ()
% Repository group: 07_dpf_parameter_studies
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: w & s vs. d over mu_s & mu_a: more photons, more combinations
% the analysis: population stat, transmission vs. reflectance

clc
close all
set(0,'DefaultFigureWindowStyle','docked')

% optical & geometry properties
% z_air = 0.0; % the thickness of air layer
Lx = 29.1; Ly = 29.1; Lz = 05.0; % The size of the computational domain (in cm) to prevent reflections from the boundaries.

% n = 1.4; g = 0.95;
set_of_mua = (0.0:.05:0.5)+eps;
set_of_mus = (000:050:500)+eps;
L = 150; dlta = 0.17;

% arrays to save variables
R_of_s_diff = nan(length(set_of_mua),length(set_of_mus),2,L);
T_of_s_trns = nan(length(set_of_mua),length(set_of_mus),2,L);

for i_a = 1:1:length(set_of_mua)
    for i_s = 1:1:length(set_of_mus)
        % read dbase
        mua = set_of_mua(i_a);
        mus = set_of_mus(i_s);
        the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_33_mua_',sprintf('%.2f',mua),'_mus_',sprintf('%.2f',mus),'.mat'];

        [y_bind, x_bind] = get_OD(the_filename, dlta);
        R_of_s_diff(i_a,i_s,1,1:length(x_bind)) = x_bind;
        R_of_s_diff(i_a,i_s,2,1:length(y_bind)) = y_bind;
        clearvars x_bind y_bind

        [y_bind, x_bind] = get_TR(the_filename, dlta, Lz);
        T_of_s_trns(i_a,i_s,1,1:length(x_bind)) = x_bind;
        T_of_s_trns(i_a,i_s,2,1:length(y_bind)) = y_bind;
        clearvars x_bind y_bind

        disp(['mua = ',num2str(mua,'%.2f'),', mus = ',num2str(mus,'%.0f')])
        TheColor = [(i_s-1)/(length(set_of_mus)-1) 0 1-(i_s-1)/(length(set_of_mus)-1)]; % make my own colormap

        figure(i_a), subplot(1,2,1)
        semilogy(squeeze(T_of_s_trns(i_a,i_s,1,:)),squeeze(T_of_s_trns(i_a,i_s,2,:)),'LineStyle','-','LineWidth',1.5,'Color',TheColor,'DisplayName',['\mu_s = ',num2str(mus,'%.0f'),' cm^{-1}'])
        set(gca,'fontsize',18), axis square, grid on, hold on
        xlabel('separation (cm)')
        ylabel('transmittance (cm^{-2})')
        title(['transmittance vs. separation (\mu_a = ',num2str(mua,'%.2f'),' cm^{-1})'])
        legend('show','Location','northeast','NumColumns',2)
        axis([0 8 10^(-8) 10^(0)])

        figure(i_a), subplot(1,2,2)
        semilogy(squeeze(R_of_s_diff(i_a,i_s,1,:)),squeeze(R_of_s_diff(i_a,i_s,2,:)),'LineStyle','-','LineWidth',1.5,'Color',TheColor,'DisplayName',['\mu_s = ',num2str(mus,'%.0f'),' cm^{-1}'])
        set(gca,'fontsize',18), axis square, grid on, hold on
        xlabel('separation (cm)')
        ylabel('reflectance (cm^{-2})')
        title(['reflectance vs. separation (\mu_a = ',num2str(mua,'%.2f'),' cm^{-1})'])
        legend('show','Location','northeast','NumColumns',2)
        axis([0 8 10^(-8) 10^(0)])

        clearvars the_filename mua mus TheColor
    end
end
end

function [y_bind, x_bind] = get_OD(the_filename, dlta_d)
% Bin exiting photon weights as a function of diffuse source-detector distance.
Lx = 29.1; Ly = 29.1; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
t_db = load(the_filename);

% 1-D radial binning of photons exiting from the top surface.
d_diff_edges = 0:dlta_d:sqrt((Lx/2).^2 + (Ly/2).^2)+dlta_d;
idx = t_db.z <= 0;
x_temp = sqrt(t_db.x(idx).^2 + ...
              t_db.y(idx).^2);
y_temp = t_db.w(idx);
[~, ~, index_in] = histcounts(x_temp, d_diff_edges);

% Mean distance and summed detected weight in each radial bin.
x_bind = accumarray(index_in, x_temp, [], @mean, nan);
y_bind = accumarray(index_in, y_temp, [], @sum,  nan);

% Convert summed photon weight to normalized diffuse reflectance.
ring_area = pi.*dlta_d.*dlta_d + 2.*pi.*x_bind.*dlta_d;
source_area = pi.*dlta_d.*dlta_d;
source_area = 1;
y_bind = (y_bind ./ ring_area) ./ (t_db.no_of_photons ./ source_area);
x_bind = (d_diff_edges(1:end-1)+d_diff_edges(2:end-0))./2; x_bind = x_bind(1:length(y_bind)).';
end
function [y_bind, x_bind] = get_TR(the_filename, dlta_d, Lz)
% Bin exiting photon weights as a function of diffuse source-detector distance.
Lx = 29.1; Ly = 29.1; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
t_db = load(the_filename);

% 1-D radial binning of photons exiting from the bottom surface.
d_diff_edges = 0:dlta_d:sqrt((Lx/2).^2 + (Ly/2).^2)+dlta_d;
idx = t_db.z >= Lz;
x_temp = sqrt(t_db.x(idx).^2 + ...
              t_db.y(idx).^2);
y_temp = t_db.w(idx);
[~, ~, index_in] = histcounts(x_temp, d_diff_edges);

% Mean distance and summed detected weight in each radial bin.
x_bind = accumarray(index_in, x_temp, [], @mean, nan);
y_bind = accumarray(index_in, y_temp, [], @sum,  nan);

% Convert summed photon weight to normalized diffuse transmittance.
ring_area = pi.*dlta_d.*dlta_d + 2.*pi.*x_bind.*dlta_d;
source_area = pi.*dlta_d.*dlta_d;
source_area = 1;
y_bind = (y_bind ./ ring_area) ./ (t_db.no_of_photons ./ source_area);
x_bind = (d_diff_edges(1:end-1)+d_diff_edges(2:end-0))./2; x_bind = x_bind(1:length(y_bind)).';
end
