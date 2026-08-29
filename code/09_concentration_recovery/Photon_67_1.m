function [] = Photon_67_1 ()
% Repository group: 09_concentration_recovery
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% This Version: TiO2 experimental phantom (India ink tumor inclusion in gelatin)
% Configuration:
%   Background: SAP hydrogel + TiO2 scattering matrix
%   Inclusion: Localized absorbing tumor (India ink embedded in gelatin)
% Purpose:
%   Validation of DPF-based absorption estimation under experimentally realistic conditions.
% Based on Version 29, adapted for TiO2 experimental dataset.

%   New tumor location at X = 1.5 cm (previously X = 0.0 cm),
%   centered between the source (X = 0 cm) and detector (X = 3.0 cm).
%
%   The inclusion is positioned within the high-sensitivity
%   (banana-shaped) photon migration region between source and detector,
%   ensuring stronger interaction with dominant photon trajectories.
%
%   This symmetric placement maximizes sensitivity to absorption changes
%   and improves robustness of DPF-based estimation.

% analysis data

clc
close all

set_of_concentration = [nan, 0.00:0.01:0.5]; % set of concentrations
set_of_ODs = nan(length(set_of_concentration),500,2);
for i_concentration = 1:length(set_of_concentration)
    the_concentration = set_of_concentration(i_concentration);
    the_color = [(i_concentration-1)./(length(set_of_concentration)-1) 0 1-(i_concentration-1)./(length(set_of_concentration)-1)];

    the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_67_concentration_',num2str(the_concentration),'.mat' ];
    [y_bind,x_bind] = get_OD (the_filename,0.17);
    set_of_ODs(i_concentration,1:length(x_bind),1) = x_bind;
    set_of_ODs(i_concentration,1:length(y_bind),2) = y_bind;

    figure(1), subplot(1,2,1), plot(x_bind,-log(y_bind),'LineStyle','-','LineWidth',2,'Color',the_color,'DisplayName',num2str(the_concentration)), hold on
    xlabel('d (cm)'), ylabel('OD (a.u.)'), set(gca,'fontsize',16)
    title('Effect of Inclusion Concentration on Optical Density'), grid off, axis tight square
    clearvars the_concentration the_color the_filename y_bind x_bind
end

dims = axis;
set_of_ODs_at_3cm = nan(length(set_of_concentration),1);
for i_concentration = 1:length(set_of_concentration)
    [~,idx] = min(abs(set_of_ODs(i_concentration,:,1)-3));
    figure(1), subplot(1,2,1), plot([set_of_ODs(i_concentration,idx,1) set_of_ODs(i_concentration,idx,1)],[dims(3) dims(4)],'LineStyle','-.','LineWidth',2,'Color','k','HandleVisibility','off'),
    set_of_ODs_at_3cm(i_concentration) = set_of_ODs(i_concentration,idx,2);
    clearvars idx
end
figure(1), subplot(1,2,2), plot(set_of_concentration,-log(set_of_ODs_at_3cm),'LineStyle','-','LineWidth',2,'Color','k'),
xlabel('concentration (%v/v)'), ylabel('OD (a.u.)'), set(gca,'fontsize',16)
title('Effect of Inclusion Concentration on Optical Density at SDS = 3 cm'), grid on, axis square
end

function [y_bind,x_bind] = get_OD (the_filename,dlta_d)
Lx = 19.1; Ly = 19.1; Lz = 6.0;  % The size of the computational domain (in cm) to prevent reflections from the boundaries.
t_db = load(the_filename); clearvars the_filename

% % removing noisy points
% u = round(-log(t_db.w)./t_db.s,4); u_unique = unique(u); freq = nan(size(u_unique));
% for i_f = 1:length(u_unique), freq(i_f) = sum(u==u_unique(i_f)); end
% [~,i_fx] = max(freq);
% t_db.x_in = t_db.x_in(u==u_unique(i_fx));
% t_db.y_in = t_db.y_in(u==u_unique(i_fx));
% t_db.z_in = t_db.z_in(u==u_unique(i_fx));
% t_db.x_ot = t_db.x_ot(u==u_unique(i_fx));
% t_db.y_ot = t_db.y_ot(u==u_unique(i_fx));
% t_db.z_ot = t_db.z_ot(u==u_unique(i_fx));
% t_db.s = t_db.s(u==u_unique(i_fx));
% t_db.w = t_db.w(u==u_unique(i_fx));
% disp(num2str([length(u) sum(u==u_unique(i_fx)) sum(u==u_unique(i_fx))/length(u)*100]))
% t_db.no_of_photons = t_db.no_of_photons-sum(u~=u_unique(i_fx));
% clearvars u u_unique freq i_f i_fx

% 1-D sorting
d_diff_edges = 0:dlta_d:sqrt((Lx/2).^2+(Ly/2).^2);
idx = t_db.z_ot<=0&-dlta_d<=t_db.y_ot&t_db.y_ot<=dlta_d&0<=t_db.x_ot;
[~,~,index_in] = histcounts(sqrt( ...
    (t_db.x_in(idx)-t_db.x_ot(idx)).^2 + ...
    (t_db.y_in(idx)-t_db.y_ot(idx)).^2 + ...
    (t_db.z_in(idx)-t_db.z_ot(idx)).^2),d_diff_edges); % d_diffuse bins
dlta_d = mean(diff(d_diff_edges));
clearvars d_diff_edges Lx Ly Lz

% I vs. d
fun_x = @mean; fun_y = @sum;
x_temp = sqrt( ...
    (t_db.x_in(idx)-t_db.x_ot(idx)).^2 + ...
    (t_db.y_in(idx)-t_db.y_ot(idx)).^2 + ...
    (t_db.z_in(idx)-t_db.z_ot(idx)).^2);
y_temp = t_db.w(idx);
x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
y_bind = accumarray(index_in,y_temp,[],fun_y,nan);
y_bind = ...
    (            y_bind./(dlta_d.*dlta_d)) ./ ...
    (t_db.no_of_photons./(dlta_d.*dlta_d));
end
