function [] = Photon_87_3 ()
% Repository group: 11_validation_and_calibration
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% This Version: TiO2 experimental phantom
% Configuration:
%   Background: SAP hydrogel + TiO2 scattering matrix
%   No Inclusion
% Purpose:
%   Validation of DPF Paper
% Based on Version 79
% data analysis: compare 79 & 58
% data analysis: compare 1e7 & 1e6

clc
close all

dlta_d = 0.17;
mus = 55; mua = 0.0275;
set_of_musphant = [0.25, 1.00]; % set of concentrations

for i_musphant = 1:length(set_of_musphant)
    the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_87_concentration_',num2str(set_of_musphant(i_musphant),'%.2f'),'.mat'];
    the_color = [(i_musphant-1)/(length(set_of_musphant)-1) 0 1-(i_musphant-1)/(length(set_of_musphant)-1)];

    [y_bind,x_bind] = get_OD (the_filename,dlta_d);
    plot(x_bind,-log(y_bind),'LineWidth',2,'LineStyle','-','Color',the_color,'DisplayName',['C_{TiO_2} = ',num2str(set_of_musphant(i_musphant),'%.2f'),', # of photons = 10^6']), hold on
    clearvars x_bind y_bind

    [y_bind,x_bind] = get_OD_prime (mua, mus*set_of_musphant(i_musphant), dlta_d);
    plot(x_bind,-log(y_bind),'LineWidth',4,'LineStyle',':','Color',the_color,'DisplayName',['C_{TiO_2} = ',num2str(set_of_musphant(i_musphant),'%.2f'),', # of photons = 10^7']), hold on
    clearvars x_bind y_bind

    clearvars the_filename the_color
end
xlabel('separation (cm)')
ylabel('OD (a.u.)')
grid on, axis square, hold off, axis([0 8 2 12])
set(gca,'fontsize',20), legend('show','Location','southeast','NumColumns',1,'orientation','horizontal')
title('OD vs source-detector separation')
axis square
end

function [y_bind,x_bind] = get_OD       (the_filename,dlta_d)
Lx = 19.1; Ly = 19.1; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
t_db = load(the_filename); clearvars the_filename

% removing noisy points
u = round(-log(t_db.w)./t_db.s,4); u_unique = unique(u); u_unique = u_unique(~isnan(u_unique)); freq = nan(size(u_unique));
for i_f = 1:length(u_unique), freq(i_f) = sum(u==u_unique(i_f)); end
[~,i_fx] = max(freq);
t_db.x_in = t_db.x_in(u==u_unique(i_fx));
t_db.y_in = t_db.y_in(u==u_unique(i_fx));
t_db.z_in = t_db.z_in(u==u_unique(i_fx));
t_db.x_ot = t_db.x_ot(u==u_unique(i_fx));
t_db.y_ot = t_db.y_ot(u==u_unique(i_fx));
t_db.z_ot = t_db.z_ot(u==u_unique(i_fx));
t_db.s = t_db.s(u==u_unique(i_fx));
t_db.w = t_db.w(u==u_unique(i_fx));
u = u(~isnan(u));
disp(num2str([length(u) sum(u==u_unique(i_fx)) sum(u==u_unique(i_fx))/length(u)*100 u_unique(i_fx)]))
t_db.no_of_photons = t_db.no_of_photons-sum(u~=u_unique(i_fx));
clearvars u u_unique freq i_f i_fx

% 1-D sorting
d_diff_edges = 0:dlta_d:sqrt((Lx/2).^2+(Ly/2).^2)+dlta_d;
idx = t_db.z_ot<=0;
[~,~,index_in] = histcounts(sqrt( ...
    (t_db.x_in(idx)-t_db.x_ot(idx)).^2 + ...
    (t_db.y_in(idx)-t_db.y_ot(idx)).^2 + ...
    (t_db.z_in(idx)-t_db.z_ot(idx)).^2),d_diff_edges); % d_diffuse bins
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
    (            y_bind./(pi.*dlta_d.*dlta_d + 2.*pi.*x_bind.*dlta_d)) ./ ...
    (t_db.no_of_photons./(pi.*dlta_d.*dlta_d));
end
function [y_bind,x_bind] = get_OD_prime (mua,mus,dlta_d)
Lx = 29.1; Ly = 29.1; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_58_mua_',sprintf('%.4f',mua),'_mus_',sprintf('%.2f',mus),'.mat'];
t_db = load(the_filename); clearvars the_filename

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
t_db.no_of_photons = t_db.no_of_photons-sum(u~=u_unique(i_fx));
clearvars u u_unique freq i_f i_fx

% 1-D sorting
d_diff_edges = 0:dlta_d:sqrt((Lx/2).^2+(Ly/2).^2); clearvars Lx Ly Lz
[~,~,index_in] = histcounts(sqrt(t_db.x(t_db.c==0).^2+t_db.y(t_db.c==0).^2+t_db.z(t_db.c==0).^2),d_diff_edges); % d_diffuse bins
clearvars d_diff_edges Lx Ly Lz

% I vs. d
fun_x = @mean; fun_y = @sum; TheCode = 0;
x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2);
y_temp = t_db.w(t_db.c==TheCode);
x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
y_bind = accumarray(index_in,y_temp,[],fun_y,nan);
y_bind = ...
    (            y_bind./(pi.*dlta_d.*dlta_d + 2.*pi.*x_bind.*dlta_d)) ./ ...
    (t_db.no_of_photons./(pi.*dlta_d.*dlta_d));
end
