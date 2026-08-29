function [] = Photon_58_1 ()
% Repository group: 07_dpf_parameter_studies
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: w & s vs. d over mu_s & mu_a: similar to 33 but for TiO2 experiment
% analysis of data

clc
close all

n = 1.33; g = 0.93;
set_of_mua = [0.0275];
set_of_mus = (30:120)./4;
dlta_d = 0.17;
for i_a = 1:length(set_of_mua)
    for i_s = 1:length(set_of_mus)
        mua = set_of_mua(i_a);
        mus = set_of_mus(i_s);
        the_color = [(i_s-1)/(length(set_of_mus)-1) 0 1-(i_s-1)/(length(set_of_mus)-1)];
        [y_bind,x_bind] = get_OD (mua,mus,dlta_d);
        plot(x_bind,-log(y_bind),'Color',the_color,'DisplayName',['\mu_s = ',num2str(mus)]), hold on
        clearvars mua mus the_color
    end
end
xlabel('separation (cm)')
ylabel('OD (a.u.)')
grid on, axis square, hold off, axis([0 8 2 14])
set(gca,'fontsize',18), legend('show','Location','southeast','NumColumns',1)
title('OD vs source-detector separation')
end
function [y_bind,x_bind] = get_OD (mua,mus,dlta_d)
Lx = 29.1; Ly = 29.1; Lz = 06.0; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
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
