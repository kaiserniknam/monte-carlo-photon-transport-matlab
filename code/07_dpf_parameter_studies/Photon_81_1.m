function [] = Photon_81_1 ()
% Repository group: 07_dpf_parameter_studies
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% This version computes w and s as functions of d across varying mu_s and mu_a,
% based on Case 33 but using a different set of optical properties.
% analysis data: make data in a single file

clc
close all

dlta_d = 0.2; % in cm
n = 1.33; g = 0.93;
set_of_mua = (0.00:0.01:0.25);
set_of_mus = ( 001: 001:0100);
set_of_ODs = nan(length(set_of_mua),length(set_of_mus),150);
set_of__ds = nan(length(set_of_mua),length(set_of_mus),150);
set_of_DoR = nan(length(set_of_mua),length(set_of_mus));

for i_a = 1:length(set_of_mua)
    for i_s = 1:length(set_of_mus)
        mua = set_of_mua(i_a);
        mus = set_of_mus(i_s);
        [y_bind,x_bind,z_bind] = get_OD (mua,mus,dlta_d);
        set_of__ds(i_a,i_s,1:length(x_bind)) = x_bind;
        set_of_ODs(i_a,i_s,1:length(x_bind)) = -log(y_bind);
        set_of_DoR(i_a,i_s) = -log(z_bind);
        disp(['N = ',num2str(length(x_bind))])
        clearvars mua mus x_bind y_bind
    end
end
save('/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_81.mat','set_of_mua','set_of_mus','set_of__ds','set_of_ODs','set_of_DoR')
end

function [y_bind,x_bind,z_bind] = get_OD (mua,mus,dlta_d)
Lx = 29.1; Ly = 29.1; Lz = 6.0;   % computational domain size in cm
the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_81_mua_', ...
                sprintf('%.2f',mua), '_mus_', sprintf('%.2f',mus), '.mat'];
t_db = load(the_filename);
clearvars the_filename

% Removing noisy points
u = round(-log(t_db.w)./t_db.s, 4);
u_unique = unique(u);
freq = nan(size(u_unique));
for i_f = 1:length(u_unique)
    freq(i_f) = sum(u == u_unique(i_f));
end
[~, i_fx] = max(freq);
valid_u = (u == u_unique(i_fx));
t_db.x = t_db.x(valid_u);
t_db.y = t_db.y(valid_u);
t_db.z = t_db.z(valid_u);
t_db.d = t_db.d(valid_u);
t_db.s = t_db.s(valid_u);
t_db.w = t_db.w(valid_u);
t_db.c = t_db.c(valid_u);
t_db.no_of_photons = t_db.no_of_photons - sum(~valid_u);
disp(['for mua=', num2str(mua), ...
      ' & mus=', num2str(mus), ...
      ' -> # of photons=', num2str(length(u)), ...
      ' & included ', num2str(sum(valid_u)), ...
      ' photons, ~', num2str(sum(valid_u)/length(u)*100), '%'])
clearvars u u_unique freq i_f i_fx valid_u

% 1-D sorting: diffuse photons only
TheCode = 0;
idx_code = (t_db.c == TheCode);
x_raw = t_db.x(idx_code);
y_raw = t_db.y(idx_code);
z_raw = t_db.z(idx_code);
w_raw = t_db.w(idx_code);
r_raw = sqrt(x_raw.^2 + y_raw.^2 + z_raw.^2);
% Maximum possible 3D distance from center to corner
dmax = sqrt((Lx/2).^2 + (Ly/2).^2 + (Lz/2).^2);
% Make sure bin edges fully cover the data
dmax = max(dmax, max(r_raw));
d_diff_edges = 0:dlta_d:(ceil(dmax/dlta_d)*dlta_d + dlta_d);
[~,~,index_in] = histcounts(r_raw, d_diff_edges);

% Remove points outside bins, if any
valid_bin = index_in > 0;
index_in = index_in(valid_bin);
x_temp   = r_raw(valid_bin);
y_temp   = w_raw(valid_bin);
clearvars Lz dmax d_diff_edges x_raw y_raw z_raw w_raw idx_code ...
    r_raw valid_bin fun_x fun_y TheCode shell_area source_area mua mus

% I vs. d
fun_x = @mean;
fun_y = @sum;
x_bind = accumarray(index_in, x_temp, [], fun_x, nan);
y_bind = accumarray(index_in, y_temp, [], fun_y, nan);
% Normalize by annular shell area
shell_area  = pi.*dlta_d.*dlta_d + 2.*pi.*x_bind.*dlta_d;
source_area = pi.*dlta_d.*dlta_d;
y_bind = (y_bind ./ shell_area) ./ (t_db.no_of_photons ./ source_area);
z_bind = (fun_y(y_temp)./(Lx*Ly)) ./ (t_db.no_of_photons ./ source_area);
clearvars fun_x fun_y shell_area source_area x_temp y_temp dlta_d
end
