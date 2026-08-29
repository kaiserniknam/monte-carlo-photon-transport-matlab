function [] = Photon_89_1 ()
% Repository group: 11_validation_and_calibration
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Analysis of Photon_89 Monte Carlo simulations
% Computes OD as a function of source-detector separation d
% over varying mu_a, mu_s, and g.
%
% Analysis data are collected into a single MAT file.

clc
close all

% Maximum number of distance bins
N_d = 150;       % number of distance bins
dlta_d = 0.17;
d_diff_edges = 0:dlta_d:N_d*dlta_d;
d_cntrs = 0.5 .* (d_diff_edges(1:end-1) + d_diff_edges(2:end));

set_of_mua = 0.00:0.02:0.60;
set_of_mus = 1:1:60;
set_of_g   = 0.70:0.05:0.90;

set_of_ODs = nan( ...
    length(set_of_mua), ...
    length(set_of_mus), ...
    length(set_of_g), ...
    N_d);
set_of_Rat = nan( ...
    length(set_of_mua), ...
    length(set_of_mus), ...
    length(set_of_g));

for i_a = 1:length(set_of_mua)
    for i_s = 1:length(set_of_mus)
        for i_g = 1:length(set_of_g)
            mua = set_of_mua(i_a);
            mus = set_of_mus(i_s);
            g   = set_of_g(i_g);

            [y_bind,alpha] = get_OD(mua,mus,g,d_diff_edges,dlta_d);
            set_of_ODs(i_a,i_s,i_g,1:length(y_bind)) = -log(y_bind);
            set_of_Rat(i_a,i_s,i_g) = alpha;
            clearvars mua mus g y_bind alpha
        end
    end
end

save( ...
    '/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_89.mat', ...
    'set_of_mua', ...
    'set_of_mus', ...
    'set_of_g', ...
    'd_cntrs', ...
    'set_of_ODs', ...
    'set_of_Rat', ...
    '-v7.3');

end

function [y_bind,alpha] = get_OD(mua,mus,g,d_diff_edges,dlta_d)
the_filename = [ ...
    '/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/', ...
    'Photon_89_mua_',sprintf('%.4f',mua), ...
    '_mus_',sprintf('%.2f',mus), ...
    '_g_',sprintf('%.2f',g), ...
    '.mat'];
t_db = load(the_filename);
clearvars the_filename

% Removing noisy points
u = round(-log(t_db.w)./t_db.s,4);
u_unique = unique(u);
freq = nan(size(u_unique));
for i_f = 1:length(u_unique)
    freq(i_f) = sum(u == u_unique(i_f));
end
[~,i_fx] = max(freq);
valid_u = (u == u_unique(i_fx));
t_db.x = t_db.x(valid_u);
t_db.y = t_db.y(valid_u);
t_db.z = t_db.z(valid_u);
t_db.d = t_db.d(valid_u);
t_db.s = t_db.s(valid_u);
t_db.w = t_db.w(valid_u);
t_db.c = t_db.c(valid_u);
t_db.no_of_photons = ...
    t_db.no_of_photons - sum(~valid_u);
alpha = sum(valid_u)/length(u);
fprintf([ ...
    'for mua = %.2f, mus = %.0f, g = %.2f', ...
    ' -> # photons = %d, included = %d, %.2f%%\n'], ...
    mua, ...
    mus, ...
    g, ...
    length(u), ...
    sum(valid_u), ...
    sum(valid_u)/length(u)*100);
clearvars u u_unique freq i_f i_fx valid_u

% 1-D sorting: diffuse photons only
TheCode = 0;
idx_code = (t_db.c == TheCode);
x_raw = t_db.x(idx_code);
y_raw = t_db.y(idx_code);
z_raw = t_db.z(idx_code);
w_raw = t_db.w(idx_code);
r_raw = sqrt(x_raw.^2 + y_raw.^2 + z_raw.^2);
[~,~,index_in] = histcounts(r_raw,d_diff_edges);
N_d = length(d_diff_edges)-1;
clearvars Lz d_diff_edges ...
    x_raw y_raw z_raw ...
    idx_code valid_bin ...
    fun_x fun_y TheCode ...
    shell_area source_area

% I vs. d
fun_x = @mean;
fun_y = @sum;
x_bind = accumarray(index_in, r_raw, [N_d 1], fun_x, nan);
y_bind = accumarray(index_in, w_raw, [N_d 1], fun_y, nan);

% Normalize by annular shell area
shell_area = ...
    pi.*dlta_d.*dlta_d + ...
    2.*pi.*x_bind.*dlta_d;
source_area = ...
    pi.*dlta_d.*dlta_d;
y_bind = ...
    (y_bind ./ shell_area) ./ ...
    (t_db.no_of_photons ./ source_area);
clearvars fun_x fun_y ...
    shell_area source_area ...
    x_temp y_temp ...
    dlta_d
end
