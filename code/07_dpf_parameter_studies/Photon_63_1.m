function [] = Photon_63_1 ()
% Repository group: 07_dpf_parameter_studies
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: studies the effect of a gradient versus a uniform dye distribution
% Same as Case 62, but with an improved parameter space and cleaner code structure.

clc
close all

N_bins = 79;
N_bins = 120;
set_of_status = sort([1, 1/2, linspace(0,-1,13)],'descend');
for i_stat = 1:length(set_of_status)
    the_status = set_of_status(i_stat);
    the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_63_status_',sprintf('%.2f',the_status),'.mat'];

    [y_bind,x_bind] = get_OD (the_filename,N_bins);
    subplot(1,2,1), plot(x_bind,      (y_bind),'DisplayName',get_legend(the_status),'LineWidth',get_linewidth(the_status),'LineStyle',get_linestyle(the_status),'Color',get_color(the_status)), hold on
    subplot(1,2,2), plot(x_bind,-log10(y_bind),'DisplayName',get_legend(the_status),'LineWidth',get_linewidth(the_status),'LineStyle',get_linestyle(the_status),'Color',get_color(the_status)), hold on
    clearvars the_filename the_status
end
subplot(1,2,1)
xlabel('separation (cm)'), set(gca,'xtick',0:2:8)
ylabel('I/I_0 (a.u.)'), set(gca,'yScale','log')
title('Relative Intensity vs. Separation for Different Scenarios')
legend('show','Location','southeast','NumColumns',3)
axis([0 8 0 0.6]), axis square, set(gca,'fontsize',14)
subplot(1,2,2)
xlabel('separation (cm)'), set(gca,'xtick',0:2:8)
ylabel('OD (a.u.)'), set(gca,'ytick',0:5:10)
title('Optical Density vs. Separation for Different Scenarios')
legend('show','Location','southeast','NumColumns',3)
axis([0 8 0 10]), axis square, set(gca,'fontsize',14)
end

function [y_bind,x_bind] = get_OD (the_filename,N_bins)
Lx = 19.1; Ly = 19.1; Lz = 6.0; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
if exist(the_filename,'file')
    t_db = load(the_filename); clearvars the_filename
else
    y_bind = [];
    x_bind = [];
    keyboard
end

% removing noisy points
u = round(-log(t_db.w)./t_db.s,4); u_unique = unique(u); freq = nan(size(u_unique));
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
disp(num2str([length(u) sum(u==u_unique(i_fx)) sum(u==u_unique(i_fx))/length(u)*100]))
t_db.no_of_photons = t_db.no_of_photons-sum(u~=u_unique(i_fx));
clearvars u u_unique freq i_f i_fx

% 1-D sorting
d_diff_edges = linspace(0,sqrt((Lx/2).^2+(Ly/2).^2),N_bins+1);
idx = t_db.z_ot<=0;
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
    (            y_bind./(2.*pi.*dlta_d.*x_bind)) ./ ...
    (t_db.no_of_photons./(   pi.*dlta_d.*dlta_d));
end
function [s_legend] = get_legend (the_status)
% 0   → no dye (high scattering)
% 1/2 → layered (gradient) distribution (low → high scattering with depth)
% 1   → uniform distribution (low scattering)
% neg → randomly layered distribution (layers randomly distributed)
if     the_status == 0
    s_legend = 'no dye';
elseif the_status == 1/2
    s_legend = 'gradient';
elseif the_status == 1
    s_legend = 'uniform';
else
    s_legend = 'random';
end
end
function [out] = get_linestyle (the_status)
% 0   → no dye (high scattering)
% 1/2 → layered (gradient) distribution (low → high scattering with depth)
% 1   → uniform distribution (low scattering)
% neg → randomly layered distribution (layers randomly distributed)
if     the_status == 0
    out = '-';
elseif the_status == 1/2
    out = '-';
elseif the_status == 1
    out = '-';
else
    out = ':';
end
end
function [out] = get_linewidth (the_status)
% 0   → no dye (high scattering)
% 1/2 → layered (gradient) distribution (low → high scattering with depth)
% 1   → uniform distribution (low scattering)
% neg → randomly layered distribution (layers randomly distributed)
if     the_status == 0
    out = 2;
elseif the_status == 1/2
    out = 2;
elseif the_status == 1
    out = 2;
else
    out = 1;
end
end
function [out] = get_color (the_status)
% 0   → no dye (high scattering)
% 1/2 → layered (gradient) distribution (low → high scattering with depth)
% 1   → uniform distribution (low scattering)
% neg → randomly layered distribution (layers randomly distributed)
if     the_status == 0
    out = 'k';
elseif the_status == 1/2
    out = 'b';
elseif the_status == 1
    out = 'r';
else
    out = rand(1,3);
end
end
