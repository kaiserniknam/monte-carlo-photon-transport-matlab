function [] = Photon_39_8 ()
% Repository group: 04_compression_beam_and_dpf
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: s vs. d for breast tissue (compression-induced hemodynamic changes in breast tissue: 2-D measurement -> beam at several distances)
% the visualization: for lab meeting & so

clc
close all
set(0,'DefaultFigureWindowStyle','docked')

rx = 1.00; ry = 1.00; rz = 1.00; % Tumor radius within the breast (in cm)(1.0 to 2.5 cm).
set_of_beam_X = [0.0,0.5,1.0,1.5]; % X deviation of beam in cm
set_of_beam_Y = [0.0,0.5,1.0,1.5]; % Y deviation of beam in cm
set_of_scal   = [0,1:0.5:3];   % set of scales
code_num = 39;
load(['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_',num2str(code_num),'.mat'])
mua_tumr = 0.0840*100;
mua_brst = 0.0300*100;

% plot DPF
for i_s = 1:length(set_of_scal)
    for i_X = 1:length(set_of_beam_X)
        for i_Y = 1:length(set_of_beam_Y)
            beam_X = set_of_beam_X(i_X);
            beam_Y = set_of_beam_Y(i_Y);
            scal   = set_of_scal(i_s);

            figure(6*0+i_s), subplot(4,4,(i_Y-1).*length(set_of_beam_Y)+i_X)
            N_lines = 55; [set_of_lines,my_colormap] = make_colormap(0,25,N_lines);
            contourf(squeeze(set_of_DPFx(i_s,i_X,i_Y,:))-beam_X, squeeze(set_of_DPFy(i_s,i_X,i_Y,:))-beam_Y, squeeze(set_of_DPFs(i_s,i_X,i_Y,:,:)).', set_of_lines, 'LineColor', 'none', 'LineWidth', 0.1), clearvars t_db p_val x_c y_c
            hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_s))-beam_X -ry*sqrt(set_of_scal(i_s))-beam_Y 2*rx*sqrt(set_of_scal(i_s))+eps 2*ry*sqrt(set_of_scal(i_s))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area
            set(gca,'fontsize',12), axis equal, axis tight, grid on; h = colorbar; ylabel(h,'OD/d (cm^{-1})')
            colormap("jet"); clim([min(set_of_lines) max(set_of_lines)])
            if i_Y==1, title(get_X_label(beam_X)); else, title("");end
            if i_Y==length(set_of_beam_Y), xlabel('x (cm)'); else, xlabel(''); end
            set(gca,'ytick',[-3         0      +3       ]), ylim([-3        +3       ])
            if i_X==1, ylabel([{[get_title(scal),', ',get_Y_label(beam_Y)]},{'y (cm)'}]); else, ylabel(""); end
            set(gca,'xtick',[-3         0      +3       ]), xlim([-3        +3       ])
            clearvars t_db h N_lines set_of_lines my_colormap

            clearvars beam_X beam_Y scal
        end
    end
end



% plot OD
for i_s = 1:length(set_of_scal)
    for i_X = 1:length(set_of_beam_X)
        for i_Y = 1:length(set_of_beam_Y)
            beam_X = set_of_beam_X(i_X);
            beam_Y = set_of_beam_Y(i_Y);
            scal   = set_of_scal(i_s);

            figure(6*1+i_s), subplot(4,4,(i_Y-1).*length(set_of_beam_Y)+i_X)
            N_lines = 55; [set_of_lines,my_colormap] = make_colormap(0,18,N_lines);
            contourf(squeeze(set_of_OD_x(i_s,i_X,i_Y,:))-beam_X, squeeze(set_of_OD_y(i_s,i_X,i_Y,:))-beam_Y, -log(squeeze(set_of_OD_s(i_s,i_X,i_Y,:,:)).'), set_of_lines, 'LineColor', 'none', 'LineWidth', 0.1), clearvars t_db p_val x_c y_c
            hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_s))-beam_X -ry*sqrt(set_of_scal(i_s))-beam_Y 2*rx*sqrt(set_of_scal(i_s))+eps 2*ry*sqrt(set_of_scal(i_s))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area
            set(gca,'fontsize',12), axis equal, axis tight, grid on; h = colorbar; ylabel(h,'OD (a.u.)')
            colormap("jet"); clim([min(set_of_lines) max(set_of_lines)])
            if i_Y==1, title(get_X_label(beam_X)); else, title("");end
            if i_Y==length(set_of_beam_Y), xlabel('x (cm)'); else, xlabel(''); end
            set(gca,'ytick',[-3         0      +3       ]), ylim([-3        +3       ])
            if i_X==1, ylabel([{[get_title(scal),', ',get_Y_label(beam_Y)]},{'y (cm)'}]); else, ylabel(""); end
            set(gca,'xtick',[-3         0      +3       ]), xlim([-3        +3       ])
            clearvars t_db h N_lines set_of_lines my_colormap

            clearvars beam_X beam_Y scal
        end
    end
end



% compare not-compressed DPF vs. true DPF
% ...







end
function [set_of_lines,my_colormap] = make_colormap(l_start,l_end,N_lines)
set_of_lines = linspace(l_start,l_end,2*N_lines+1);
my_colormap = [1, 1, 1];
if     l_end<=0
    for idx = 1:2*N_lines
        my_colormap = [...
            [1-idx/2/N_lines, 1-idx/2/N_lines, 1];...
            my_colormap];
    end
elseif l_start>=0
    for idx = 1:2*N_lines
        my_colormap = [...
            my_colormap;...
            [1, 1-idx/2/N_lines, 1-idx/2/N_lines]];
    end
else
    Np = sum(set_of_lines>0);
    Nn = sum(set_of_lines<0);
    for idx = 1:Np
        my_colormap = [...
            my_colormap;...
            [1, 1-idx/Np, 1-idx/Np]];
    end
    for idx = 1:Nn
        my_colormap = [...
            [1-idx/Nn, 1-idx/Nn, 1];...
            my_colormap];
    end
end
end
function [out] = get_title(scal)
if scal == 0
    out = 'no tumor';
else
    out = ['press. ',sprintf('%.1f',scal),' times'];
end
end
function [out] = get_X_label(beam_X)
out = ['X = ',num2str(beam_X)];
end
function [out] = get_Y_label(beam_Y)
out = ['Y = ',num2str(beam_Y)];
end

function [u_out,p_val,x_c,y_c] = do_compare_mean (x_sg,y_sg,u_sg,x_rf,y_rf,u_rf,Lx,Ly,Nx_bin,Ny_bin)
x_edges = linspace(-Lx/2,+Lx/2,Nx_bin+1); % # of x bins
y_edges = linspace(-Ly/2,+Ly/2,Ny_bin+1); % # of y bins
[N_rf,~,~,ind_x_rf,ind_y_rf] = histcounts2(x_rf,y_rf,x_edges,y_edges);
[N_sg,~,~,ind_x_sg,ind_y_sg] = histcounts2(x_sg,y_sg,x_edges,y_edges);
x_c = 1/2*(x_edges(1:end-1)+x_edges(2:end-0)).';
y_c = 1/2*(y_edges(1:end-1)+y_edges(2:end-0)).';
clearvars x_edges y_edges

u_out =  nan(Nx_bin,Ny_bin);
p_val = ones(Nx_bin,Ny_bin);
U = [[ind_x_rf;ind_x_sg],[ind_y_rf;ind_y_sg]];
[unique_U,~,~] = unique(U,'rows'); clearvars U
for i_U = 1:size(unique_U,1)
    idx_from_rf = ind_x_rf==unique_U(i_U,1)&ind_y_rf==unique_U(i_U,2);
    idx_from_sg = ind_x_sg==unique_U(i_U,1)&ind_y_sg==unique_U(i_U,2);
    if sum(idx_from_rf)>0&&sum(idx_from_sg)>0
        p_val(unique_U(i_U,1),unique_U(i_U,2)) = ranksum(u_sg(idx_from_sg),      u_rf(idx_from_rf));
        u_out(unique_U(i_U,1),unique_U(i_U,2)) =    mean(u_sg(idx_from_sg))-mean(u_rf(idx_from_rf));
    end
    clearvars idx_from_rf idx_from_sg
end
end
function [u_out,p_val,x_c,y_c] = do_compare_sum  (x_sg,y_sg,u_sg,x_rf,y_rf,u_rf,Lx,Ly,Nx_bin,Ny_bin,TheFun,n_sg,n_rf,beam_X,beam_Y)
x_edges = linspace(-Lx/2,+Lx/2,Nx_bin+1); % # of x bins
y_edges = linspace(-Ly/2,+Ly/2,Ny_bin+1); % # of y bins
[~,~,~,ind_x_rf,ind_y_rf] = histcounts2(x_rf,y_rf,x_edges,y_edges);
[~,~,~,ind_x_sg,ind_y_sg] = histcounts2(x_sg,y_sg,x_edges,y_edges);
x_c = 1/2*(x_edges(1:end-1)+x_edges(2:end-0)).';
y_c = 1/2*(y_edges(1:end-1)+y_edges(2:end-0)).';
N_min_sample = 25; % define the minimum number of datapoints (in each of two rf and sg vectors) required for statistical analysis
N_no_sample = min(round(N_min_sample*2/4),100); % determine the number of resampling iterations (bootstrapping), ensuring it does not exceed 250
prct_sample = 0.5; % specify the fraction of data points to be sampled in each resampling iteration
clearvars x_edges y_edges

u_out =  nan(Nx_bin,Ny_bin);
p_val = ones(Nx_bin,Ny_bin);
clearvars Nx_bin Ny_bin
U = [[ind_x_rf;ind_x_sg],[ind_y_rf;ind_y_sg]];
[unique_U,~,~] = unique(U,'rows'); clearvars U
for i_U = 1:size(unique_U,1)
    idx_from_rf = ind_x_rf==unique_U(i_U,1)&ind_y_rf==unique_U(i_U,2);
    idx_from_sg = ind_x_sg==unique_U(i_U,1)&ind_y_sg==unique_U(i_U,2);
    if sum(idx_from_rf)>N_min_sample&&sum(idx_from_sg)>N_min_sample
        u_rf_loop = u_rf(idx_from_rf); x_rf_loop = x_rf(idx_from_rf); y_rf_loop = y_rf(idx_from_rf); t_rf = nan(N_no_sample,1); d_rf_loop = sqrt((x_rf_loop-beam_X).^2+(y_rf_loop-beam_Y).^2);
        u_sg_loop = u_sg(idx_from_sg); x_sg_loop = x_sg(idx_from_sg); y_sg_loop = y_sg(idx_from_sg); t_sg = nan(N_no_sample,1); d_sg_loop = sqrt((x_sg_loop-beam_X).^2+(y_sg_loop-beam_Y).^2);
        for i_counter = 1:N_no_sample
            i_rf_perm = randperm(length(u_rf_loop),round(prct_sample*length(u_rf_loop)));
            i_sg_perm = randperm(length(u_sg_loop),round(prct_sample*length(u_sg_loop)));
            t_rf(i_counter) = TheFun(u_rf_loop(i_rf_perm),mean(d_rf_loop(i_rf_perm)),n_rf);
            t_sg(i_counter) = TheFun(u_sg_loop(i_sg_perm),mean(d_sg_loop(i_sg_perm)),n_sg);
            clearvars i_rf_perm i_sg_perm
        end
        p_val(unique_U(i_U,1),unique_U(i_U,2)) = ranksum(t_sg,t_rf);
        u_out(unique_U(i_U,1),unique_U(i_U,2)) = TheFun(u_sg_loop,mean(d_sg_loop),n_sg)-TheFun(u_rf_loop,mean(d_rf_loop),n_rf);
        clearvars u_rf_loop u_sg_loop x_rf_loop x_sg_loop y_rf_loop y_sg_loop d_rf_loop d_sg_loop i_counter N t_rf t_sg idx_from_rf idx_from_sg
    end
    clearvars idx_from_rf idx_from_sg
end
end
function [u_out,x_c,y_c] = do_sum  (x_sg,y_sg,u_sg,Lx,Ly,Nx_bin,Ny_bin,TheFun,n_sg,beam_X,beam_Y)
x_edges = linspace(-Lx/2,+Lx/2,Nx_bin+1); % # of x bins
y_edges = linspace(-Ly/2,+Ly/2,Ny_bin+1); % # of y bins
[~,~,~,ind_x_sg,ind_y_sg] = histcounts2(x_sg,y_sg,x_edges,y_edges);
x_c = 1/2*(x_edges(1:end-1)+x_edges(2:end-0)).';
y_c = 1/2*(y_edges(1:end-1)+y_edges(2:end-0)).';
clearvars x_edges y_edges

u_out =  nan(Nx_bin,Ny_bin);
clearvars Nx_bin Ny_bin
U = [ind_x_sg,ind_y_sg];
[unique_U,~,~] = unique(U,'rows'); clearvars U
for i_U = 1:size(unique_U,1)
    idx_from_sg = ind_x_sg==unique_U(i_U,1)&ind_y_sg==unique_U(i_U,2);
    u_sg_loop = u_sg(idx_from_sg); x_sg_loop = x_sg(idx_from_sg); y_sg_loop = y_sg(idx_from_sg); d_sg_loop = sqrt((x_sg_loop-beam_X).^2+(y_sg_loop-beam_Y).^2);
    u_out(unique_U(i_U,1),unique_U(i_U,2)) = TheFun(u_sg_loop,mean(d_sg_loop),n_sg);
    clearvars u_sg_loop x_sg_loop y_sg_loop d_sg_loop i_counter N t_sg idx_from_sg idx_from_sg
end
end
