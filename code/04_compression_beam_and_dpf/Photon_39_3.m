function [] = Photon_39_3 ()
% Repository group: 04_compression_beam_and_dpf
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: s vs. d for breast tissue (compression-induced hemodynamic changes in breast tissue: 2-D measurement -> beam at several distances)
% the analysis: population stat / w_bar

clc
close all
set(0,'DefaultFigureWindowStyle','docked')

Lx = 19.1; Ly = 19.1; Lz = 10.0; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
cx = 0.00; cy = 0.00; cz = 1.50; % Tumor depth from the skin (in cm)(0.5 to 5.0 cm)
rx = 1.00; ry = 1.00; rz = 1.00; % Tumor radius within the breast (in cm)(1.0 to 2.5 cm).
beam_phi = 0;     % polar angle of beam in cm
beam_tht = 0;     % azimuthal angle of beam in cm
set_of_beam_X = [0.0,0.5,1.0,1.5]; % X deviation of beam in cm
set_of_beam_Y = [0.0,0.5,1.0,1.5]; % Y deviation of beam in cm
set_of_scal   = [0,1:0.5:3];   % set of scales
p_threshold = 0.05;
code_num = 39;

% load database
for i_X = 1:length(set_of_beam_X)
    for i_Y = 1:length(set_of_beam_Y)
        for i_s = 1:length(set_of_scal)
            beam_X = set_of_beam_X(i_X);
            beam_Y = set_of_beam_Y(i_Y);
            scal   = set_of_scal(i_s);
            eval(['db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',scal),'.','_'),' = load(''/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_',num2str(code_num),'_X_',num2str(beam_X),'_Y_',num2str(beam_Y),'_phi_0_tht_0_sclae_',sprintf('%.1f',scal),'.mat'');'])
        end
    end
end
clearvars i_X i_Y i_s beam_X beam_Y scal



% plot w_bar
N_lines = 50; [set_of_lines,my_colormap] = make_colormap(0,0.5,N_lines);
for i_s = 1:length(set_of_scal)
    figure(6*0+i_s)
    for i_X = 1:length(set_of_beam_X)
        for i_Y = 1:length(set_of_beam_Y)
            beam_X = set_of_beam_X(i_X);
            beam_Y = set_of_beam_Y(i_Y);
            scal   = set_of_scal(i_s);
            subplot(4,4,(i_Y-1).*length(set_of_beam_Y)+i_X)
            eval(['t_db = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',scal),'.','_'),'.w_mean)./(db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',scal),'.','_'),'.N).*((db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',scal),'.','_'),'.no_of_surf_photons)+(db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',scal),'.','_'),'.no_of_dpth_photons));']);
            contourf(db_0_0_0_0_0_0.x_c,db_0_0_0_0_0_0.y_c,-log10(t_db).',set_of_lines,'LineColor','none','LineWidth',0.1), clearvars t_db
            hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_s)) -ry*sqrt(set_of_scal(i_s)) 2*rx*sqrt(set_of_scal(i_s))+eps 2*ry*sqrt(set_of_scal(i_s))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area
            hold on, plot(beam_X,beam_Y,'o','LineWidth',1,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y'), hold off
            set(gca,'fontsize',12), axis equal, axis tight, grid on; h = colorbar; ylabel(h,'w (a.u.)')
            colormap(my_colormap); clim([min(set_of_lines) max(set_of_lines)])

            if i_Y==1, title(get_X_label(beam_X)); else, title("");end

            if i_Y==length(set_of_beam_Y), xlabel('x (cm)'); else, xlabel(''); end
            set(gca,'ytick',[-3+beam_Y +beam_Y +3+beam_Y]), ylim([-3+beam_Y +3+beam_Y])
            set(gca,'ytick',[-3         0      +3       ]), ylim([-3        +3       ])

            if i_X==1, ylabel([{[get_title(scal),', ',get_Y_label(beam_Y)]},{'y (cm)'}]); else, ylabel(""); end
            set(gca,'xtick',[-3+beam_X +beam_X +3+beam_X]), xlim([-3+beam_X +3+beam_X])
            set(gca,'xtick',[-3         0      +3       ]), xlim([-3        +3       ])

            clearvars t_db h beam_X beam_Y scal
        end
    end
end
clearvars N_lines set_of_lines my_colormap i_czs i_scale



% plot w_bar - w_bar_{no-tumor}
N_lines = 75; [set_of_lines,my_colormap] = make_colormap(-0.25,0.50,N_lines);
for i_s = 2:length(set_of_scal)
    figure(6*1+i_s)
    for i_X = 1:length(set_of_beam_X)
        for i_Y = 1:length(set_of_beam_Y)
            beam_X = set_of_beam_X(i_X);
            beam_Y = set_of_beam_Y(i_Y);
            scal   = set_of_scal(i_s);
            subplot(4,4,(i_Y-1).*length(set_of_beam_Y)+i_X)

            eval(['x_sg = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',scal),'.','_'),'.x);']);
            eval(['x_rf = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',0   ),'.','_'),'.x);']);
            eval(['y_sg = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',scal),'.','_'),'.y);']);
            eval(['y_rf = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',0   ),'.','_'),'.y);']);
            eval(['u_sg = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',scal),'.','_'),'.w);']);
            eval(['u_rf = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',0   ),'.','_'),'.w);']);

            Nx_bin = size(db_0_0_0_0_0_0.N,1); Ny_bin = size(db_0_0_0_0_0_0.N,2);
            [t_db,p_val,x_c,y_c] = do_compare_mean (x_sg,y_sg,u_sg,x_rf,y_rf,u_rf,Lx,Ly,Nx_bin,Ny_bin);
            t_db(p_val>p_threshold) = nan;
            clearvars x_sg y_sg u_sg x_rf y_rf u_rf Nx_bin Ny_bin

            contourf(x_c, y_c, t_db.', set_of_lines, 'LineColor', 'none', 'LineWidth', 0.1), clearvars t_db p_val x_c y_c
            hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_s)) -ry*sqrt(set_of_scal(i_s)) 2*rx*sqrt(set_of_scal(i_s))+eps 2*ry*sqrt(set_of_scal(i_s))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area
            hold on, plot(beam_X,beam_Y,'o','LineWidth',1,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y'), hold off
            set(gca,'fontsize',12), axis equal, axis tight, grid on; h = colorbar; ylabel(h,'w-w_{no-tumor} (a.u.)')
            colormap(my_colormap); clim([min(set_of_lines) max(set_of_lines)])

            if i_Y==1, title(get_X_label(beam_X)); else, title("");end

            if i_Y==length(set_of_beam_Y), xlabel('x (cm)'); else, xlabel(''); end
            set(gca,'ytick',[-3+beam_Y +beam_Y +3+beam_Y]), ylim([-3+beam_Y +3+beam_Y])
            set(gca,'ytick',[-3         0      +3       ]), ylim([-3        +3       ])

            if i_X==1, ylabel([{[get_title(scal),', ',get_Y_label(beam_Y)]},{'y (cm)'}]); else, ylabel(""); end
            set(gca,'xtick',[-3+beam_X +beam_X +3+beam_X]), xlim([-3+beam_X +3+beam_X])
            set(gca,'xtick',[-3         0      +3       ]), xlim([-3        +3       ])

            clearvars t_db h beam_X beam_Y scal
        end
    end
end
clearvars N_lines set_of_lines my_colormap i_czs i_scale



% plot w_bar - w_bar_{not-compressed}
N_lines = 75; [set_of_lines,my_colormap] = make_colormap(-0.25,0.50,N_lines);
for i_s = 3:length(set_of_scal)
    figure(6*2+i_s)
    for i_X = 1:length(set_of_beam_X)
        for i_Y = 1:length(set_of_beam_Y)
            beam_X = set_of_beam_X(i_X);
            beam_Y = set_of_beam_Y(i_Y);
            scal   = set_of_scal(i_s);
            subplot(4,4,(i_Y-1).*length(set_of_beam_Y)+i_X)

            eval(['x_sg = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',scal),'.','_'),'.x);']);
            eval(['x_rf = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',1   ),'.','_'),'.x);']);
            eval(['y_sg = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',scal),'.','_'),'.y);']);
            eval(['y_rf = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',1   ),'.','_'),'.y);']);
            eval(['u_sg = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',scal),'.','_'),'.w);']);
            eval(['u_rf = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',1   ),'.','_'),'.w);']);

            Nx_bin = size(db_0_0_0_0_0_0.N,1); Ny_bin = size(db_0_0_0_0_0_0.N,2);
            [t_db,p_val,x_c,y_c] = do_compare_mean (x_sg,y_sg,u_sg,x_rf,y_rf,u_rf,Lx,Ly,Nx_bin,Ny_bin);
            t_db(p_val>p_threshold) = nan;
            clearvars x_sg y_sg u_sg x_rf y_rf u_rf Nx_bin Ny_bin

            contourf(x_c, y_c, t_db.', set_of_lines, 'LineColor', 'none', 'LineWidth', 0.1), clearvars t_db p_val x_c y_c
            hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_s)) -ry*sqrt(set_of_scal(i_s)) 2*rx*sqrt(set_of_scal(i_s))+eps 2*ry*sqrt(set_of_scal(i_s))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area
            hold on, plot(beam_X,beam_Y,'o','LineWidth',1,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y'), hold off
            set(gca,'fontsize',12), axis equal, axis tight, grid on; h = colorbar; ylabel(h,'w-w_{not-comp} (a.u)')
            colormap(my_colormap); clim([min(set_of_lines) max(set_of_lines)])

            if i_Y==1, title(get_X_label(beam_X)); else, title("");end

            if i_Y==length(set_of_beam_Y), xlabel('x (cm)'); else, xlabel(''); end
            set(gca,'ytick',[-3+beam_Y +beam_Y +3+beam_Y]), ylim([-3+beam_Y +3+beam_Y])
            set(gca,'ytick',[-3         0      +3       ]), ylim([-3        +3       ])

            if i_X==1, ylabel([{[get_title(scal),', ',get_Y_label(beam_Y)]},{'y (cm)'}]); else, ylabel(""); end
            set(gca,'xtick',[-3+beam_X +beam_X +3+beam_X]), xlim([-3+beam_X +3+beam_X])
            set(gca,'xtick',[-3         0      +3       ]), xlim([-3        +3       ])

            clearvars t_db h beam_X beam_Y scal
        end
    end
end
clearvars N_lines set_of_lines my_colormap i_czs i_scale



for i_fig = 18:-1:1, figure(i_fig), saveas(gcf,['Photon_',num2str(code_num),'_',num2str(3),'_fig_',sprintf('%2.0f',i_fig),'.jpg' ],'jpeg'); end
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
function [u_out,p_val,x_c,y_c] = do_compare_sum  (x_sg,y_sg,u_sg,x_rf,y_rf,u_rf,Lx,Ly,Nx_bin,Ny_bin,TheFun)
x_edges = linspace(-Lx/2,+Lx/2,Nx_bin+1); % # of x bins
y_edges = linspace(-Ly/2,+Ly/2,Ny_bin+1); % # of y bins
[N_rf,~,~,ind_x_rf,ind_y_rf] = histcounts2(x_rf,y_rf,x_edges,y_edges);
[N_sg,~,~,ind_x_sg,ind_y_sg] = histcounts2(x_sg,y_sg,x_edges,y_edges);
x_c = 1/2*(x_edges(1:end-1)+x_edges(2:end-0)).';
y_c = 1/2*(y_edges(1:end-1)+y_edges(2:end-0)).';
N_min_sample = 25; % define the minimum number of datapoints (in each of two rf and sg vectors) required for statistical analysis
N_no_sample = min(round(N_min_sample*2/4),100); % determine the number of resampling iterations (bootstrapping), ensuring it does not exceed 250
prct_sample = 0.5; % specify the fraction of data points to be sampled in each resampling iteration
clearvars x_edges y_edges

u_out =  nan(Nx_bin,Ny_bin);
p_val = ones(Nx_bin,Ny_bin);
clearvars Nx_bin Ny_bin N_rf N_sg
U = [[ind_x_rf;ind_x_sg],[ind_y_rf;ind_y_sg]];
[unique_U,~,~] = unique(U,'rows'); clearvars U
for i_U = 1:size(unique_U,1)
    idx_from_rf = ind_x_rf==unique_U(i_U,1)&ind_y_rf==unique_U(i_U,2);
    idx_from_sg = ind_x_sg==unique_U(i_U,1)&ind_y_sg==unique_U(i_U,2);
    if sum(idx_from_rf)>N_min_sample&&sum(idx_from_sg)>N_min_sample
        u_rf_loop = u_rf(idx_from_rf); t_rf = nan(N_no_sample,1);
        u_sg_loop = u_sg(idx_from_sg); t_sg = nan(N_no_sample,1);
        for i_counter = 1:N_no_sample
            t_rf(i_counter) = TheFun(sum(u_rf_loop(randperm(length(u_rf_loop),round(prct_sample*length(u_rf_loop))))));
            t_sg(i_counter) = TheFun(sum(u_sg_loop(randperm(length(u_sg_loop),round(prct_sample*length(u_sg_loop))))));
        end
        p_val(unique_U(i_U,1),unique_U(i_U,2)) = ranksum(t_sg,t_rf);
        u_out(unique_U(i_U,1),unique_U(i_U,2)) = TheFun(sum(u_sg_loop))-TheFun(sum(u_rf_loop));
        clearvars u_rf_loop u_sg_loop i_counter N t_rf t_sg
    end
    clearvars idx_from_rf idx_from_sg
end
end
