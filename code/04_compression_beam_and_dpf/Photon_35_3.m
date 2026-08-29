function [] = Photon_35_3 ()
% Repository group: 04_compression_beam_and_dpf
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: s vs. d for breast tissue (compression-induced hemodynamic changes in breast tissue: beam at center, blood is emitting)
% the analysis: average w

clc
close all
set(0,'DefaultFigureWindowStyle','docked')

Lx = 19.1; Ly = 19.1; Lz = 10.0; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
cx = 0.00; cy = 0.00;            % Tumor depth from the skin (in cm)(0.5 to 5.0 cm)
rx = 1.00; ry = 1.00; rz = 1.00; % Tumor radius within the breast (in cm)(1.0 to 2.5 cm).
beam_phi = 0;     % polar angle of beam in cm
beam_tht = 0;     % azimuthal angle of beam in cm
beam_X = 0.0;     % X deviation of beam in cm
beam_Y = 0.0;     % Y deviation of beam in cm

set_of_czs = (1.0:0.5:2.5); % Z deviation of beam in cm
set_of_scal = 1:1:3; % set of scales
p_threshold = 0.05;
code_num = 35;

% load database
db_0_0_0 = load('/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_31_X_0_Y_0_phi_0_tht_0_sclae_0.0.mat');
for i_scale = 1:length(set_of_scal)
    for i_czs = 1:length(set_of_czs)
        cz = set_of_czs(i_czs);
        eval(['db_',sprintf('%.0f',set_of_scal(i_scale)),'_',replace(sprintf('%.1f',set_of_czs(i_czs)),'.','_'),' = load(''/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_',num2str(code_num),'_cz_',sprintf('%.1f',cz),'_sclae_',sprintf('%.1f',set_of_scal(i_scale)),'.mat'');'])
    end
end
clearvars i_scale i_czs



% plot sources
figure(1),
cubeHalfSize_x = Lx/2;
cubeHalfSize_y = Ly/2;
cubeHalfSize_z = Lz/2;
for i_czs = 1:length(set_of_czs)
    cz = set_of_czs(i_czs);
    for i_scale = 1:length(set_of_scal)
        subplot(3,4,(i_scale-1).*length(set_of_czs)+i_czs),
        eval(['t_db = (db_',sprintf('%.0f',set_of_scal(i_scale)),'_',replace(sprintf('%.1f',set_of_czs(i_czs)),'.','_'),');']);
        % Draw the cube
        % Cube vertices
        vertices = [
            -cubeHalfSize_x -cubeHalfSize_y  0;
             cubeHalfSize_x -cubeHalfSize_y  0;
             cubeHalfSize_x  cubeHalfSize_y  0;
            -cubeHalfSize_x  cubeHalfSize_y  0;
            -cubeHalfSize_x -cubeHalfSize_y  -2*cubeHalfSize_z;
             cubeHalfSize_x -cubeHalfSize_y  -2*cubeHalfSize_z;
             cubeHalfSize_x  cubeHalfSize_y  -2*cubeHalfSize_z;
            -cubeHalfSize_x  cubeHalfSize_y  -2*cubeHalfSize_z;
        ];
        % Cube faces
        faces = [
            1 2 3 4; % Bottom
            5 6 7 8; % Top
            1 2 6 5; % Front
            2 3 7 6; % Right
            3 4 8 7; % Back
            4 1 5 8; % Left
        ];
        patch('Vertices', vertices, 'Faces', faces, ...
              'FaceColor', 'b', 'EdgeColor', 'black', 'LineWidth', 1,'FaceAlpha', 0.1); hold on
        % Draw sources
        plot3(t_db.st(:,1),t_db.st(:,2),-t_db.st(:,3),'r.'), hold on
        set(gca,'xtick',[round(min(t_db.x_c),1) 0 round(max(t_db.x_c),1)]), xlabel(['x-axis'])
        set(gca,'ytick',[round(min(t_db.x_c),1) 0 round(max(t_db.x_c),1)]), ylabel(['y-axis'])
        set(gca,'ztick',[-Lz 0]), zlabel(['z-axis']), title([get_title(cz,set_of_scal(i_scale))])
        set(gca,'fontsize',12), axis equal, axis tight, grid on; view([35 5])
        clearvars faces vertices t_db
    end
end
clearvars N_lines set_of_lines my_colormap h i_scale cubeHalfSize_x cubeHalfSize_y cubeHalfSize_z



% plot w_bar
figure(2), N_lines = 50; [set_of_lines,my_colormap] = make_colormap(0,0.75,N_lines);
for i_czs = 1:length(set_of_czs)
    cz = set_of_czs(i_czs);
    for i_scale = 1:length(set_of_scal)
        subplot(3,4,(i_scale-1).*length(set_of_czs)+i_czs),
        eval(['t_db = (db_',sprintf('%.0f',set_of_scal(i_scale)),'_',replace(sprintf('%.1f',set_of_czs(i_czs)),'.','_'),'.w_mean)./(db_',sprintf('%.0f',set_of_scal(i_scale)),'_',replace(sprintf('%.1f',set_of_czs(i_czs)),'.','_'),'.N).*((db_',sprintf('%.0f',set_of_scal(i_scale)),'_',replace(sprintf('%.1f',set_of_czs(i_czs)),'.','_'),'.no_of_surf_photons)+(db_',sprintf('%.0f',set_of_scal(i_scale)),'_',replace(sprintf('%.1f',set_of_czs(i_czs)),'.','_'),'.no_of_dpth_photons));']);
        contourf(db_1_1_0.x_c,db_1_1_0.y_c,t_db,set_of_lines,'LineColor','none','LineWidth',0.1), clearvars t_db
        hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_scale)) -ry*sqrt(set_of_scal(i_scale)) 2*rx*sqrt(set_of_scal(i_scale))+eps 2*ry*sqrt(set_of_scal(i_scale))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area
        xlabel('x (cm)'), ylabel('y (cm)'), title([get_title(cz,set_of_scal(i_scale))])
        set(gca,'fontsize',12), axis equal, axis tight, grid on; h = colorbar; ylabel(h,'w (a.u.)')
        colormap(my_colormap); clim([min(set_of_lines) max(set_of_lines)])
        axis(2.*[-2+beam_X/2 +2+beam_X/2 -2+beam_Y/2 +2+beam_Y/2]);
        clearvars t_db h
    end
end
clearvars N_lines set_of_lines my_colormap i_czs i_scale



% plot w - w_{tumor}
i_fig = 3; N_lines = 75; [set_of_lines,my_colormap] = make_colormap(-0.25,0.50,N_lines);
for i_czs = 1:length(set_of_czs)
    cz = set_of_czs(i_czs);
    for i_scale = 1:length(set_of_scal)
        eval(['x_sg = (db_',sprintf('%.0f',set_of_scal(i_scale)),'_',replace(sprintf('%.1f',set_of_czs(i_czs)),'.','_'),'.x);']); eval(['y_sg = (db_',sprintf('%.0f',set_of_scal(i_scale)),'_',replace(sprintf('%.1f',set_of_czs(i_czs)),'.','_'),'.y);']); eval(['u_sg = (db_',sprintf('%.0f',set_of_scal(i_scale)),'_',replace(sprintf('%.1f',set_of_czs(i_czs)),'.','_'),'.w);']);
        eval(['x_rf = (db_',sprintf('%.0f',0                   ),'_',replace(sprintf('%.1f',0.0              ),'.','_'),'.x);']); eval(['y_rf = (db_',sprintf('%.0f',0                   ),'_',replace(sprintf('%.1f',0.0              ),'.','_'),'.y);']); eval(['u_rf = (db_',sprintf('%.0f',0                   ),'_',replace(sprintf('%.1f',0.0              ),'.','_'),'.w);']);
        Nx_bin = size(db_0_0_0.N,1); Ny_bin = size(db_0_0_0.N,2);
        [t_db,p_val,x_c,y_c] = do_compare_mean (x_sg,y_sg,u_sg,x_rf,y_rf,u_rf,Lx,Ly,Nx_bin,Ny_bin);
        t_db(p_val>p_threshold) = nan;
        clearvars x_sg y_sg u_sg x_rf y_rf u_rf Nx_bin Ny_bin
        figure(i_fig), subplot(3,4,(i_scale-1).*length(set_of_czs)+i_czs), contourf(x_c, y_c, t_db, set_of_lines, 'LineColor', 'none', 'LineWidth', 0.1), clearvars t_db p_val x_c y_c
        hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_scale)) -ry*sqrt(set_of_scal(i_scale)) 2*rx*sqrt(set_of_scal(i_scale))+eps 2*ry*sqrt(set_of_scal(i_scale))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area
        xlabel('x (cm)'), ylabel('y (cm)'), title([get_title(cz,set_of_scal(i_scale))])
        set(gca,'fontsize',12), axis equal, axis tight, grid on; h = colorbar; ylabel(h,'w - w_{no-tumor}(a.u.)')
        colormap(my_colormap); clim([min(set_of_lines) max(set_of_lines)])
        axis(2.*[-2+beam_X/2 +2+beam_X/2 -2+beam_Y/2 +2+beam_Y/2]);
    end
end
clearvars N_lines set_of_lines my_colormap i_czs i_scale i_fig



% plot w - w_{cz = 1}
i_fig = 4; N_lines = 25; [set_of_lines,my_colormap] = make_colormap(-0.25,+0.0,N_lines);
for i_czs = 2:length(set_of_czs)
    cz = set_of_czs(i_czs);
    for i_scale = 1:length(set_of_scal)
        eval(['x_sg = (db_',sprintf('%.0f',set_of_scal(i_scale)),'_',replace(sprintf('%.1f',set_of_czs(i_czs)),'.','_'),'.x);']); eval(['y_sg = (db_',sprintf('%.0f',set_of_scal(i_scale)),'_',replace(sprintf('%.1f',set_of_czs(i_czs)),'.','_'),'.y);']); eval(['u_sg = (db_',sprintf('%.0f',set_of_scal(i_scale)),'_',replace(sprintf('%.1f',set_of_czs(i_czs)),'.','_'),'.w);']);
        eval(['x_rf = (db_',sprintf('%.0f',set_of_scal(i_scale)),'_',replace(sprintf('%.1f',set_of_czs(1    )),'.','_'),'.x);']); eval(['y_rf = (db_',sprintf('%.0f',set_of_scal(i_scale)),'_',replace(sprintf('%.1f',set_of_czs(1    )),'.','_'),'.y);']); eval(['u_rf = (db_',sprintf('%.0f',set_of_scal(i_scale)),'_',replace(sprintf('%.1f',set_of_czs(1    )),'.','_'),'.w);']);
        Nx_bin = size(db_0_0_0.N,1); Ny_bin = size(db_0_0_0.N,2);
        [t_db,p_val,x_c,y_c] = do_compare_mean (x_sg,y_sg,u_sg,x_rf,y_rf,u_rf,Lx,Ly,Nx_bin,Ny_bin);
        t_db(p_val>p_threshold) = nan;
        clearvars x_sg y_sg u_sg x_rf y_rf u_rf Nx_bin Ny_bin
        figure(i_fig), subplot(3,4,(i_scale-1).*length(set_of_czs)+i_czs), contourf(x_c, y_c, t_db, set_of_lines, 'LineColor', 'none', 'LineWidth', 0.1), clearvars t_db p_val x_c y_c
        hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_scale)) -ry*sqrt(set_of_scal(i_scale)) 2*rx*sqrt(set_of_scal(i_scale))+eps 2*ry*sqrt(set_of_scal(i_scale))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area
        xlabel('x (cm)'), ylabel('y (cm)'), title([get_title(cz,set_of_scal(i_scale))])
        set(gca,'fontsize',12), axis equal, axis tight, grid on; h = colorbar; ylabel(h,'w - w_{cz = 1}(a.u.)')
        colormap(my_colormap); clim([min(set_of_lines) max(set_of_lines)])
        axis(2.*[-2+beam_X/2 +2+beam_X/2 -2+beam_Y/2 +2+beam_Y/2]);
    end
end
clearvars N_lines set_of_lines my_colormap i_czs i_scale i_fig



% plot w - w_{scale = 1}
i_fig = 5; N_lines = 50; [set_of_lines,my_colormap] = make_colormap(-0.0,+0.5,N_lines);
for i_czs = 1:length(set_of_czs)
    cz = set_of_czs(i_czs);
    for i_scale = 2:length(set_of_scal)
        eval(['x_sg = (db_',sprintf('%.0f',set_of_scal(i_scale)),'_',replace(sprintf('%.1f',set_of_czs(i_czs)),'.','_'),'.x);']); eval(['y_sg = (db_',sprintf('%.0f',set_of_scal(i_scale)),'_',replace(sprintf('%.1f',set_of_czs(i_czs)),'.','_'),'.y);']); eval(['u_sg = (db_',sprintf('%.0f',set_of_scal(i_scale)),'_',replace(sprintf('%.1f',set_of_czs(i_czs)),'.','_'),'.w);']);
        eval(['x_rf = (db_',sprintf('%.0f',set_of_scal(1      )),'_',replace(sprintf('%.1f',set_of_czs(i_czs)),'.','_'),'.x);']); eval(['y_rf = (db_',sprintf('%.0f',set_of_scal(1      )),'_',replace(sprintf('%.1f',set_of_czs(i_czs)),'.','_'),'.y);']); eval(['u_rf = (db_',sprintf('%.0f',set_of_scal(1      )),'_',replace(sprintf('%.1f',set_of_czs(i_czs)),'.','_'),'.w);']);
        Nx_bin = size(db_0_0_0.N,1); Ny_bin = size(db_0_0_0.N,2);
        [t_db,p_val,x_c,y_c] = do_compare_mean (x_sg,y_sg,u_sg,x_rf,y_rf,u_rf,Lx,Ly,Nx_bin,Ny_bin);
        t_db(p_val>p_threshold) = nan;
        clearvars x_sg y_sg u_sg x_rf y_rf u_rf Nx_bin Ny_bin
        figure(i_fig), subplot(3,4,(i_scale-1).*length(set_of_czs)+i_czs), contourf(x_c, y_c, t_db, set_of_lines, 'LineColor', 'none', 'LineWidth', 0.1), clearvars t_db p_val x_c y_c
        hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_scale)) -ry*sqrt(set_of_scal(i_scale)) 2*rx*sqrt(set_of_scal(i_scale))+eps 2*ry*sqrt(set_of_scal(i_scale))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area
        xlabel('x (cm)'), ylabel('y (cm)'), title([get_title(cz,set_of_scal(i_scale))])
        set(gca,'fontsize',12), axis equal, axis tight, grid on; h = colorbar; ylabel(h,'w - w_{scale = 1}(a.u.)')
        colormap(my_colormap); clim([min(set_of_lines) max(set_of_lines)])
        axis(2.*[-2+beam_X/2 +2+beam_X/2 -2+beam_Y/2 +2+beam_Y/2]);
    end
end
clearvars N_lines set_of_lines my_colormap i_czs i_scale i_fig

for i_fig = 5:-1:1, figure(i_fig), saveas(gcf,['Photon_',num2str(code_num),'_',num2str(3),'_fig_',sprintf('%1.0f',i_fig),'.jpg' ],'jpeg'); end
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
function [out] = get_title(cz,scal)
out = ['c_z = ',num2str(cz),' cm, compressed ',sprintf('%.1f',scal),' times'];
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
N_min_sample = 25; N_no_sample = round(N_min_sample*2/4); prct_sample = 0.5;
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
