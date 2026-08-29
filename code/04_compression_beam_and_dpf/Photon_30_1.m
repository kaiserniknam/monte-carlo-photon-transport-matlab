function [] = Photon_30_1 ()
% Repository group: 04_compression_beam_and_dpf
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: s vs. d for breast tissue (compression-induced hemodynamic changes in breast tissue: beam angle, 2-D measurement -> beam at distance and tilted)
% ref. #1: Assessment of the size, position, and optical properties of breast tumors in vivo by noninvasive optical methods
% ref. #2: Quantitative Absorption and Scattering Spectra in Thick Tissues Using Broadband Diffuse Optical Spectroscopy
% ref. #3: Tumor location of the central and nipple portion is associated with impaired survival for women with breast cancer
% ref. #4: Broadband Optical Mammography: Chromophore Concentration and Hemoglobin Saturation Contrast in Breast Cancer

clc
close all

Lx = 19.1; Ly = 19.1; Lz = 10.0; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
cx = 0.00; cy = 0.00; cz = 1.50; % Tumor depth from the skin (in cm)(0.5 to 5.0 cm)
rx = 1.00; ry = 1.00; rz = 1.00; % Tumor radius within the breast (in cm)(1.0 to 2.5 cm).
code_num = 30;      % code number
beam_phi = -1.5708; % polar angle of beam in cm
beam_tht = +0.7854; % azimuthal angle of beam in cm
beam_X = 0.0;       % X deviation of beam in cm
beam_Y = 1.0;       % Y deviation of beam in cm

N_smooth = 500;
p_threshold = 1.05;
% load database
set_of_scal = [0,1:0.5:3]; % set of scales
for i_scale = 1:length(set_of_scal)
    eval(['db_',replace(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),' = load(''/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_',num2str(code_num),'_X_',num2str(beam_X),'_Y_',num2str(beam_Y),'_phi_',num2str(beam_phi),'_tht_',num2str(beam_tht),'_sclae_',sprintf('%.1f',set_of_scal(i_scale)),'.mat'');'])
end

% plot s
N_lines = 15; [set_of_lines,my_colormap] = make_colormap(0,14,N_lines);
for i_scale = 1:length(set_of_scal)
    eval(['t_db = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.s_mean);']);
    [xq,yq,out] = do_smooth(db_0_0.x_c,db_0_0.y_c,t_db,N_smooth);
    figure(1), subplot(2,3,i_scale), contourf(xq,yq,out,set_of_lines,'LineColor','none','LineWidth',0.1), clearvars t_db xq yq out
    hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_scale)) -ry*sqrt(set_of_scal(i_scale)) 2*rx*sqrt(set_of_scal(i_scale))+eps 2*ry*sqrt(set_of_scal(i_scale))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area
    hold on, plot(beam_X,beam_Y,'o','LineWidth',2,'MarkerSize',16,'MarkerEdgeColor','k','MarkerFaceColor','y'), hold off
    xlabel('x (cm)'), ylabel('y (cm)'), title([get_title(set_of_scal(i_scale))]),
    set(gca,'fontsize',16), axis equal, axis tight, grid on; h = colorbar; axis(2.*[-1+beam_X/2 +1+beam_X/2 -1+beam_Y/2 +1+beam_Y/2]); colormap jet; clim([min(set_of_lines) max(set_of_lines)])
    ylabel(h,'s (cm)','FontSize',16)
end
clearvars N_lines set_of_lines my_colormap h i_scale

% plot diff_s
N_lines = 15; [set_of_lines,my_colormap] = make_colormap(-5,+5,N_lines);
for i_scale = 2:length(set_of_scal)
    % eval(['t_db = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.s_mean.'')-(db_',strrep(sprintf('%.1f',set_of_scal(1)),'.','_'),'.s_mean.'');']);
    % eval(['[t_db,p_val,x_c,y_c] = do_difference (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),',db_',strrep(sprintf('%.1f',set_of_scal(1)),'.','_'),',Lx,Ly,''s'',0);']); t_db(p_val>p_threshold) = nan;
    eval(['x_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.x);']); eval(['y_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.y);']); eval(['u_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.s);']);
    eval(['x_rf = (db_',strrep(sprintf('%.1f',set_of_scal(1      )),'.','_'),'.x);']); eval(['y_rf = (db_',strrep(sprintf('%.1f',set_of_scal(1      )),'.','_'),'.y);']); eval(['u_rf = (db_',strrep(sprintf('%.1f',set_of_scal(1      )),'.','_'),'.s);']);
    Nx_bin = size(db_0_0.N,1); Ny_bin = size(db_0_0.N,2);
    [t_db,p_val,x_c,y_c] = do_compare_mean (x_sg,y_sg,u_sg,x_rf,y_rf,u_rf,Lx,Ly,Nx_bin,Ny_bin);
    t_db(p_val>p_threshold) = nan; pause(0.05)
    clearvars x_sg y_sg u_sg x_rf y_rf u_rf Nx_bin Ny_bin
    [xq,yq,out] = do_smooth(x_c,y_c,t_db,N_smooth);
    figure(2), subplot(2,3,i_scale), contourf(xq, yq, out, set_of_lines, 'LineColor', 'none', 'LineWidth', 0.1), clearvars t_db p_val x_c y_c xq yq out
    hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_scale)) -ry*sqrt(set_of_scal(i_scale)) 2*rx*sqrt(set_of_scal(i_scale))+eps 2*ry*sqrt(set_of_scal(i_scale))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area
    hold on, plot(beam_X,beam_Y,'o','LineWidth',2,'MarkerSize',16,'MarkerEdgeColor','k','MarkerFaceColor','y'), hold off
    xlabel('x (cm)'), ylabel('y (cm)'), title([get_title(set_of_scal(i_scale))]),
    set(gca,'fontsize',16), axis equal, axis tight, grid on; h = colorbar; axis(2.*[-1+beam_X/2 +1+beam_X/2 -1+beam_Y/2 +1+beam_Y/2]); colormap jet; clim([min(set_of_lines) max(set_of_lines)])
    ylabel(h,'s - s_{no-tumor} (cm)','FontSize',16)
end
clearvars N_lines set_of_lines my_colormap h i_scale

% plot diff_s (2)
N_lines = 15; [set_of_lines,my_colormap] = make_colormap(-5,+5,N_lines);
for i_scale = 3:length(set_of_scal)
    % eval(['t_db = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.s_mean.'')-(db_',strrep(sprintf('%.1f',set_of_scal(2)),'.','_'),'.s_mean.'');']);
    % eval(['[t_db,p_val,x_c,y_c] = do_difference (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),',db_',strrep(sprintf('%.1f',set_of_scal(2)),'.','_'),',Lx,Ly,''s'',0);']); t_db(p_val>p_threshold) = nan;
    eval(['x_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.x);']); eval(['y_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.y);']); eval(['u_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.s);']);
    eval(['x_rf = (db_',strrep(sprintf('%.1f',set_of_scal(2      )),'.','_'),'.x);']); eval(['y_rf = (db_',strrep(sprintf('%.1f',set_of_scal(2      )),'.','_'),'.y);']); eval(['u_rf = (db_',strrep(sprintf('%.1f',set_of_scal(2      )),'.','_'),'.s);']);
    Nx_bin = size(db_0_0.N,1); Ny_bin = size(db_0_0.N,2);
    [t_db,p_val,x_c,y_c] = do_compare_mean (x_sg,y_sg,u_sg,x_rf,y_rf,u_rf,Lx,Ly,Nx_bin,Ny_bin);
    t_db(p_val>p_threshold) = nan; pause(0.05)
    clearvars x_sg y_sg u_sg x_rf y_rf u_rf Nx_bin Ny_bin
    [xq,yq,out] = do_smooth(x_c,y_c,t_db,N_smooth);
    figure(3), subplot(2,3,i_scale), contourf(xq, yq, out, set_of_lines, 'LineColor', 'none', 'LineWidth', 0.1), clearvars t_db p_val x_c y_c xq yq out
    hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_scale)) -ry*sqrt(set_of_scal(i_scale)) 2*rx*sqrt(set_of_scal(i_scale))+eps 2*ry*sqrt(set_of_scal(i_scale))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area
    hold on, plot(beam_X,beam_Y,'o','LineWidth',2,'MarkerSize',16,'MarkerEdgeColor','k','MarkerFaceColor','y'), hold off
    xlabel('x (cm)'), ylabel('y (cm)'), title([get_title(set_of_scal(i_scale))]),
    set(gca,'fontsize',16), axis equal, axis tight, grid on; h = colorbar; axis(2.*[-1+beam_X/2 +1+beam_X/2 -1+beam_Y/2 +1+beam_Y/2]); colormap jet; clim([min(set_of_lines) max(set_of_lines)])
    ylabel(h,'s - s_{scale = 1} (cm)','FontSize',16)
end
clearvars N_lines set_of_lines my_colormap h i_scale idx



% plot DPF
N_lines = 15; [set_of_lines,my_colormap] = make_colormap(1,25,N_lines);
for i_scale = 1:length(set_of_scal)
    eval(['t_db = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.s_mean)./(db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.d_mean);']);
    [xq,yq,out] = do_smooth(db_0_0.x_c,db_0_0.y_c,t_db,N_smooth);
    figure(4), subplot(2,3,i_scale), contourf(xq, yq, out, set_of_lines, 'LineColor', 'none', 'LineWidth', 0.1), clearvars t_db xq yq out
    hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_scale)) -ry*sqrt(set_of_scal(i_scale)) 2*rx*sqrt(set_of_scal(i_scale))+eps 2*ry*sqrt(set_of_scal(i_scale))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area
    hold on, plot(beam_X,beam_Y,'o','LineWidth',2,'MarkerSize',16,'MarkerEdgeColor','k','MarkerFaceColor','y'), hold off
    xlabel('x (cm)'), ylabel('y (cm)'), title([get_title(set_of_scal(i_scale))]),
    set(gca,'fontsize',16), axis equal, axis tight, grid on; h = colorbar; axis(2.*[-1+beam_X/2 +1+beam_X/2 -1+beam_Y/2 +1+beam_Y/2]); colormap jet; clim([min(set_of_lines) max(set_of_lines)])
    ylabel(h,'DPF','FontSize',16)
end
clearvars N_lines set_of_lines my_colormap h i_scale idx

% plot diff_DPF
N_lines = 15; [set_of_lines,my_colormap] = make_colormap(-10,+10,N_lines);
for i_scale = 2:length(set_of_scal)
    % eval(['t_db = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.s_mean.'')./(db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.d_mean.'');']);
    % eval(['t_db = t_db - (db_',strrep(sprintf('%.1f',set_of_scal(1)),'.','_'),'.s_mean.'')./(db_',strrep(sprintf('%.1f',set_of_scal(1)),'.','_'),'.d_mean.'');']);
    % eval(['[t_db,p_val,x_c,y_c] = do_difference (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),',db_',strrep(sprintf('%.1f',set_of_scal(1)),'.','_'),',Lx,Ly,''s'',1);']); t_db(p_val>p_threshold) = nan;
    eval(['x_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.x);']); eval(['y_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.y);']); eval(['u_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.s)./(db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.d);']);
    eval(['x_rf = (db_',strrep(sprintf('%.1f',set_of_scal(1      )),'.','_'),'.x);']); eval(['y_rf = (db_',strrep(sprintf('%.1f',set_of_scal(1      )),'.','_'),'.y);']); eval(['u_rf = (db_',strrep(sprintf('%.1f',set_of_scal(1      )),'.','_'),'.s)./(db_',strrep(sprintf('%.1f',set_of_scal(1      )),'.','_'),'.d);']);
    Nx_bin = size(db_0_0.N,1); Ny_bin = size(db_0_0.N,2);
    [t_db,p_val,x_c,y_c] = do_compare_mean (x_sg,y_sg,u_sg,x_rf,y_rf,u_rf,Lx,Ly,Nx_bin,Ny_bin);
    t_db(p_val>p_threshold) = nan; pause(0.05)
    clearvars x_sg y_sg u_sg x_rf y_rf u_rf Nx_bin Ny_bin
    [xq,yq,out] = do_smooth(x_c,y_c,t_db,N_smooth);
    figure(5), subplot(2,3,i_scale), contourf(xq, yq, out, set_of_lines, 'LineColor', 'none', 'LineWidth', 0.1), clearvars t_db p_val x_c y_c xq yq out
    hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_scale)) -ry*sqrt(set_of_scal(i_scale)) 2*rx*sqrt(set_of_scal(i_scale))+eps 2*ry*sqrt(set_of_scal(i_scale))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area
    hold on, plot(beam_X,beam_Y,'o','LineWidth',2,'MarkerSize',16,'MarkerEdgeColor','k','MarkerFaceColor','y'), hold off
    xlabel('x (cm)'), ylabel('y (cm)'), title([get_title(set_of_scal(i_scale))]),
    set(gca,'fontsize',16), axis equal, axis tight, grid on; h = colorbar; axis(2.*[-1+beam_X/2 +1+beam_X/2 -1+beam_Y/2 +1+beam_Y/2]); colormap jet; clim([min(set_of_lines) max(set_of_lines)])
    ylabel(h,'DPF - DPF_{no-tumor}','FontSize',16)
end
clearvars N_lines set_of_lines my_colormap h i_scale idx

% plot diff_DPF (2)
N_lines = 15; [set_of_lines,my_colormap] = make_colormap(-10,+10,N_lines);
for i_scale = 3:length(set_of_scal)
    % eval(['t_db = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.s_mean.'')./(db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.d_mean.'');']);
    % eval(['t_db = t_db - (db_',strrep(sprintf('%.1f',set_of_scal(2)),'.','_'),'.s_mean.'')./(db_',strrep(sprintf('%.1f',set_of_scal(2)),'.','_'),'.d_mean.'');']);
    % eval(['[t_db,p_val,x_c,y_c] = do_difference (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),',db_',strrep(sprintf('%.1f',set_of_scal(2)),'.','_'),',Lx,Ly,''s'',1);']); t_db(p_val>p_threshold) = nan;
    eval(['x_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.x);']); eval(['y_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.y);']); eval(['u_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.s)./(db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.d);']);
    eval(['x_rf = (db_',strrep(sprintf('%.1f',set_of_scal(2      )),'.','_'),'.x);']); eval(['y_rf = (db_',strrep(sprintf('%.1f',set_of_scal(2      )),'.','_'),'.y);']); eval(['u_rf = (db_',strrep(sprintf('%.1f',set_of_scal(2      )),'.','_'),'.s)./(db_',strrep(sprintf('%.1f',set_of_scal(2      )),'.','_'),'.d);']);
    Nx_bin = size(db_0_0.N,1); Ny_bin = size(db_0_0.N,2);
    [t_db,p_val,x_c,y_c] = do_compare_mean (x_sg,y_sg,u_sg,x_rf,y_rf,u_rf,Lx,Ly,Nx_bin,Ny_bin);
    t_db(p_val>p_threshold) = nan; pause(0.05)
    clearvars x_sg y_sg u_sg x_rf y_rf u_rf Nx_bin Ny_bin
    [xq,yq,out] = do_smooth(x_c,y_c,t_db,N_smooth);
    figure(6), subplot(2,3,i_scale), contourf(xq, yq, out, set_of_lines, 'LineColor', 'none', 'LineWidth', 0.1), clearvars t_db p_val x_c y_c xq yq out
    hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_scale)) -ry*sqrt(set_of_scal(i_scale)) 2*rx*sqrt(set_of_scal(i_scale))+eps 2*ry*sqrt(set_of_scal(i_scale))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area
    hold on, plot(beam_X,beam_Y,'o','LineWidth',2,'MarkerSize',16,'MarkerEdgeColor','k','MarkerFaceColor','y'), hold off
    xlabel('x (cm)'), ylabel('y (cm)'), title([get_title(set_of_scal(i_scale))]),
    set(gca,'fontsize',16), axis equal, axis tight, grid on; h = colorbar; axis(2.*[-1+beam_X/2 +1+beam_X/2 -1+beam_Y/2 +1+beam_Y/2]); colormap jet; clim([min(set_of_lines) max(set_of_lines)])
    ylabel(h,'DPF - DPF_{scale = 1}','FontSize',16)
end
clearvars N_lines set_of_lines my_colormap h i_scale idx



% plot OD
N_lines = 15; [set_of_lines,my_colormap] = make_colormap(0,5,N_lines);
for i_scale = 1:length(set_of_scal)
    eval(['t_db = -log10(db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.w_mean);']);
    [xq,yq,out] = do_smooth(db_0_0.x_c,db_0_0.y_c,t_db,N_smooth);
    figure(7), subplot(2,3,i_scale), contourf(xq, yq, out, set_of_lines, 'LineColor', 'none', 'LineWidth', 0.1), clearvars t_db xq yq out
    hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_scale)) -ry*sqrt(set_of_scal(i_scale)) 2*rx*sqrt(set_of_scal(i_scale))+eps 2*ry*sqrt(set_of_scal(i_scale))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area
    hold on, plot(beam_X,beam_Y,'o','LineWidth',2,'MarkerSize',16,'MarkerEdgeColor','k','MarkerFaceColor','y'), hold off
    xlabel('x (cm)'), ylabel('y (cm)'), title([get_title(set_of_scal(i_scale))]),
    set(gca,'fontsize',16), axis equal, axis tight, grid on; h = colorbar; axis(2.*[-1+beam_X/2 +1+beam_X/2 -1+beam_Y/2 +1+beam_Y/2]); colormap jet; clim([min(set_of_lines) max(set_of_lines)])
    ylabel(h,'OD','FontSize',16)
end
clearvars N_lines set_of_lines my_colormap h i_scale idx

% plot diff_OD
N_lines = 15; [set_of_lines,my_colormap] = make_colormap(-0.05,+0.1,N_lines);
for i_scale = 2:length(set_of_scal)
    % eval(['t_db = -log10(db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.w_mean.'');']);
    % eval(['t_db = t_db + log10(db_',strrep(sprintf('%.1f',set_of_scal(1)),'.','_'),'.w_mean.'');']);
    % eval(['[t_db,p_val,x_c,y_c] = do_difference (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),',db_',strrep(sprintf('%.1f',set_of_scal(1)),'.','_'),',Lx,Ly,''w'',0);']); t_db(p_val>p_threshold) = nan;
    eval(['x_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.x);']); eval(['y_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.y);']); eval(['u_sg = ((db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.w)./(db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.N_samples));']);
    eval(['x_rf = (db_',strrep(sprintf('%.1f',set_of_scal(1      )),'.','_'),'.x);']); eval(['y_rf = (db_',strrep(sprintf('%.1f',set_of_scal(1      )),'.','_'),'.y);']); eval(['u_rf = ((db_',strrep(sprintf('%.1f',set_of_scal(1      )),'.','_'),'.w)./(db_',strrep(sprintf('%.1f',set_of_scal(1      )),'.','_'),'.N_samples));']);
    Nx_bin = size(db_0_0.N,1); Ny_bin = size(db_0_0.N,2);
    [t_db,p_val,x_c,y_c] = do_compare_sum (x_sg,y_sg,u_sg,x_rf,y_rf,u_rf,Lx,Ly,Nx_bin,Ny_bin,@(x)(-log10(x)));
    t_db(p_val>p_threshold) = nan; pause(0.05)
    clearvars x_sg y_sg u_sg x_rf y_rf u_rf Nx_bin Ny_bin
    [xq,yq,out] = do_smooth(x_c,y_c,t_db,N_smooth);
    figure(8), subplot(2,3,i_scale), contourf(xq, yq, out, set_of_lines, 'LineColor', 'none', 'LineWidth', 0.1), clearvars t_db p_val x_c y_c xq yq out
    hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_scale)) -ry*sqrt(set_of_scal(i_scale)) 2*rx*sqrt(set_of_scal(i_scale))+eps 2*ry*sqrt(set_of_scal(i_scale))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area
    hold on, plot(beam_X,beam_Y,'o','LineWidth',2,'MarkerSize',16,'MarkerEdgeColor','k','MarkerFaceColor','y'), hold off
    xlabel('x (cm)'), ylabel('y (cm)'), title([get_title(set_of_scal(i_scale))]),
    set(gca,'fontsize',16), axis equal, axis tight, grid on; h = colorbar; axis(2.*[-1+beam_X/2 +1+beam_X/2 -1+beam_Y/2 +1+beam_Y/2]); colormap(my_colormap); clim([min(set_of_lines) max(set_of_lines)])
    ylabel(h,'OD - OD_{no-tumor}','FontSize',16)
end
clearvars N_lines set_of_lines my_colormap h i_scale idx

% plot diff_OD (2)
N_lines = 15; [set_of_lines,my_colormap] = make_colormap(-0.05,+0.1,N_lines);
for i_scale = 3:length(set_of_scal)
    % eval(['t_db = -log10(db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.w_mean.'');']);
    % eval(['t_db = t_db + log10(db_',strrep(sprintf('%.1f',set_of_scal(2)),'.','_'),'.w_mean.'');']);
    % eval(['[t_db,p_val,x_c,y_c] = do_difference (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),',db_',strrep(sprintf('%.1f',set_of_scal(2)),'.','_'),',Lx,Ly,''w'',0);']); t_db(p_val>p_threshold) = nan;
    eval(['x_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.x);']); eval(['y_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.y);']); eval(['u_sg = ((db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.w)./(db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.N_samples));']);
    eval(['x_rf = (db_',strrep(sprintf('%.1f',set_of_scal(2      )),'.','_'),'.x);']); eval(['y_rf = (db_',strrep(sprintf('%.1f',set_of_scal(2      )),'.','_'),'.y);']); eval(['u_rf = ((db_',strrep(sprintf('%.1f',set_of_scal(2      )),'.','_'),'.w)./(db_',strrep(sprintf('%.1f',set_of_scal(2      )),'.','_'),'.N_samples));']);
    Nx_bin = size(db_0_0.N,1); Ny_bin = size(db_0_0.N,2);
    [t_db,p_val,x_c,y_c] = do_compare_sum (x_sg,y_sg,u_sg,x_rf,y_rf,u_rf,Lx,Ly,Nx_bin,Ny_bin,@(x)(-log10(x)));
    t_db(p_val>p_threshold) = nan; pause(0.05)
    clearvars x_sg y_sg u_sg x_rf y_rf u_rf Nx_bin Ny_bin
    [xq,yq,out] = do_smooth(x_c,y_c,t_db,N_smooth);
    figure(9), subplot(2,3,i_scale), contourf(xq, yq, out, set_of_lines, 'LineColor', 'none', 'LineWidth', 0.1), clearvars t_db p_val x_c y_c xq yq out
    hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_scale)) -ry*sqrt(set_of_scal(i_scale)) 2*rx*sqrt(set_of_scal(i_scale))+eps 2*ry*sqrt(set_of_scal(i_scale))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area
    hold on, plot(beam_X,beam_Y,'o','LineWidth',2,'MarkerSize',16,'MarkerEdgeColor','k','MarkerFaceColor','y'), hold off
    xlabel('x (cm)'), ylabel('y (cm)'), title([get_title(set_of_scal(i_scale))]),
    set(gca,'fontsize',16), axis equal, axis tight, grid on; h = colorbar; axis(2.*[-1+beam_X/2 +1+beam_X/2 -1+beam_Y/2 +1+beam_Y/2]); colormap(my_colormap); clim([min(set_of_lines) max(set_of_lines)])
    ylabel(h,'OD - OD_{scale = 1}','FontSize',16)
end
clearvars N_lines set_of_lines my_colormap h i_scale idx



% plot N
N_lines = 15; [set_of_lines,my_colormap] = make_colormap(0,5,N_lines);
for i_scale = 1:length(set_of_scal)
    eval(['t_db = log10(db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.N);']);
    [xq,yq,out] = do_smooth(db_0_0.x_c,db_0_0.y_c,t_db,N_smooth);
    figure(10), subplot(2,3,i_scale), contourf(xq,yq,out,set_of_lines,'LineColor','none','LineWidth',0.1), clearvars t_db xq yq out
    hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_scale)) -ry*sqrt(set_of_scal(i_scale)) 2*rx*sqrt(set_of_scal(i_scale))+eps 2*ry*sqrt(set_of_scal(i_scale))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area % breast area
    hold on, plot(beam_X,beam_Y,'+','LineWidth',2,'MarkerSize',16,'MarkerEdgeColor','k','MarkerFaceColor','y'), hold off
    xlabel('x (cm)'), ylabel('y (cm)'), title([get_title(set_of_scal(i_scale))]),
    set(gca,'fontsize',16), axis equal, axis tight, grid on; h = colorbar; axis(2.*[-1+beam_X/2 +1+beam_X/2 -1+beam_Y/2 +1+beam_Y/2]); colormap jet; clim([min(set_of_lines) max(set_of_lines)])
    ylabel(h,'log_{10}(N) (#)','FontSize',16)
end
clearvars N_lines set_of_lines my_colormap h i_scale

% plot diff_N
N_lines = 15; [set_of_lines,my_colormap] = make_colormap(-5,+5,N_lines);
for i_scale = 2:length(set_of_scal)
    eval(['x_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.x);']); eval(['y_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.y);']); eval(['u_sg = ones(size(db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.w));']);
    eval(['x_rf = (db_',strrep(sprintf('%.1f',set_of_scal(1      )),'.','_'),'.x);']); eval(['y_rf = (db_',strrep(sprintf('%.1f',set_of_scal(1      )),'.','_'),'.y);']); eval(['u_rf = ones(size(db_',strrep(sprintf('%.1f',set_of_scal(1      )),'.','_'),'.w));']);
    Nx_bin = size(db_0_0.N,1); Ny_bin = size(db_0_0.N,2);
    [t_db,p_val,x_c,y_c] = do_compare_sum (x_sg,y_sg,u_sg,x_rf,y_rf,u_rf,Lx,Ly,Nx_bin,Ny_bin,@(x)(x));
    t_db(p_val>p_threshold) = nan; pause(0.05)
    clearvars x_sg y_sg u_sg x_rf y_rf u_rf Nx_bin Ny_bin
    [xq,yq,out] = do_smooth(db_0_0.x_c,db_0_0.y_c,t_db,N_smooth);
    figure(11), subplot(2,3,i_scale), contourf(xq, yq, out, set_of_lines, 'LineColor', 'none', 'LineWidth', 0.1), clearvars t_db p_val x_c y_c xq yq out
    hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_scale)) -ry*sqrt(set_of_scal(i_scale)) 2*rx*sqrt(set_of_scal(i_scale))+eps 2*ry*sqrt(set_of_scal(i_scale))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area
    hold on, plot(beam_X,beam_Y,'o','LineWidth',2,'MarkerSize',16,'MarkerEdgeColor','k','MarkerFaceColor','y'), hold off
    xlabel('x (cm)'), ylabel('y (cm)'), title([get_title(set_of_scal(i_scale))]),
    set(gca,'fontsize',16), axis equal, axis tight, grid on; h = colorbar; axis(2.*[-1+beam_X/2 +1+beam_X/2 -1+beam_Y/2 +1+beam_Y/2]); colormap(my_colormap); clim([min(set_of_lines) max(set_of_lines)])
    ylabel(h,'N - N_{no-tumor} (#)','FontSize',16)
end
clearvars N_lines set_of_lines my_colormap h i_scale

% plot diff_N (2)
N_lines = 15; [set_of_lines,my_colormap] = make_colormap(-5,+5,N_lines);
for i_scale = 3:length(set_of_scal)
    eval(['x_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.x);']); eval(['y_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.y);']); eval(['u_sg = ones(size(db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.w));']);
    eval(['x_rf = (db_',strrep(sprintf('%.1f',set_of_scal(2      )),'.','_'),'.x);']); eval(['y_rf = (db_',strrep(sprintf('%.1f',set_of_scal(2      )),'.','_'),'.y);']); eval(['u_rf = ones(size(db_',strrep(sprintf('%.1f',set_of_scal(2      )),'.','_'),'.w));']);
    Nx_bin = size(db_0_0.N,1); Ny_bin = size(db_0_0.N,2);
    [t_db,p_val,x_c,y_c] = do_compare_sum (x_sg,y_sg,u_sg,x_rf,y_rf,u_rf,Lx,Ly,Nx_bin,Ny_bin,@(x)(x));
    t_db(p_val>p_threshold) = nan; pause(0.05)
    clearvars x_sg y_sg u_sg x_rf y_rf u_rf Nx_bin Ny_bin
    [xq,yq,out] = do_smooth(db_0_0.x_c,db_0_0.y_c,t_db,N_smooth);
    figure(12), subplot(2,3,i_scale), contourf(xq, yq, out, set_of_lines, 'LineColor', 'none', 'LineWidth', 0.1), clearvars t_db p_val x_c y_c xq yq out
    hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_scale)) -ry*sqrt(set_of_scal(i_scale)) 2*rx*sqrt(set_of_scal(i_scale))+eps 2*ry*sqrt(set_of_scal(i_scale))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area
    hold on, plot(beam_X,beam_Y,'o','LineWidth',2,'MarkerSize',16,'MarkerEdgeColor','k','MarkerFaceColor','y'), hold off
    xlabel('x (cm)'), ylabel('y (cm)'), title([get_title(set_of_scal(i_scale))]),
    set(gca,'fontsize',16), axis equal, axis tight, grid on; h = colorbar; axis(2.*[-1+beam_X/2 +1+beam_X/2 -1+beam_Y/2 +1+beam_Y/2]); colormap(my_colormap); clim([min(set_of_lines) max(set_of_lines)])
    ylabel(h,'N - N_{scale = 1} (#)','FontSize',16)
end
clearvars N_lines set_of_lines my_colormap h i_scale

% plot w_bar
N_lines = 15; [set_of_lines,my_colormap] = make_colormap(0,+1,N_lines);
for i_scale = 1:length(set_of_scal)
    eval(['t_db = -log10(db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.w_mean)./(db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.N);']);
    [xq,yq,out] = do_smooth(db_0_0.x_c,db_0_0.y_c,t_db,N_smooth);
    figure(13), subplot(2,3,i_scale), contourf(xq,yq,out,set_of_lines,'LineColor','none','LineWidth',0.1), clearvars t_db xq yq out
    hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_scale)) -ry*sqrt(set_of_scal(i_scale)) 2*rx*sqrt(set_of_scal(i_scale))+eps 2*ry*sqrt(set_of_scal(i_scale))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area % breast area
    hold on, plot(beam_X,beam_Y,'o','LineWidth',2,'MarkerSize',16,'MarkerEdgeColor','k','MarkerFaceColor','y'), hold off
    xlabel('x (cm)'), ylabel('y (cm)'), title([get_title(set_of_scal(i_scale))]),
    set(gca,'fontsize',16), axis equal, axis tight, grid on; h = colorbar; axis(2.*[-1+beam_X/2 +1+beam_X/2 -1+beam_Y/2 +1+beam_Y/2]); colormap jet; clim([min(set_of_lines) max(set_of_lines)])
    ylabel(h,'w_{avg} (a.u.)','FontSize',16)
end
clearvars N_lines set_of_lines my_colormap h i_scale

% plot diff_w_bar
N_lines = 15; [set_of_lines,my_colormap] = make_colormap(-0.5,+0.5,N_lines);
for i_scale = 2:length(set_of_scal)
    eval(['x_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.x);']); eval(['y_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.y);']); eval(['u_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.w);']);
    eval(['x_rf = (db_',strrep(sprintf('%.1f',set_of_scal(1      )),'.','_'),'.x);']); eval(['y_rf = (db_',strrep(sprintf('%.1f',set_of_scal(1      )),'.','_'),'.y);']); eval(['u_rf = (db_',strrep(sprintf('%.1f',set_of_scal(1      )),'.','_'),'.w);']);
    Nx_bin = size(db_0_0.N,1); Ny_bin = size(db_0_0.N,2);
    [t_db,p_val,x_c,y_c] = do_compare_mean (x_sg,y_sg,u_sg,x_rf,y_rf,u_rf,Lx,Ly,Nx_bin,Ny_bin);
    t_db(p_val>p_threshold) = nan; pause(0.05)
    clearvars x_sg y_sg u_sg x_rf y_rf u_rf Nx_bin Ny_bin
    [xq,yq,out] = do_smooth(x_c,y_c,t_db,N_smooth);
    figure(14), subplot(2,3,i_scale), contourf(xq, yq, out, set_of_lines, 'LineColor', 'none', 'LineWidth', 0.1), clearvars t_db p_val x_c y_c xq yq out
    hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_scale)) -ry*sqrt(set_of_scal(i_scale)) 2*rx*sqrt(set_of_scal(i_scale))+eps 2*ry*sqrt(set_of_scal(i_scale))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area
    hold on, plot(beam_X,beam_Y,'o','LineWidth',2,'MarkerSize',16,'MarkerEdgeColor','k','MarkerFaceColor','y'), hold off
    xlabel('x (cm)'), ylabel('y (cm)'), title([get_title(set_of_scal(i_scale))]),
    set(gca,'fontsize',16), axis equal, axis tight, grid on; h = colorbar; axis(2.*[-1+beam_X/2 +1+beam_X/2 -1+beam_Y/2 +1+beam_Y/2]); colormap(my_colormap); clim([min(set_of_lines) max(set_of_lines)])
    ylabel(h,'w - w_{avg,no-tumor} (a.u.)','FontSize',16)
end
clearvars N_lines set_of_lines my_colormap h i_scale

% plot diff_w_bar
N_lines = 15; [set_of_lines,my_colormap] = make_colormap(-0.5,+0.5,N_lines);
for i_scale = 3:length(set_of_scal)
    eval(['x_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.x);']); eval(['y_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.y);']); eval(['u_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.w);']);
    eval(['x_rf = (db_',strrep(sprintf('%.1f',set_of_scal(2      )),'.','_'),'.x);']); eval(['y_rf = (db_',strrep(sprintf('%.1f',set_of_scal(2      )),'.','_'),'.y);']); eval(['u_rf = (db_',strrep(sprintf('%.1f',set_of_scal(2      )),'.','_'),'.w);']);
    Nx_bin = size(db_0_0.N,1); Ny_bin = size(db_0_0.N,2);
    [t_db,p_val,x_c,y_c] = do_compare_mean (x_sg,y_sg,u_sg,x_rf,y_rf,u_rf,Lx,Ly,Nx_bin,Ny_bin);
    t_db(p_val>p_threshold) = nan; pause(0.05)
    clearvars x_sg y_sg u_sg x_rf y_rf u_rf Nx_bin Ny_bin
    [xq,yq,out] = do_smooth(x_c,y_c,t_db,N_smooth);
    figure(15), subplot(2,3,i_scale), contourf(xq, yq, out, set_of_lines, 'LineColor', 'none', 'LineWidth', 0.1), clearvars t_db p_val x_c y_c xq yq out
    hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_scale)) -ry*sqrt(set_of_scal(i_scale)) 2*rx*sqrt(set_of_scal(i_scale))+eps 2*ry*sqrt(set_of_scal(i_scale))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area
    hold on, plot(beam_X,beam_Y,'o','LineWidth',2,'MarkerSize',16,'MarkerEdgeColor','k','MarkerFaceColor','y'), hold off
    xlabel('x (cm)'), ylabel('y (cm)'), title([get_title(set_of_scal(i_scale))]),
    set(gca,'fontsize',16), axis equal, axis tight, grid on; h = colorbar; axis(2.*[-1+beam_X/2 +1+beam_X/2 -1+beam_Y/2 +1+beam_Y/2]); colormap(my_colormap) ; clim([min(set_of_lines) max(set_of_lines)])
    ylabel(h,'w - w_{avg,scale = 1} (a.u.)','FontSize',16)
end
clearvars N_lines set_of_lines my_colormap h i_scale



for i_fig = 15:-1:1, figure(i_fig), saveas(gcf,['Photon_',num2str(code_num),'_X_',num2str(beam_X),'_Y_',num2str(beam_Y),'_tht_',num2str(beam_tht),'_fig_',sprintf('%1.0f',i_fig),'.jpg' ],'jpeg'); end
figure(16)
for i_scale = 2:length(set_of_scal)
    eval(['M_row = db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.M_raw;']);
    subplot(2,5,i_scale-1)
    imagesc(linspace(-Lx/2,+Lx/2,size(M_row,1)),linspace(-Ly/2,+Ly/2,size(M_row,2)),squeeze(M_row(:,:,round(cz/set_of_scal(i_scale)/0.1))).')
    xlabel('x (cm)'), ylabel('y (cm)'), title([get_title(set_of_scal(i_scale))]),
    set(gca,'fontsize',12), axis equal, axis tight, axis(2.*[-1 +1 -1 +1]);
    subplot(2,5,5+i_scale-1)
    imagesc(linspace(-Lx/2,+Lx/2,size(M_row,1)),linspace(0,Lz,size(M_row,3)),squeeze(M_row(:,round(Ly/2/0.1),:)).')
    xlabel('x (cm)'), ylabel('z (cm)'), title([get_title(set_of_scal(i_scale))]),
    set(gca,'fontsize',12), axis equal, axis tight, axis(2.*[-1 +1 0 +3/2]);
end
end
function [xq,yq,out] = do_smooth(xrow,ycol,in,newpoints)
% in(isnan(in)) = 0;
% [xq,yq] = meshgrid(...
%             linspace(min(xrow), max(xrow), newpoints),...
%             linspace(min(ycol), max(ycol), newpoints)...
%           );
% out = interp2(xrow,ycol,in.',xq,yq,'cubic');
xq = xrow; yq = ycol; out = in.';
end
function [out] = get_title(in)
if in == 0
    out = ['no tumor'];
else
    out = ['tumor compressed ',sprintf('%.1f',in),' times'];
end
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

% function [v_out,p_val,x_c,y_c] = do_difference (db_sg,db_rf,Lx,Ly,the_letter,do_normalize)
% Nx_bin = size(db_sg.d_mean,1);
% Ny_bin = size(db_rf.d_mean,2);
%
% norm_rf = db_rf.d;
% norm_sg = db_sg.d;
% if do_normalize~=1
%     norm_rf = ones(size(norm_rf));
%     norm_sg = ones(size(norm_sg));
% end
%
% x_edges = linspace(-Lx/2,+Lx/2,Nx_bin+1); % # of x bins
% y_edges = linspace(-Ly/2,+Ly/2,Ny_bin+1); % # of y bins
% [N_rf,~,~,ind_x_rf,ind_y_rf] = histcounts2(db_rf.x,db_rf.y,x_edges,y_edges);
% [N_sg,~,~,ind_x_sg,ind_y_sg] = histcounts2(db_sg.x,db_sg.y,x_edges,y_edges);
% x_c = 1/2*(x_edges(1:end-1)+x_edges(2:end-0)).';
% y_c = 1/2*(y_edges(1:end-1)+y_edges(2:end-0)).';
% clearvars x_edges y_edges
%
% v_out =  nan(Nx_bin,Ny_bin);
% p_val = ones(Nx_bin,Ny_bin);
% unique_idx = unique([ind_x_rf;ind_x_sg]);
% unique_idy = unique([ind_y_rf;ind_y_sg]);
%
% for i_x = 1:length(unique_idx)
%     for i_y = 1:length(unique_idy)
%         idx_from_rf = ind_x_rf==unique_idx(i_x)&ind_y_rf==unique_idy(i_y);
%         idx_from_sg = ind_x_sg==unique_idx(i_x)&ind_y_sg==unique_idy(i_y);
%         if sum(idx_from_rf)>0&&sum(idx_from_sg)>0
%             eval(['p_val(unique_idx(i_x),unique_idy(i_y)) = ranksum(db_sg.',the_letter,'(idx_from_sg)./norm_sg(idx_from_sg),      db_rf.',the_letter,'(idx_from_rf)./norm_rf(idx_from_rf));'])
%             eval(['v_out(unique_idx(i_x),unique_idy(i_y)) = mean(   db_sg.',the_letter,'(idx_from_sg)./norm_sg(idx_from_sg))-mean(db_rf.',the_letter,'(idx_from_rf)./norm_rf(idx_from_rf));'])
%         end
%         clearvars idx_from_rf idx_from_sg
%     end
% end
% end
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
