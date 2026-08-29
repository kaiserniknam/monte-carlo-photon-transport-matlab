function [] = Photon_32_2 ()
% Repository group: 04_compression_beam_and_dpf
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: s vs. d for breast tissue (compression-induced hemodynamic changes in breast tissue: beam angle, 2-D measurement -> beam at center, tilted)
% ref. #1: Assessment of the size, position, and optical properties of breast tumors in vivo by noninvasive optical methods
% ref. #2: Quantitative Absorption and Scattering Spectra in Thick Tissues Using Broadband Diffuse Optical Spectroscopy
% ref. #3: Tumor location of the central and nipple portion is associated with impaired survival for women with breast cancer
% ref. #4: Broadband Optical Mammography: Chromophore Concentration and Hemoglobin Saturation Contrast in Breast Cancer
% analysis. P-value included

clc
close all

Lx = 19.1; Ly = 19.1; Lz = 10.0; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
cx = 0.00; cy = 0.00; cz = 1.50; % Tumor depth from the skin (in cm)(0.5 to 5.0 cm)
rx = 1.00; ry = 1.00; rz = 1.00; % Tumor radius within the breast (in cm)(1.0 to 2.5 cm).
code_num = 32;      % code number
beam_phi = -1.5708; % polar angle of beam in cm
beam_tht = 0.7854;  % azimuthal angle of beam in cm
beam_X = 0.0;       % X deviation of beam in cm
beam_Y = 0.0;       % Y deviation of beam in cm

N_smooth = 500;
p_threshold = 0.05;
% load database
set_of_scal = [0,1:0.5:3]; % set of scales
for i_scale = 1:length(set_of_scal)
    eval(['db_',replace(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),' = load(''/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_',num2str(code_num),'_X_',num2str(beam_X),'_Y_',num2str(beam_Y),'_phi_',num2str(beam_phi),'_tht_',num2str(beam_tht),'_sclae_',sprintf('%.1f',set_of_scal(i_scale)),'.mat'');'])
    eval(['rf_',replace(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),' = load(''/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_',num2str(31),      '_X_',num2str(beam_X),'_Y_',num2str(beam_Y),'_phi_',num2str(0),       '_tht_',num2str(0),       '_sclae_',sprintf('%.1f',set_of_scal(i_scale)),'.mat'');'])
end

% plot s
N_lines = 15; [set_of_lines,my_colormap] = make_colormap(-20,+20,N_lines);
for i_scale = 1:length(set_of_scal)
    % eval(['t_db = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.s_mean);']);
    % eval(['t_db = t_db - (rf_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.s_mean);']);
    % eval(['[t_db,p_val,x_c,y_c] = do_difference (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),',rf_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),',Lx,Ly,''s'',0);']); t_db(p_val>0.05) = nan;
    eval(['x_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.x);']); eval(['y_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.y);']); eval(['u_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.s);']);
    eval(['x_rf = (rf_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.x);']); eval(['y_rf = (rf_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.y);']); eval(['u_rf = (rf_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.s);']);
    Nx_bin = size(db_0_0.N,1); Ny_bin = size(db_0_0.N,2);
    [t_db,p_val,x_c,y_c] = do_compare (x_sg,y_sg,u_sg,x_rf,y_rf,u_rf,Lx,Ly,Nx_bin,Ny_bin,@(x)(mean(x)));
    t_db(p_val>p_threshold) = nan; pause(0.05)
    clearvars x_sg y_sg u_sg x_rf y_rf u_rf Nx_bin Ny_bin
    [xq,yq,out] = do_smooth(x_c,y_c,t_db,N_smooth);
    figure(1), subplot(2,3,i_scale), contourf(xq, yq, out, set_of_lines, 'LineColor', 'none', 'LineWidth', 0.1), clearvars t_db p_val x_c y_c xq yq out
    hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_scale)) -ry*sqrt(set_of_scal(i_scale)) 2*rx*sqrt(set_of_scal(i_scale))+eps 2*ry*sqrt(set_of_scal(i_scale))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area
    % hold on, plot(beam_Y,beam_X,'o','LineWidth',2,'MarkerSize',16,'MarkerEdgeColor','k','MarkerFaceColor','y'), hold off
    xlabel('x (cm)'), ylabel('y (cm)'), title([get_title(set_of_scal(i_scale))]),
    set(gca,'fontsize',16), axis equal, axis tight, grid on; h = colorbar; axis(2.*[-1+beam_X/2 +1+beam_X/2 -1+beam_Y/2 +1+beam_Y/2]); colormap("jet"); clim([min(set_of_lines) max(set_of_lines)])
    ylabel(h,'s_{tilted} - s (cm)','FontSize',16)
end
clearvars N_lines set_of_lines my_colormap h i_scale

% plot DPF
N_lines = 15; [set_of_lines,my_colormap] = make_colormap(-10,+10,N_lines);
for i_scale = 1:length(set_of_scal)
    % eval(['t_db = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.s_mean.'')./(db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.d_mean);']);
    % eval(['t_db = t_db - (rf_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.s_mean.'')./(db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.d_mean);']);
    % eval(['[t_db,p_val,x_c,y_c] = do_difference (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),',rf_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),',Lx,Ly,''s'',1);']); t_db(p_val>0.05) = nan;
    eval(['x_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.x);']); eval(['y_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.y);']); eval(['u_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.s)./(db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.d);']);
    eval(['x_rf = (rf_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.x);']); eval(['y_rf = (rf_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.y);']); eval(['u_rf = (rf_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.s)./(rf_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.d);']);
    Nx_bin = size(db_0_0.N,1); Ny_bin = size(db_0_0.N,2);
    [t_db,p_val,x_c,y_c] = do_compare (x_sg,y_sg,u_sg,x_rf,y_rf,u_rf,Lx,Ly,Nx_bin,Ny_bin,@(x)(mean(x)));
    t_db(p_val>p_threshold) = nan; pause(0.05)
    clearvars x_sg y_sg u_sg x_rf y_rf u_rf Nx_bin Ny_bin
    [xq,yq,out] = do_smooth(x_c,y_c,t_db,N_smooth);
    figure(2), subplot(2,3,i_scale), contourf(xq, yq, out, set_of_lines, 'LineColor', 'none', 'LineWidth', 0.1), clearvars t_db p_val x_c y_c xq yq out
    hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_scale)) -ry*sqrt(set_of_scal(i_scale)) 2*rx*sqrt(set_of_scal(i_scale))+eps 2*ry*sqrt(set_of_scal(i_scale))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area
    % hold on, plot(beam_Y,beam_X,'o','LineWidth',2,'MarkerSize',16,'MarkerEdgeColor','k','MarkerFaceColor','y'), hold off
    xlabel('x (cm)'), ylabel('y (cm)'), title([get_title(set_of_scal(i_scale))]),
    set(gca,'fontsize',16), axis equal, axis tight, grid on; h = colorbar; axis(2.*[-1+beam_X/2 +1+beam_X/2 -1+beam_Y/2 +1+beam_Y/2]); colormap("jet"); clim([min(set_of_lines) max(set_of_lines)])
    ylabel(h,'DPF_{tilted} - DPF','FontSize',16)
end
clearvars N_lines set_of_lines my_colormap h i_scale idx

% plot OD
N_lines = 15; [set_of_lines,my_colormap] = make_colormap(-1,+1,N_lines);
for i_scale = 1:length(set_of_scal)
    subplot(2,3,i_scale)
    % eval(['t_db = -log10(db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.w_mean);']);
    % eval(['t_db = t_db + log10(rf_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.w_mean);']);
    % eval(['[t_db,p_val,x_c,y_c] = do_difference (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),',rf_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),',Lx,Ly,''w'',0);']); t_db(p_val>0.05) = nan;
    eval(['x_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.x);']); eval(['y_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.y);']); eval(['u_sg = ((db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.w)./(db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.N_samples));']);
    eval(['x_rf = (rf_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.x);']); eval(['y_rf = (rf_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.y);']); eval(['u_rf = ((rf_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.w)./(rf_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.N_samples));']);
    Nx_bin = size(db_0_0.N,1); Ny_bin = size(db_0_0.N,2);
    [t_db,p_val,x_c,y_c] = do_compare (x_sg,y_sg,u_sg,x_rf,y_rf,u_rf,Lx,Ly,Nx_bin,Ny_bin,@(x)(-log10(sum(x))));
    t_db(p_val>p_threshold) = nan; pause(0.05)
    clearvars x_sg y_sg u_sg x_rf y_rf u_rf Nx_bin Ny_bin
    [xq,yq,out] = do_smooth(x_c,y_c,t_db,N_smooth);
    figure(3), subplot(2,3,i_scale), contourf(xq, yq, out, set_of_lines, 'LineColor', 'none', 'LineWidth', 0.1), clearvars t_db p_val x_c y_c xq yq out
    hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_scale)) -ry*sqrt(set_of_scal(i_scale)) 2*rx*sqrt(set_of_scal(i_scale))+eps 2*ry*sqrt(set_of_scal(i_scale))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area
    % hold on, plot(beam_Y,beam_X,'o','LineWidth',2,'MarkerSize',16,'MarkerEdgeColor','k','MarkerFaceColor','y'), hold off
    xlabel('x (cm)'), ylabel('y (cm)'), title([get_title(set_of_scal(i_scale))]),
    set(gca,'fontsize',16), axis equal, axis tight, grid on; h = colorbar; axis(2.*[-1+beam_X/2 +1+beam_X/2 -1+beam_Y/2 +1+beam_Y/2]); colormap("jet"); clim([min(set_of_lines) max(set_of_lines)])
    ylabel(h,'OD_{tilted} - OD','FontSize',16)
end
clearvars N_lines set_of_lines my_colormap h i_scale idx

% plot N
N_lines = 15; [set_of_lines,my_colormap] = make_colormap(-1,+1,N_lines);
for i_scale = 1:length(set_of_scal)
    eval(['x_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.x);']); eval(['y_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.y);']); eval(['u_sg = ones(size(db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.w));']);
    eval(['x_rf = (rf_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.x);']); eval(['y_rf = (rf_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.y);']); eval(['u_rf = ones(size(rf_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.w));']);
    Nx_bin = size(db_0_0.N,1); Ny_bin = size(db_0_0.N,2);
    [t_db,p_val,x_c,y_c] = do_compare (x_sg,y_sg,u_sg,x_rf,y_rf,u_rf,Lx,Ly,Nx_bin,Ny_bin,@(x)(log10(sum(x))));
    t_db(p_val>p_threshold) = nan; pause(0.05)
    clearvars x_sg y_sg u_sg x_rf y_rf u_rf Nx_bin Ny_bin
    [xq,yq,out] = do_smooth(x_c,y_c,t_db,N_smooth);
    figure(4), subplot(2,3,i_scale), contourf(xq, yq, out, set_of_lines, 'LineColor', 'none', 'LineWidth', 0.1), clearvars t_db p_val x_c y_c xq yq out
    hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_scale)) -ry*sqrt(set_of_scal(i_scale)) 2*rx*sqrt(set_of_scal(i_scale))+eps 2*ry*sqrt(set_of_scal(i_scale))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area
    % hold on, plot(beam_Y,beam_X,'o','LineWidth',2,'MarkerSize',16,'MarkerEdgeColor','k','MarkerFaceColor','y'), hold off
    xlabel('x (cm)'), ylabel('y (cm)'), title([get_title(set_of_scal(i_scale))]),
    set(gca,'fontsize',16), axis equal, axis tight, grid on; h = colorbar; axis(2.*[-1+beam_X/2 +1+beam_X/2 -1+beam_Y/2 +1+beam_Y/2]); colormap("jet"); clim([min(set_of_lines) max(set_of_lines)])
    ylabel(h,'N_{tilted} - N','FontSize',16)
end
clearvars N_lines set_of_lines my_colormap h i_scale idx
for i_fig = 4:-1:1, figure(i_fig); saveas(gcf,['Photon_',num2str(code_num),'_',num2str(i_fig),'.jpg' ],'jpeg'); end
end
function [xq,yq,out] = do_smooth(xrow,ycol,in,newpoints)
% % in(isnan(in)) = 0;
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
function [u_out,p_val,x_c,y_c] = do_compare (x_sg,y_sg,u_sg,x_rf,y_rf,u_rf,Lx,Ly,Nx_bin,Ny_bin,TheFun)
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
        p_val(unique_U(i_U,1),unique_U(i_U,2)) = ranksum(u_sg(idx_from_sg),        u_rf(idx_from_rf));
        u_out(unique_U(i_U,1),unique_U(i_U,2)) =  TheFun(u_sg(idx_from_sg))-TheFun(u_rf(idx_from_rf));
    end
    clearvars idx_from_rf idx_from_sg
end
end
