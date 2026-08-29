function [] = Photon_29_4 ()
% Repository group: 04_compression_beam_and_dpf
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: s vs. d for breast tissue (compression-induced hemodynamic changes in breast tissue: beam angle, 2-D measurement -> beam at distance)
% ref. #1: Assessment of the size, position, and optical properties of breast tumors in vivo by noninvasive optical methods
% ref. #2: Quantitative Absorption and Scattering Spectra in Thick Tissues Using Broadband Diffuse Optical Spectroscopy
% ref. #3: Tumor location of the central and nipple portion is associated with impaired survival for women with breast cancer
% ref. #4: Broadband Optical Mammography: Chromophore Concentration and Hemoglobin Saturation Contrast in Breast Cancer

clc
close all

Lx = 19.1; Ly = 19.1; Lz = 10.0; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
cx = 0.00; cy = 0.00; cz = 1.50; % Tumor depth from the skin (in cm)(0.5 to 5.0 cm)
rx = 1.00; ry = 1.00; rz = 1.00; % Tumor radius within the breast (in cm)(1.0 to 2.5 cm).
code_num = 29;    % code number
beam_phi = 0;     % polar angle of beam in cm
beam_tht = 0;     % azimuthal angle of beam in cm
beam_X = 0.0;     % X deviation of beam in cm
beam_Y = 1.0;     % Y deviation of beam in cm

N_smooth = 500;
p_threshold = 0.05;
% load database
set_of_scal = [0,1:0.5:3]; % set of scales
for i_scale = 1:length(set_of_scal)
    eval(['db_',replace(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),' = load(''/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_',num2str(code_num),'_X_',num2str(beam_X),'_Y_',num2str(beam_Y),'_phi_',num2str(beam_phi),'_tht_',num2str(beam_tht),'_sclae_',sprintf('%.1f',set_of_scal(i_scale)),'.mat'');'])
end

% % plot schematics
% figure(1)
% for i_scale = 1:length(set_of_scal)
%     eval(['M_raw = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.M_raw);']);
%     subplot(4,6,i_scale)
%     if set_of_scal(i_scale) > 0
%         [xq,yq,out] = do_smooth(linspace(-Lx/2,+Lx/2,size(M_raw,1)),linspace(0,Lz,size(M_raw,3)),double(squeeze(M_raw(:,round(Ly/2/0.1),:))),5000);
%         imagesc(xq(1,:),yq(:,1),out)
%     else
%         imagesc(linspace(-Lx/2,+Lx/2,size(M_raw,1)),linspace(-Ly/2,+Ly/2,size(M_raw,2)),squeeze(M_raw(:,:,round(cz/2/0.1))).')
%     end
%     xlabel('x (cm)'), ylabel('z (cm)'), title([get_title(set_of_scal(i_scale))]),
%     set(gca,'fontsize',12), axis equal, axis tight, axis(2.*[-1 +1 0 +3/2]); clim([1 2])
%     clearvars M_raw
% end

% plot diff_OD
figure(1), N_lines = 15; [set_of_lines,my_colormap] = make_colormap(-0.5,+0.5,N_lines);
for i_scale = 3:length(set_of_scal)
    subplot(4,6,18+i_scale)
    eval(['x_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.x);']); eval(['y_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.y);']); eval(['u_sg = (db_',strrep(sprintf('%.1f',set_of_scal(i_scale)),'.','_'),'.w);']);
    eval(['x_rf = (db_',strrep(sprintf('%.1f',set_of_scal(2      )),'.','_'),'.x);']); eval(['y_rf = (db_',strrep(sprintf('%.1f',set_of_scal(2      )),'.','_'),'.y);']); eval(['u_rf = (db_',strrep(sprintf('%.1f',set_of_scal(2      )),'.','_'),'.w);']);
    Nx_bin = size(db_0_0.N,1); Ny_bin = size(db_0_0.N,2);
    [t_db,p_val,x_c,y_c] = do_compare_sum (x_sg,y_sg,u_sg,x_rf,y_rf,u_rf,Lx,Ly,Nx_bin,Ny_bin,@(x)(-log10(x)));
    t_db(p_val>p_threshold) = nan; pause(0.05)
    clearvars x_sg y_sg u_sg x_rf y_rf u_rf Nx_bin Ny_bin
    [xq,yq,out] = do_smooth(x_c,y_c,t_db,N_smooth);
    contourf(xq, yq, out, set_of_lines, 'LineColor', 'none', 'LineWidth', 0.1), clearvars t_db p_val x_c y_c xq yq out
    hold on, rectangle('Position',[-rx*sqrt(set_of_scal(i_scale)) -ry*sqrt(set_of_scal(i_scale)) 2*rx*sqrt(set_of_scal(i_scale))+eps 2*ry*sqrt(set_of_scal(i_scale))+eps],'Curvature',[1 1],'LineWidth',1,'LineStyle','-.'), hold off % breast area
    hold on, plot(beam_X,beam_Y,'h','LineWidth',1,'MarkerSize',16,'MarkerEdgeColor','k','MarkerFaceColor','w'), hold off
    if i_scale==1, ylabel('y (cm)'); else, set(gca,'ytick',[]); end; xlabel('x (cm)'),
    set(gca,'fontsize',16), axis equal, axis tight, grid on; axis(2.*[-1+beam_X/2 +1+beam_X/2 -1+beam_Y/2 +1+beam_Y/2]); colormap jet; clim([min(set_of_lines) max(set_of_lines)])
    if i_scale==length(set_of_scal), h = colorbar; ylabel(h,'OD - OD_{no-tumor}','FontSize',16); set(h,'ytick',[min(set_of_lines) mean(set_of_lines) max(set_of_lines)]); end
end
clearvars N_lines set_of_lines my_colormap h i_scale idx
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
    out = ['compressed ',sprintf('%.1f',in)];
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
