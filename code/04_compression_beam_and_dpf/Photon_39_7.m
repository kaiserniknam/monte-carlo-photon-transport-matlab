function [] = Photon_39_7 ()
% Repository group: 04_compression_beam_and_dpf
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: s vs. d for breast tissue (compression-induced hemodynamic changes in breast tissue: 2-D measurement -> beam at several distances)
% the analysis: save 2-D I for all scenarios

clc
close all
set(0,'DefaultFigureWindowStyle','docked')

Lx = 19.1; Ly = 19.1; Lz = 10.0; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
% cx = 0.00; cy = 0.00; cz = 1.50; % Tumor depth from the skin (in cm)(0.5 to 5.0 cm)
% rx = 1.00; ry = 1.00; rz = 1.00; % Tumor radius within the breast (in cm)(1.0 to 2.5 cm).
% beam_phi = 0;     % polar angle of beam in cm
% beam_tht = 0;     % azimuthal angle of beam in cm
set_of_beam_X = [0.0,0.5,1.0,1.5]; % X deviation of beam in cm
set_of_beam_Y = [0.0,0.5,1.0,1.5]; % Y deviation of beam in cm
set_of_scal   = [0,1:0.5:3];   % set of scales
% p_threshold = 0.05;
code_num = 39;


set_of_DPFs = nan(length(set_of_scal),length(set_of_beam_X),length(set_of_beam_Y),250,250);
set_of_DPFx = nan(length(set_of_scal),length(set_of_beam_X),length(set_of_beam_Y),250);
set_of_DPFy = nan(length(set_of_scal),length(set_of_beam_X),length(set_of_beam_Y),250);
set_of_OD_s = nan(length(set_of_scal),length(set_of_beam_X),length(set_of_beam_Y),250,250);
set_of_OD_x = nan(length(set_of_scal),length(set_of_beam_X),length(set_of_beam_Y),250);
set_of_OD_y = nan(length(set_of_scal),length(set_of_beam_X),length(set_of_beam_Y),250);


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



% plot DPF
for i_s = 1:length(set_of_scal)
    for i_X = 1:length(set_of_beam_X)
        for i_Y = 1:length(set_of_beam_Y)
            beam_X = set_of_beam_X(i_X);
            beam_Y = set_of_beam_Y(i_Y);
            scal   = set_of_scal(i_s);

            eval(['x_sg = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',scal),'.','_'),'.x);']);
            eval(['y_sg = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',scal),'.','_'),'.y);']);
            eval(['n_sg = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',scal),'.','_'),'.no_of_surf_photons)+(db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',scal),'.','_'),'.no_of_dpth_photons);']);
            eval(['u_sg = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',scal),'.','_'),'.w);']);

            % DPF
            Nx_bin = size(db_0_0_0_0_0_0.N,1); Ny_bin = size(db_0_0_0_0_0_0.N,2);
            [t_db,x_c,y_c] = do_sum (x_sg,y_sg,u_sg,Lx,Ly,Nx_bin,Ny_bin,@(w,d,n)(-log(sum(w)./n)./d),n_sg,beam_X,beam_Y);
            set_of_DPFs(i_s,i_X,i_Y,:,:) = t_db;
            set_of_DPFx(i_s,i_X,i_Y,:) = x_c;
            set_of_DPFy(i_s,i_X,i_Y,:) = y_c;
            clearvars x_sg y_sg n_sg u_sg x_rf y_rf n_rf u_rf Nx_bin Ny_bin
            clearvars t_db h x_c y_c

            % OD
            eval(['t_db = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',scal),'.','_'),'.w_mean);']);
            set_of_OD_s(i_s,i_X,i_Y,:,:) = t_db;
            set_of_OD_x(i_s,i_X,i_Y,:) = db_0_0_0_0_0_0.x_c;
            set_of_OD_y(i_s,i_X,i_Y,:) = db_0_0_0_0_0_0.y_c;
            clearvars t_db beam_X beam_Y scal
        end
    end
end
save('Photon_39.mat','set_of_DPFs','set_of_DPFx','set_of_DPFy','set_of_OD_s','set_of_OD_x','set_of_OD_y')
clearvars N_lines set_of_lines my_colormap i_czs i_scale
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
