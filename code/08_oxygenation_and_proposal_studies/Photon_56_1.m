function [] = Photon_56_1 ()
% Repository group: 08_oxygenation_and_proposal_studies
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor and Signal on DPF/PPF Simulations
% This version: Comparison of source–detector separation (s vs. d) in breast tissue
% Context: Includes an embedded tumor with varying depths and oxygenation levels
% Dr. Das’s Grant Proposal - The same as 55 but simple plot
% Notes, simulations, and related analysis for proposal development
% ref: Monte Carlo investigation of the effect of blood volume and oxygen saturation on optical path in reflectance pulse oximetry
% S Chatterjee et. al. 2016 Biomed. Phys. Eng. Express 2 065018
% Population Analysis: Optical Density (OD) and Δ (Depth Difference) Results

clc
close all
code_num = 56;
set(0,'DefaultFigureWindowStyle','docked')

% z_air = 0.0; % the thickness of air layer
% Optical properties of breast and tumor [1-4]
Lx = 19.1; Ly = 19.1; Lz = 10.0; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
dx = 2.00; cy = 0.00;            % Tumor depth from the skin (in cm)
% rx = 0.65; ry = 0.65; rz = 0.65; % Tumor radius within the breast (in cm)
% beam_phi = 0;     % polar angle of beam in cm
% beam_tht = 0;     % azimuthal angle of beam in cm

set_of_beam_X = -2:2:+2;               % X deviation of beam in cm
set_of_beam_Y = -2:2:+2;               % Y deviation of beam in cm
set_of_czs = [1, 2];                   % tumor depth in cm
set_of_wvl = [660 940];                % wavelength

iBeamX = 2;
iBeamY = 2;
N_bin = 101;
set_of_data = nan(length(set_of_wvl),length(set_of_czs),N_bin,N_bin);

% load database
for i_BeamX = iBeamX %1:length(set_of_beam_X)
    for i_BeamY = iBeamY %1:length(set_of_beam_Y)
        for i_czs = 1:length(set_of_czs)
            for i_wvl = 1:length(set_of_wvl)
                beam_X = set_of_beam_X(i_BeamX);
                beam_Y = set_of_beam_Y(i_BeamY);
                cz = set_of_czs(i_czs);
                wvlngth = set_of_wvl(i_wvl);
                eval(['t_db = load(''/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_',num2str(code_num),'_X_',num2str(beam_X),'_Y_',num2str(beam_Y),'_cz_',num2str(cz),'_wvlngth_',num2str(wvlngth),'.mat'');'])
                [x_c,y_c,set_of_data(i_wvl,i_czs,:,:)] = make_2d(t_db,Lx,Ly,N_bin);
                clearvars beam_X beam_Y cz wvlngth
            end
        end
    end
end
clearvars i_BeamX i_BeamY i_csz i_vwl code_num

% plot OD
N_lines = 51; [set_of_lines,my_colormap] = make_colormap(2,18,N_lines);
for i_BeamX = iBeamX %1:length(set_of_beam_X)
    for i_BeamY = iBeamY %1:length(set_of_beam_Y)
        for i_czs = 1:length(set_of_czs)
            for i_wvl = 1:length(set_of_wvl)
                beam_X = set_of_beam_X(i_BeamX);
                beam_Y = set_of_beam_Y(i_BeamY);
                cz = set_of_czs(i_czs);
                wvlngth = set_of_wvl(i_wvl);
                figure(1), subplot(2,2,(i_czs-1).*length(set_of_czs)+i_wvl)

                contourf(x_c,y_c,-log(squeeze(set_of_data(i_czs,i_wvl,:,:))).',N_lines,'LineColor','none','LineWidth',0.1)
                hold on, rectangle('Position',[-t_db.rx-dx/2 -t_db.ry  2*t_db.rx+eps 2*t_db.ry+eps],'Curvature',[1 1],'LineWidth',2,'LineStyle','-','EdgeColor','r'), hold off % breast area
                hold on, rectangle('Position',[-t_db.rx+dx/2 -t_db.ry  2*t_db.rx+eps 2*t_db.ry+eps],'Curvature',[1 1],'LineWidth',2,'LineStyle',':','EdgeColor','k'), hold off % breast area
                hold on, plot(t_db.beam_X,t_db.beam_Y,'h','LineWidth',1,'MarkerSize',16,'MarkerEdgeColor','k','MarkerFaceColor','w'), hold off
                set(gca,'fontsize',12), axis equal, axis tight, grid on; h = colorbar; ylabel(h,'OD (a.u.)')
                % colormap(my_colormap); clim([min(set_of_lines) max(set_of_lines)])

                title(['depth = ',num2str(cz), ' cm, \lambda = ',num2str(wvlngth),' nm; source at (',num2str(beam_X),',',num2str(beam_Y),') cm']); colormap jet
                set(gca,'xtick',[-3 0 +3]), xlim([-3 +3]); xlabel('x (cm)');
                set(gca,'ytick',[-3 0 +3]), ylim([-3 +3]); ylabel('y (cm)');
                clearvars h beam_X beam_Y cz wvlngth
            end
        end
    end
end
clearvars N_lines set_of_lines my_colormap i_BeamX i_BeamY i_czs i_wvl



% plot OD_{2cm} - OD_{1cm}
N_lines = 51; [set_of_lines,my_colormap] = make_colormap(-2,+5,N_lines);
for i_BeamX = iBeamX %1:length(set_of_beam_X)
    for i_BeamY = iBeamY %1:length(set_of_beam_Y)
        for i_wvl = 1:length(set_of_wvl)
            beam_X = set_of_beam_X(i_BeamX);
            beam_Y = set_of_beam_Y(i_BeamY);
            wvlngth = set_of_wvl(i_wvl);
            figure(2), subplot(1,2,i_wvl)

            contourf(x_c,y_c,-log(squeeze(set_of_data(2,i_wvl,:,:))).'+log(squeeze(set_of_data(1,i_wvl,:,:))).',N_lines,'LineColor','none','LineWidth',0.1)
            hold on, rectangle('Position',[-t_db.rx-dx/2 -t_db.ry  2*t_db.rx+eps 2*t_db.ry+eps],'Curvature',[1 1],'LineWidth',2,'LineStyle','-','EdgeColor','r'), hold off % breast area
            hold on, rectangle('Position',[-t_db.rx+dx/2 -t_db.ry  2*t_db.rx+eps 2*t_db.ry+eps],'Curvature',[1 1],'LineWidth',2,'LineStyle',':','EdgeColor','k'), hold off % breast area
            hold on, plot(t_db.beam_X,t_db.beam_Y,'h','LineWidth',1,'MarkerSize',16,'MarkerEdgeColor','k','MarkerFaceColor','w'), hold off
            set(gca,'fontsize',12), axis equal, axis tight, grid on; h = colorbar; ylabel(h,'OD (a.u.)')
            clim([min(set_of_lines) max(set_of_lines)])

            title(['OD_{2 cm} - OD_{1 cm}, \lambda = ',num2str(wvlngth),' nm; source at (',num2str(beam_X),',',num2str(beam_Y),') cm']); colormap jet
            set(gca,'xtick',[-3 0 +3]), xlim([-3 +3]); xlabel('x (cm)');
            set(gca,'ytick',[-3 0 +3]), ylim([-3 +3]); ylabel('y (cm)');
            clearvars h beam_X beam_Y cz wvlngth
        end
    end
end
clearvars N_lines set_of_lines my_colormap i_BeamX i_BeamY i_czs i_wvl

figure(1), saveas(gcf,[ 'OD_beamX_',num2str(set_of_beam_X(iBeamX)),'_beamY_',num2str(set_of_beam_Y(iBeamY)),'.jpg'],'jpeg')
figure(2), saveas(gcf,['dOD_beamX_',num2str(set_of_beam_X(iBeamX)),'_beamY_',num2str(set_of_beam_Y(iBeamY)),'.jpg'],'jpeg')
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

function [x_c,y_c,out] = make_2d(t_db,Lx,Ly,N_bin)
x_edges = linspace(-Lx/2,+Lx/2,N_bin+1); % # of x bins
y_edges = linspace(-Ly/2,+Ly/2,N_bin+1); % # of y bins
[N,x_edges,y_edges,ind_x,ind_y] = histcounts2(t_db.p_ot(:,1),t_db.p_ot(:,2),x_edges,y_edges);
x_c = 1/2*(x_edges(1:end-1)+x_edges(2:end-0)).';
y_c = 1/2*(y_edges(1:end-1)+y_edges(2:end-0)).';
clearvars x_edges y_edges
N_samples = length(t_db.s);

out = zeros(size(N,1),size(N,2));
for i_sample = 1:N_samples
    out(ind_x(i_sample),ind_y(i_sample)) = out(ind_x(i_sample),ind_y(i_sample)) + t_db.w(i_sample);
end
out = out./t_db.no_of_photons;
end
