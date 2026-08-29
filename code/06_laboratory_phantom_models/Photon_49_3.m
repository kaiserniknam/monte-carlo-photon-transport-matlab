function [] = Photon_49_3 ()
% Repository group: 06_laboratory_phantom_models
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% Current version: Optics Letters DPF Test (based on dataset from 39)
% Objective: Analyze how tumor-related signal influences the new DPF model

clc
close all
set(0,'DefaultFigureWindowStyle','docked')

Lx = 19.1; Ly = 19.1; Lz = 10.0; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
cx = 0.00; cy = 0.00; cz = 1.50; % Tumor depth from the skin (in cm)(0.5 to 5.0 cm)
rx = 1.00; ry = 1.00; rz = 1.00; % Tumor radius within the breast (in cm)(1.0 to 2.5 cm).
set_of_beam_X = [0.0,0.5,1.0,1.5]; % X deviation of beam in cm
set_of_beam_Y = [0.0,0.5,1.0,1.5]; % Y deviation of beam in cm
set_of_scal   = [0,1];   % set of scales
code_num = 39; L = 6; % piecewise path length

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



% plot I/d
for i_X = 1:length(set_of_beam_X)
    for i_Y = 1:length(set_of_beam_Y)
        if (i_X~=1||i_Y~=1)
            beam_X = set_of_beam_X(i_X);
            beam_Y = set_of_beam_Y(i_Y);
            detc_X = -(L./sqrt(beam_X.^2+beam_Y.^2)-1).*beam_X;
            detc_Y = -(L./sqrt(beam_X.^2+beam_Y.^2)-1).*beam_Y;
            % disp(num2str([beam_X beam_Y detc_X detc_Y norm([beam_X-detc_X beam_Y-detc_Y])]))
            % figure(1), plot([beam_X detc_X],[beam_Y detc_Y],'o-'), hold on

            eval(['x_rf = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',0),'.','_'),'.x);']);
            eval(['y_rf = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',0),'.','_'),'.y);']);
            eval(['w_rf = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',0),'.','_'),'.w);']);
            eval(['n_rf = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',0),'.','_'),'.no_of_surf_photons)+' ...
                         '(db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',0),'.','_'),'.no_of_dpth_photons);']);
            eval(['x_sg = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',1),'.','_'),'.x);']);
            eval(['y_sg = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',1),'.','_'),'.y);']);
            eval(['w_sg = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',1),'.','_'),'.w);']);
            eval(['n_sg = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',1),'.','_'),'.no_of_surf_photons)+' ...
                         '(db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',1),'.','_'),'.no_of_dpth_photons);']);

            N_bin = 150;
            od_rf = nan(1,N_bin);
            od_sg = nan(1,N_bin);
            d__rf = nan(1,N_bin);
            d__sg = nan(1,N_bin);
            d_edges = linspace(-L,0,N_bin+1);
            x_edges = linspace(detc_X,beam_X,N_bin+1);
            y_edges = linspace(detc_Y,beam_Y,N_bin+1);
            for i_d = 1:N_bin
                if     detc_X==beam_X
                    idx_rf = y_edges(i_d)<=y_rf&y_rf<=y_edges(i_d+1);
                    idx_sg = y_edges(i_d)<=y_sg&y_sg<=y_edges(i_d+1);
                elseif detc_Y==beam_Y
                    idx_rf = x_edges(i_d)<=x_rf&x_rf<=x_edges(i_d+1);
                    idx_sg = x_edges(i_d)<=x_sg&x_sg<=x_edges(i_d+1);
                else
                    idx_sg = x_edges(i_d)<=x_sg&x_sg<=x_edges(i_d+1)&y_edges(i_d)<=y_sg&y_sg<=y_edges(i_d+1);
                    idx_rf = x_edges(i_d)<=x_rf&x_rf<=x_edges(i_d+1)&y_edges(i_d)<=y_rf&y_rf<=y_edges(i_d+1);
                end
                od_rf(i_d) = sum(w_rf(idx_rf))./n_rf;
                od_sg(i_d) = sum(w_sg(idx_sg))./n_sg;
                d__rf(i_d) = mean([d_edges(i_d) d_edges(i_d+1)]);
                d__sg(i_d) = mean([d_edges(i_d) d_edges(i_d+1)]);
                clearvars idx_rf idx_sg
            end
            mua = 0.03; mus = 240;
            p = 0.493.*(mua.^(-0.507)).*(mus.^(+0.506)); % from model
            q = 3.112.*(mua.^(-1.370)).*(mus.^(-0.198)); % from model

            X = d__sg; Y = -log(od_sg);
            X = X(~isnan(Y)&~isinf(Y)); Y = Y(~isnan(Y)&~isinf(Y));
            coeffs = polyfit(-1./X, -Y./X, 1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % coeffs = polyfit(-1./X(round([1 length(X)])), -Y(round([1 length(X)]))./X(round([1 length(X)])), 1);
            % coeffs = polyfit(-1./X(round([1 length(X)/2 length(X)])), -Y(round([1 length(X)/2 length(X)]))./X(round([1 length(X)/2 length(X)])), 1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            q_esp = coeffs(1);  % slope
            p_esp = coeffs(2);  % intercept
            ft_og = p_esp - q_esp./X;
            r_squared = calc_rsquared(Y,ft_og);
            clearvars SS_res SS_tot X Y

            figure(2), subplot(4,4,(i_Y-1).*length(set_of_beam_Y)+i_X)
            plot(-d__rf,-log(od_rf)./(-d__rf),'r','LineWidth',2,'DisplayName','no tumor'), hold on
            plot(-d__sg,-log(od_sg)./(-d__sg),'b','LineWidth',1,'DisplayName','wt tumor'), hold on
            plot(-d__rf,p_esp-q_esp./d__sg,'g--','LineWidth',2,'DisplayName',['model: r^2~',num2str(r_squared,'%.2f')]), hold on
            axis([0 5 0 50])
            if i_Y==1, title(get_X_label(beam_X)); else, title("");end
            if i_Y==length(set_of_beam_Y), xlabel('d (cm)'); else, xlabel(''); end
            if i_X==1, ylabel([{[get_Y_label(beam_Y)]},{'OD/d (cm^{-1})'}]); else, ylabel(""); end
            set(gca,'xtick',0:1:5)
            set(gca,'ytick',0:25:50)
            set(gca,'fontsize',16), grid on, legend('show','Location','northeast')
            clearvars beam_X beam_Y detc_X detc_Y
            clearvars x_rf y_rf w_rf n_rf d__rf
            clearvars x_sg y_sg w_sg n_sg d__sg
            clearvars N_bin od_rf od_sg
            clearvars i_d d_edges x_edges y_edges
        end
    end
end
clearvars N_lines set_of_lines my_colormap i_czs i_scale
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
function [r_squared] = calc_rsquared (Y,Y_fit)
SS_res = sum((Y -   Y_fit).^2); % Residual sum of squares
SS_tot = sum((Y - mean(Y)).^2); % Total sum of squares
r_squared = 1 - (SS_res / SS_tot);
end
