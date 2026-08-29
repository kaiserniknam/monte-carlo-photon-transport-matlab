function [] = Photon_49_6 ()
% Repository group: 06_laboratory_phantom_models
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% Current version: Optics Letters DPF Test (based on dataset from 39)
% Objective: Determine the minimum number of measurement points required to experimentally estimate the empirical Differential Pathlength Factor (DPF)
% the same as Photon_49_5, but only 3-7 data points are selected to estimate parameters p and q; and also I0 is not used.

clc
close all
set(0,'DefaultFigureWindowStyle','docked')

set_of_beam_X = [0.0,0.5,1.0,1.5]; % X deviation of beam in cm
set_of_beam_Y = [0.0,0.5,1.0,1.5]; % Y deviation of beam in cm
set_of_scal   = [0,1];   % set of scales
code_num = 39; L = 6; % piecewise path length
N_bin = 250; N_rnd = 100;
set_of_N_points = 3:1:N_bin;

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



% calc p, and q of I/d
for i_X = 1:length(set_of_beam_X)
    for i_Y = 1:length(set_of_beam_Y)
        if (i_X~=1||i_Y~=1)
            beam_X = set_of_beam_X(i_X);
            beam_Y = set_of_beam_Y(i_Y);
            detc_X = -(L./sqrt(beam_X.^2+beam_Y.^2)-1).*beam_X;
            detc_Y = -(L./sqrt(beam_X.^2+beam_Y.^2)-1).*beam_Y;
            set_of_pq__sets = nan(length(set_of_N_points)+1,2+1);
            set_of_p___sets = nan(length(set_of_N_points)+1,N_rnd);

            eval(['x_sg = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',1),'.','_'),'.x);']);
            eval(['y_sg = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',1),'.','_'),'.y);']);
            eval(['w_sg = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',1),'.','_'),'.w);']);
            eval(['n_sg = (db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',1),'.','_'),'.no_of_surf_photons)+' ...
                         '(db_',replace(sprintf('%.1f',beam_X),'.','_'),'_',replace(sprintf('%.1f',beam_Y),'.','_'),'_',replace(sprintf('%.1f',1),'.','_'),'.no_of_dpth_photons);']);

            od_sg = nan(1,N_bin);
            d__sg = nan(1,N_bin);
            d_edges = linspace(-L,0,N_bin+1);
            x_edges = linspace(detc_X,beam_X,N_bin+1);
            y_edges = linspace(detc_Y,beam_Y,N_bin+1);
            for i_d = 1:N_bin
                if     detc_X==beam_X
                    idx_sg = y_edges(i_d)<=y_sg&y_sg<=y_edges(i_d+1);
                elseif detc_Y==beam_Y
                    idx_sg = x_edges(i_d)<=x_sg&x_sg<=x_edges(i_d+1);
                else
                    idx_sg = x_edges(i_d)<=x_sg&x_sg<=x_edges(i_d+1)&y_edges(i_d)<=y_sg&y_sg<=y_edges(i_d+1);
                end
                od_sg(i_d) = sum(w_sg(idx_sg))/n_sg;
                d__sg(i_d) = mean([d_edges(i_d) d_edges(i_d+1)]);
                clearvars idx_rf idx_sg
            end
            clearvars i_d

            X = d__sg; Y = -log(od_sg);
            X = X(~isnan(Y)&~isinf(Y)); Y = Y(~isnan(Y)&~isinf(Y));
            % coeffs = polyfit(-1./X, -Y./X, 1);
            coeffs = polyfit(-X, Y, 1); coeffs = coeffs(end:-1:1);
            q = coeffs(1);  % slope
            p = coeffs(2);  % intercept
            ft_og = p - q./X;
            r_squared = calc_rsquared(-Y./X,ft_og);
            set_of_pq__sets(end,1) = p;
            set_of_pq__sets(end,2) = q;
            set_of_pq__sets(end,3) = r_squared;
            clearvars p q ft_og r_squared coeffs

            for i_p = 1:length(set_of_N_points)
                idx = round(linspace(1,length(X),set_of_N_points(i_p)));
                for i_rnd = 1:N_rnd
                    Xt = X(idx); Yt = Y(idx);
                    irnd = randperm(length(Xt),1);
                    Xt = Xt-Xt(irnd); Xt(irnd) = [];
                    Yt = Yt-Yt(irnd); Yt(irnd) = [];
                    coeffs = polyfit(Xt, Yt, 1);
                    set_of_p___sets(i_p,i_rnd) = coeffs(2);
                    clearvars irnd coeffs
                end
                q = mean(set_of_p___sets(i_p,:));
                p = nan; % intercept
                ft_og = p - q./X;
                r_squared = calc_rsquared(-Y./X,ft_og);
                set_of_pq__sets(i_p,:) = [p q r_squared];
                disp(num2str([set_of_N_points(i_p) p set_of_pq__sets(end,1)]))
                clearvars idx i_rnd p q ft_og r_squared Xt Yt
            end
            clearvars i_p

            figure(3),
            X  = set_of_N_points;
            subplot(2,2,1), plot(X,set_of_pq__sets(end,2).*ones(size(X)),'LineWidth',2,'Color','b','DisplayName','true p'); hold on
            subplot(2,2,1), plot(X,set_of_pq__sets(1:end-1,2),'LineWidth',2,'Color','r','DisplayName','estimated p'); hold on
            xlabel('# of samples'); ylabel('p'); title('the evolution of p vs. # of samples');
            set(gca,'fontsize',18), hold off, axis tight, grid on
            clearvars X Y sd
            %
            % figure(3),
            % X  = set_of_N_points;
            % Y  = set_of_pq__sets(:,2);
            % subplot(2,2,2), plot(X,Y(end).*ones(size(X)),'LineWidth',2,'Color','b','DisplayName','true q'); hold on
            % subplot(2,2,2), plot(X,Y(1:end-1),'LineWidth',2,'Color','r','DisplayName','estimated q'); hold on
            % xlabel('# of samples'); ylabel('p'); title('the evolution of q vs. # of samples');
            % set(gca,'fontsize',18), hold off, axis tight, grid on
            % clearvars X Y sd
            %
            % figure(3),
            % X  = set_of_N_points;
            % Y  = set_of_pq__sets(:,3);
            % subplot(2,1,2), plot(X,Y(end).*ones(size(X)),'LineWidth',2,'Color','b','DisplayName','true q'); hold on
            % subplot(2,1,2), plot(X,Y(1:end-1),'LineWidth',2,'Color','r','DisplayName','estimated q'); hold on
            % xlabel('# of samples'); ylabel('r^2'); title('the evolution of r^2 vs. # of samples');
            % set(gca,'fontsize',18), hold off, axis tight, grid on
            % clearvars X Y sd



            clearvars beam_X beam_Y detc_X detc_Y
            clearvars x_rf y_rf w_rf n_rf d__rf
            clearvars x_sg y_sg w_sg n_sg d__sg
            clearvars ft_og set_of_p___sets
            clearvars i_d d_edges x_edges y_edges
        end
    end
end
clearvars N_lines set_of_lines my_colormap i_czs i_scale
end

function [r_squared] = calc_rsquared (Y,Y_fit)
SS_res = sum((Y -   Y_fit).^2); % Residual sum of squares
SS_tot = sum((Y - mean(Y)).^2); % Total sum of squares
r_squared = 1 - (SS_res / SS_tot);
end
