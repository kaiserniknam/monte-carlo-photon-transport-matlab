function [] = Photon_57_1 ()
% Repository group: 07_dpf_parameter_studies
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% This version is based on version #33 but uses a smaller computational domain
% to test the method for estimating μa and μs within the simulation.
% the analysis: population stat

clc
close all
set(0,'DefaultFigureWindowStyle','docked')

% optical & geometry properties
Lx = 3.1; Ly = 3.1; Lz = 1.5; % The size of the computational domain (in cm)
set_of_mua = (0.0:.05:0.5);
set_of_mus = (000:050:500);
Pd = nan(length(set_of_mua)*length(set_of_mus),1);
Ps = nan(length(set_of_mua)*length(set_of_mus),1);
Pt = nan(length(set_of_mua)*length(set_of_mus),1);
Pa = nan(length(set_of_mua)*length(set_of_mus),1);
mav= nan(length(set_of_mua)*length(set_of_mus),1);
msv= nan(length(set_of_mua)*length(set_of_mus),1);

% arrays to save variables
dl = 0.1; % spatial resolution (cm)
N_photon = nan(length(set_of_mua),length(set_of_mus),1+4+1+1); % top + laterals + bottom

for i_a = 1:1:length(set_of_mua)
    for i_s = 1:1:length(set_of_mus)
        % read dbase
        mua = set_of_mua(i_a);
        mus = set_of_mus(i_s);
        the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_57_mua_',sprintf('%.2f',mua),'_mus_',sprintf('%.2f',mus),'.mat'];
        t_db = load(the_filename); clearvars the_filename
        disp(['mua = ',num2str(mua,'%.2f'),', mus = ',num2str(mus,'%.0f')])

        % removing noisy points
        u = round(-log(t_db.w)./t_db.s,4); u_unique = unique(u); freq = nan(size(u_unique));
        for i_f = 1:length(u_unique), freq(i_f) = sum(u==u_unique(i_f)); end
        [~,i_fx] = max(freq);
        t_db.x_in = t_db.x_in(u==u_unique(i_fx));
        t_db.y_in = t_db.y_in(u==u_unique(i_fx));
        t_db.z_in = t_db.z_in(u==u_unique(i_fx));
        t_db.x_ot = t_db.x_ot(u==u_unique(i_fx));
        t_db.y_ot = t_db.y_ot(u==u_unique(i_fx));
        t_db.z_ot = t_db.z_ot(u==u_unique(i_fx));
        t_db.s = t_db.s(u==u_unique(i_fx));
        t_db.w = t_db.w(u==u_unique(i_fx));
        disp(num2str([length(u) sum(u==u_unique(i_fx)) sum(u==u_unique(i_fx))/length(u)*100]))
        t_db.no_of_photons = t_db.no_of_photons-sum(u~=u_unique(i_fx));
        clearvars u u_unique freq i_f i_fx

        the_limit = [];
        % Count number of photons on each side
        idx = -Lx/2<t_db.x_ot & t_db.x_ot<+Lx/2 & -Ly/2<t_db.y_ot & t_db.y_ot<+Ly/2 & t_db.z_ot<=00; % -Z-axis
        N_photon(i_a,i_s,1) = sum(idx)./t_db.no_of_photons;
        x_edges = linspace(-Lx/2,+Lx/2,round(Lx./dl)+1);
        y_edges = linspace(-Ly/2,+Ly/2,round(Ly./dl)+1);
        [A,x_c,y_c] = make2DMap(x_edges,y_edges,t_db.x_ot(idx),t_db.y_ot(idx),t_db.w(idx),t_db.no_of_photons);
        subplot(2,2,2), contourf(x_c,y_c,-log(A.')), shading flat
        xlabel('x (cm)'), set(gca,'xtick',[-Lx/2 0 +Lx/2])
        ylabel('y (cm)'), set(gca,'ytick',[-Ly/2 0 +Ly/2])
        title(['\mu_a = ',num2str(mua),' cm^{-1}, \mu_s = ',num2str(mus),' cm^{-1}'],'top surface'), set(gca,'fontsize',14), axis equal tight, colormap jet; colorbar;
        the_limit = [ the_limit ; [min(A(:)) max(A(:))]]; clim([4 14])
        clearvars A x_c y_c idx x_edges y_edges

        idx = -Lx/2<t_db.x_ot & t_db.x_ot<+Lx/2 & -Ly/2<t_db.y_ot & t_db.y_ot<+Ly/2 & Lz<=t_db.z_ot; % +Z-axis
        N_photon(i_a,i_s,6) = sum(idx)./t_db.no_of_photons;
        x_edges = linspace(-Lx/2,+Lx/2,round(Lx./dl)+1);
        y_edges = linspace(-Ly/2,+Ly/2,round(Ly./dl)+1);
        [A,x_c,y_c] = make2DMap(x_edges,y_edges,t_db.x_ot(idx),t_db.y_ot(idx),t_db.w(idx),t_db.no_of_photons);
        subplot(2,2,4), contourf(x_c,y_c,-log(A.')), shading flat
        xlabel('x (cm)'), set(gca,'xtick',[-Lx/2 0 +Lx/2])
        ylabel('y (cm)'), set(gca,'ytick',[-Ly/2 0 +Ly/2])
        title('bottom surface'), set(gca,'fontsize',14), axis equal tight, colormap jet; colorbar;
        the_limit = [ the_limit ; [min(A(:)) max(A(:))]]; clim([4 14])
        clearvars A x_c y_c idx x_edges y_edges



        idx = 0<t_db.z_ot & t_db.z_ot<Lz & -Ly/2<t_db.y_ot & t_db.y_ot<+Ly/2 & +Lx/2<=t_db.x_ot; % +X-axis
        N_photon(i_a,i_s,2) = sum(idx)./t_db.no_of_photons;
        x_edges = linspace(-Ly/2,+Ly/2,round(Ly./dl)+1);
        y_edges = linspace(0,    +Lz,  round(Lz./dl)+1);
        [A,x_c,y_c] = make2DMap(x_edges,y_edges,t_db.y_ot(idx),t_db.z_ot(idx),t_db.w(idx),t_db.no_of_photons);
        subplot(4,2,1), contourf(x_c,-y_c,-log(A.')), shading flat
        set(gca,'xtick',[-Ly/2 0 +Ly/2]), xlim([-Ly/2 +Ly/2])
        ylabel('z (cm)'), set(gca,'ytick',[-Lz         0]), ylim([-Lz    0   ])
        title('+X surface'), set(gca,'fontsize',12), colormap jet; colorbar;
        the_limit = [ the_limit ; [min(A(:)) max(A(:))]]; clim([4 14])
        clearvars A x_c y_c idx x_edges y_edges

        idx = 0<t_db.z_ot & t_db.z_ot<Lz & -Ly/2<t_db.y_ot & t_db.y_ot<+Ly/2 & t_db.x_ot<=-Lx/2; % -X-axis
        N_photon(i_a,i_s,3) = sum(idx)./t_db.no_of_photons;
        x_edges = linspace(-Ly/2,+Ly/2,round(Ly./dl)+1);
        y_edges = linspace(0,    +Lz,  round(Lz./dl)+1);
        [A,x_c,y_c] = make2DMap(x_edges,y_edges,t_db.y_ot(idx),t_db.z_ot(idx),t_db.w(idx),t_db.no_of_photons);
        subplot(4,2,3), contourf(x_c,-y_c,-log(A.')), shading flat
        xlabel('y (cm)'), set(gca,'xtick',[-Ly/2 0 +Ly/2]), xlim([-Ly/2 +Ly/2])
        ylabel('z (cm)'), set(gca,'ytick',[-Lz         0]), ylim([-Lz    0   ])
        title('-X surface'), set(gca,'fontsize',12), colormap jet; colorbar;
        the_limit = [ the_limit ; [min(A(:)) max(A(:))]]; clim([4 14])
        clearvars A x_c y_c idx x_edges y_edges


        idx = 0<t_db.z_ot & t_db.z_ot<Lz & -Lx/2<t_db.x_ot & t_db.x_ot<+Lx/2 & +Ly/2<=t_db.y_ot; % +Y-axis
        N_photon(i_a,i_s,4) = sum(idx)./t_db.no_of_photons;
        x_edges = linspace(-Lx/2,+Lx/2,round(Lx./dl)+1);
        y_edges = linspace(0,    +Lz,  round(Lz./dl)+1);
        [A,x_c,y_c] = make2DMap(x_edges,y_edges,t_db.x_ot(idx),t_db.z_ot(idx),t_db.w(idx),t_db.no_of_photons);
        subplot(4,2,5), contourf(x_c,-y_c,-log(A.')), shading flat
        set(gca,'xtick',[-Lx/2 0 +Lx/2]), xlim([-Lx/2 +Lx/2])
        ylabel('z (cm)'), set(gca,'ytick',[-Lz         0]), ylim([-Lz    0   ])
        title('+Y surface'), set(gca,'fontsize',12), colormap jet; colorbar;
        the_limit = [ the_limit ; [min(A(:)) max(A(:))]]; clim([4 14])
        clearvars A x_c y_c idx x_edges y_edges

        idx = 0<t_db.z_ot & t_db.z_ot<Lz & -Lx/2<t_db.x_ot & t_db.x_ot<+Lx/2 & t_db.y_ot<=-Ly/2; % -Y-axis
        N_photon(i_a,i_s,5) = sum(idx)./t_db.no_of_photons;
        x_edges = linspace(-Lx/2,+Lx/2,round(Lx./dl)+1);
        y_edges = linspace(0,    +Lz,  round(Lz./dl)+1);
        [A,x_c,y_c] = make2DMap(x_edges,y_edges,t_db.x_ot(idx),t_db.z_ot(idx),t_db.w(idx),t_db.no_of_photons);
        subplot(4,2,7), contourf(x_c,-y_c,-log(A.')), shading flat
        xlabel('x (cm)'), set(gca,'xtick',[-Lx/2 0 +Lx/2]), xlim([-Lx/2 +Lx/2])
        ylabel('z (cm)'), set(gca,'ytick',[-Lz         0]), ylim([-Lz    0   ])
        title('-Y surface'), set(gca,'fontsize',12), colormap jet; colorbar;
        the_limit = [ the_limit ; [min(A(:)) max(A(:))]]; clim([4 14])
        clearvars A x_c y_c idx x_edges y_edges

        idx = -Lx/2<t_db.x_ot & t_db.x_ot<+Lx/2 & -Ly/2<t_db.y_ot & t_db.y_ot<+Ly/2 & 0<t_db.z_ot & t_db.z_ot<Lz;
        N_photon(i_a,i_s,7) = sum(idx)./t_db.no_of_photons; clearvars idx % inside cube

        % mua, mus -> column vectors [N×1]
        % Pd, Ps, Pt, Pa  -> column vectors [N×1]
        idx = sub2ind([length(set_of_mua),length(set_of_mus)],i_a,i_s);
        Pd(idx) = N_photon(i_a,i_s,1);
        Ps(idx) = N_photon(i_a,i_s,2) + N_photon(i_a,i_s,3) + N_photon(i_a,i_s,4) + N_photon(i_a,i_s,5);
        Pt(idx) = N_photon(i_a,i_s,6);
        Pa(idx) = N_photon(i_a,i_s,7);
        mav(idx)= set_of_mua(i_a);
        msv(idx)= set_of_mus(i_s);



        clearvars t_db my_colormap mua mus idx the_limit
    end
end
clearvars i_a i_s

p_diffuse = squeeze(N_photon(:,:,1)).';
p_sidedph = squeeze(sum(N_photon(:,:,2:5),3)).';
p_trnsmit = squeeze(N_photon(:,:,6)).';
p_absorbd = squeeze(N_photon(:,:,7)).';

figure(201)
subplot(2,2,1), contourf(set_of_mua,set_of_mus,p_diffuse), shading flat
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
title('percent of diffused photons'),
set(gca,'fontsize',16), axis square, axis tight,
colormap jet; colorbar; % clim([0 1]); set(gca,'colorscale','log');
subplot(2,2,2), contourf(set_of_mua,set_of_mus,p_sidedph), shading flat
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
title('Percentage of photons transmitted to the sides')
set(gca,'fontsize',16), axis square, axis tight,
colormap jet; colorbar; % clim([0 1]); set(gca,'colorscale','log');
subplot(2,2,3), contourf(set_of_mua,set_of_mus,p_trnsmit), shading flat
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
title('Percentage of photons transmitted to the bottom side')
set(gca,'fontsize',16), axis square, axis tight,
colormap jet; colorbar; % clim([0 1]); set(gca,'colorscale','log');
subplot(2,2,4), contourf(set_of_mua,set_of_mus,p_absorbd), shading flat
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
title('Percentage of absorbed photons')
set(gca,'fontsize',16), axis square, axis tight,
colormap jet; colorbar; % clim([0 1]); set(gca,'colorscale','log');



% ============================================================
% Inverse model: estimate μa and μs from photon percentages
% ============================================================
% Predictor and response matrices
X = [Pd, Ps, Pt, Pa];
Y = [mav, msv];
% 1. Fit linear models
mdl_mua_lin = fitlm(X, Y(:,1));
mdl_mus_lin = fitlm(X, Y(:,2));
Ypred_lin_mua = predict(mdl_mua_lin, X);
Ypred_lin_mus = predict(mdl_mus_lin, X);
% 2. Fit nonlinear models (train + evaluate on same data)
mdl_mua_nl = fitrnet(X, Y(:,1), "LayerSizes", [10,10], "Standardize", true);
mdl_mus_nl = fitrnet(X, Y(:,2), "LayerSizes", [10,10], "Standardize", true);
Ypred_nl_mua = predict(mdl_mua_nl, X);
Ypred_nl_mus = predict(mdl_mus_nl, X);
% 3. Compute goodness-of-fit metrics
goodness = @(y, ypred) struct( ...
    'R2', 1 - var(y - ypred)/var(y), ...
    'RMSE', sqrt(mean((y - ypred).^2)), ...
    'r', corr(y, ypred) );
g_lin_mua = goodness(Y(:,1), Ypred_lin_mua);
g_lin_mus = goodness(Y(:,2), Ypred_lin_mus);
g_nl_mua  = goodness(Y(:,1), Ypred_nl_mua);
g_nl_mus  = goodness(Y(:,2), Ypred_nl_mus);
% 4. Display results
clc
fprintf('\nGoodness-of-fit (Training Data Only)\n');
fprintf('------------------------------------------------------------\n');
fprintf('Model     Parameter     R^2        RMSE        Corr(r)\n');
fprintf('Linear    μa        %.4f     %.4f     %.4f\n', g_lin_mua.R2, g_lin_mua.RMSE, g_lin_mua.r);
fprintf('Linear    μs        %.4f     %.4f     %.4f\n', g_lin_mus.R2, g_lin_mus.RMSE, g_lin_mus.r);
fprintf('Nonlinear μa        %.4f     %.4f     %.4f\n', g_nl_mua.R2,  g_nl_mua.RMSE,  g_nl_mua.r);
fprintf('Nonlinear μs        %.4f     %.4f     %.4f\n', g_nl_mus.R2,  g_nl_mus.RMSE,  g_nl_mus.r);
fprintf('------------------------------------------------------------\n');
% 5. Visualization #1
figure(202)
subplot(1,2,1)
plot(Y(:,1), Ypred_lin_mua,'o','MarkerSize',12,'MarkerFaceColor','b','MarkerEdgeColor','k'); hold on
plot(Y(:,1), Ypred_nl_mua, 'o','MarkerSize',12,'MarkerFaceColor','r','MarkerEdgeColor','k'); refline(1,0)
xlabel('True \mu_a (cm^{-1})'); ylabel('Predicted \mu_a (cm^{-1})');
legend('Linear','Nonlinear','1:1 Line','Location','best')
title(sprintf('\\mu_a: Linear R^2=%.3f | Nonlinear R^2=%.3f', g_lin_mua.R2, g_nl_mua.R2))
grid on; axis equal tight
set(gca,'fontsize',18)
subplot(1,2,2)
plot(Y(:,2), Ypred_lin_mus,'o','MarkerSize',12,'MarkerFaceColor','b','MarkerEdgeColor','k'); hold on
plot(Y(:,2), Ypred_nl_mus, 'o','MarkerSize',12,'MarkerFaceColor','r','MarkerEdgeColor','k'); refline(1,0)
xlabel('True \mu_s (cm^{-1})'); ylabel('Predicted \mu_s (cm^{-1})');
legend('Linear','Nonlinear','1:1 Line','Location','best')
title(sprintf('\\mu_s: Linear R^2=%.3f | Nonlinear R^2=%.3f', g_lin_mus.R2, g_nl_mus.R2))
grid on; axis equal tight
set(gca,'fontsize',18)


% 5. Visualization #2
figure(203)
subplot(1,2,1)
mesh(set_of_mua,set_of_mus,ConvToVector(Ypred_lin_mua,set_of_mua,set_of_mus).','FaceColor','none','EdgeColor','b','DisplayName','linear'), hold on
mesh(set_of_mua,set_of_mus,ConvToVector(Ypred_nl_mua, set_of_mua,set_of_mus).','FaceColor','none','EdgeColor','r','DisplayName','nonlinear'), hold on
xlabel('\mu_a (cm^{-1})'); ylabel('\mu_s (cm^{-1})'); zlabel('Predicted \mu_a (cm^{-1})');
legend('show','Location','best')
title(sprintf('\\mu_a: Linear R^2=%.3f | Nonlinear R^2=%.3f', g_lin_mua.R2, g_nl_mua.R2))
grid on; axis tight
set(gca,'fontsize',18)
subplot(1,2,2)
mesh(set_of_mua,set_of_mus,ConvToVector(Ypred_lin_mus,set_of_mua,set_of_mus).','FaceColor','none','EdgeColor','b','DisplayName','linear'), hold on
mesh(set_of_mua,set_of_mus,ConvToVector(Ypred_nl_mus, set_of_mua,set_of_mus).','FaceColor','none','EdgeColor','r','DisplayName','nonlinear'), hold on
xlabel('\mu_a (cm^{-1})'); ylabel('\mu_s (cm^{-1})'); zlabel('Predicted \mu_s (cm^{-1})');
legend('show','Location','best')
title(sprintf('\\mu_s: Linear R^2=%.3f | Nonlinear R^2=%.3f', g_lin_mus.R2, g_nl_mus.R2))
grid on; axis tight
set(gca,'fontsize',18)
end

function [Y_2D] = ConvToVector (Y, set_of_mua, set_of_mus)
Y_2D = nan(length(set_of_mua),length(set_of_mus));
for i_a = 1:1:length(set_of_mua)
    for i_s = 1:1:length(set_of_mus)
        idx = sub2ind([length(set_of_mua),length(set_of_mus)],i_a,i_s);
        Y_2D(i_a,i_s) = Y(idx); clearvars idx
    end
end
end
function [A,x_c,y_c] = make2DMap(x_edges,y_edges,x,y,z,no_of_photons)
% Bin indices
[~, ~, ~, ind_x, ind_y] = histcounts2(x, y, x_edges, y_edges);
% Remove invalid bins
mask  = (ind_x > 0) & (ind_y > 0);
ind_x = ind_x(mask);
ind_y = ind_y(mask);
z     = z(mask); % values to accumulate
Nx = length(x_edges)-1;
Ny = length(y_edges)-1;
% Accumulate z-values into bins
A = accumarray([ind_x(:) ind_y(:)], z(:)./no_of_photons, [Nx Ny], @sum, nan);
% Bin centers
x_c = 0.5.*(x_edges(1:end-1) + x_edges(2:end));
y_c = 0.5.*(y_edges(1:end-1) + y_edges(2:end));
end
