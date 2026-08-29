function [] = Photon_69_1 ()
% Repository group: 09_concentration_recovery
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% This Version: TiO2 experimental phantom (India ink tumor inclusion in gelatin)
% Configuration:
%   Background: SAP hydrogel + TiO2 scattering matrix
%   Inclusion: Localized absorbing tumor (India ink embedded in gelatin)
% Purpose:
%   Validation of DPF-based absorption estimation under experimentally realistic conditions.
%
%   Based on Version 68, adapted for the TiO2 experimental dataset.
%
%   New tumor depth = 1.00 cm / sds = 2 cm (previously 1.2 cm / sds = 3 cm).
% analysis data


SDS = 2; depth = 1.0;     % depth from top surface to circle center
% ================= Plot =================
clc
close all
set(0,'DefaultFigureWindowStyle','docked')

W = 8;           % rectangle width
H = 5.3;         % rectangle height
d_circle = 0.75; % circle diameter
r_circle = d_circle / 2;
x_rect = 0;
y_rect = 0;
xc = W/2;               % centered horizontally
yc = H - depth;         % circle center

figure(1), subplot(1,2,1)
hold on;
% Rectangle (phantom)
rectangle('Position', [x_rect, y_rect, W, H], ...
          'EdgeColor', 'k', 'LineWidth', 2);
% Circle (inclusion)
theta = linspace(0, 2*pi, 300);
x_circle = xc + r_circle*cos(theta);
y_circle = yc + r_circle*sin(theta);
plot(x_circle, y_circle, 'r', 'LineWidth', 2);
% Center point
plot(xc, yc, 'bo', 'MarkerFaceColor', 'b');
% ================= Source =================
xs = xc;     % directly above inclusion
ys = H;      % top surface
% draw small rectangle for source
src_size = 0.15;
rectangle('Position', [xs - src_size/2, ys, src_size, src_size], ...
          'FaceColor', 'r', 'EdgeColor', 'k');
% ================= Detector =================
xd = xs + SDS;     % to the right
yd = H;
% draw small rectangle for detector
det_size = 0.15;
rectangle('Position', [xd - det_size/2, yd, det_size, det_size], ...
          'FaceColor', 'b', 'EdgeColor', 'k');
% ================= Optional: line showing distance =================
plot([xs xd], [ys yd], '--k');
% ================= Final settings =================
title('Phantom with Inclusion, Source, and Detector');
xlim([0, max(W, xd) + 0.5]);
ylim([0, H + 1]);
axis equal tight
grid on, hold off
xlabel(''), xticks(0:10), set(gca,'xticklabel',-4:5)
ylabel(''), yticks(0:10)
set(gca,'fontsize',18)
clearvars -except SDS depth
% ============= End Plot =================

L = 100;
set_of_concentration = [nan, 0.00:0.01:0.5]; % set of concentrations
set_of_mua = nan(size(set_of_concentration)); % set of mu_a's
set_of_ODs = nan(length(set_of_concentration),L,2);
set_of_ODs_d = nan(length(set_of_concentration),L);
set_of_ODs_d_mua = nan(length(set_of_concentration),L);
clearvars L
for i_concentration = 1:length(set_of_concentration)
    the_concentration = set_of_concentration(i_concentration);
    the_color = [(i_concentration-1)./(length(set_of_concentration)-1) 0 1-(i_concentration-1)./(length(set_of_concentration)-1)];
    if isnan(the_concentration)
        set_of_mua(i_concentration) = nan;
    else
        if the_concentration <= 0.06
            set_of_mua(i_concentration) = max(0.0275,-0.55+67.15.*the_concentration);
        else
            set_of_mua(i_concentration) = 3.63;
        end
    end
    the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_69_concentration_',num2str(the_concentration),'.mat' ];
    [y_bind,x_bind] = get_OD (the_filename,0.17);
    G = -log(y_bind(1));
    set_of_ODs(i_concentration,1:length(x_bind),1) = x_bind;
    set_of_ODs(i_concentration,1:length(y_bind),2) = y_bind;
    set_of_ODs_d(i_concentration,1:length(y_bind)) = (-log(y_bind)-G)./x_bind;
    set_of_ODs_d_mua(i_concentration,1:length(y_bind)) = (-log(y_bind)-G)./x_bind./set_of_mua(i_concentration);

    figure(1), subplot(1,2,2), plot(x_bind,-log(y_bind),'LineStyle','-','LineWidth',2,'Color',the_color,'DisplayName',num2str(the_concentration)), hold on
    xlabel('d (cm)'), ylabel('OD (a.u.)'), set(gca,'fontsize',16)
    title('Effect of Inclusion Concentration on Optical Density'), grid off, axis tight square
    clearvars the_concentration the_color the_filename y_bind x_bind G
end
clearvars i_concentration

dims = axis;
set_of_ODs_at_SDS = nan(length(set_of_concentration),1);
for i_concentration = 1:length(set_of_concentration)
    [~,idx] = min(abs(set_of_ODs(i_concentration,:,1)-SDS));
    figure(1), subplot(1,2,2), plot([set_of_ODs(i_concentration,idx,1) set_of_ODs(i_concentration,idx,1)],[dims(3) dims(4)],'LineStyle','-.','LineWidth',2,'Color','k','HandleVisibility','off'),
    set_of_ODs_at_SDS(i_concentration) = set_of_ODs(i_concentration,idx,2);
    clearvars idx
end
figure(2), plot(set_of_concentration,-log(set_of_ODs_at_SDS),'LineStyle','-','LineWidth',2,'Color','k'),
xlabel('concentration (%v/v)'), ylabel('OD (a.u.)'), set(gca,'fontsize',16)
title(['Effect of Inclusion Concentration on Optical Density at SDS = ',num2str(SDS),' cm']), grid on
clearvars dims set_of_ODs_at_3cm i_concentration



set_of_params = nan(length(set_of_concentration),5);
for i_concentration = 2:length(set_of_concentration)
    X = squeeze(set_of_ODs(i_concentration,:,1)).';
    Y = squeeze(set_of_ODs_d_mua(i_concentration,:)).';
    idx = ~isnan(X)&~isnan(Y); X = X(idx); Y = Y(idx); clearvars idx
    the_color = [(i_concentration-1)./(length(set_of_concentration)-1) 0 1-(i_concentration-1)./(length(set_of_concentration)-1)];
    figure(3), subplot(1,2,1), plot(X,Y,'LineStyle','-','LineWidth',2,'Color',the_color,'DisplayName',num2str(set_of_concentration(i_concentration))), hold on
    invX = 1 ./ X;                                     % predictor 1
    logX_over_X = log(X) ./ X;                         % predictor 2
    T = table(Y, invX, logX_over_X, 'VariableNames', {'Y','invX','logX_over_X'});
    mdl = fitlm(T, 'Y ~ 1 + invX + logX_over_X');      % A, B, C
    clearvars T invX logX_over_X
    set_of_params(i_concentration,1) = mdl.Coefficients.Estimate(1);
    set_of_params(i_concentration,2) = mdl.Coefficients.Estimate(2);
    set_of_params(i_concentration,3) = mdl.Coefficients.Estimate(3);
    set_of_params(i_concentration,4) = mdl.Coefficients.pValue(2);
    Yhat = set_of_params(i_concentration,1) + set_of_params(i_concentration,2)./X + set_of_params(i_concentration,3).*log(X)./X;
    figure(3), subplot(1,2,1), plot(X,Yhat,'LineStyle','-.','LineWidth',2,'Color',the_color,'HandleVisibility','off'), hold on
    set_of_params(i_concentration,5) = rsquared(Yhat,Y);
    clearvars Yhat Y X mdl the_color
end
xlabel('d (cm)'), ylabel('DPF (a.u.)'), set(gca,'fontsize',16)
title('Effect of Inclusion Concentration on DPF'), grid off, axis tight square
figure(3), subplot(2,2,2)
plot(set_of_concentration,set_of_params(:,4),'LineStyle','-','LineWidth',2,'Color','k','DisplayName','')
xlabel('d (cm)'), ylabel('p-value (a.u)'), set(gca,'fontsize',16)
title({'P-value of Inverse Distance Model'},{' '}), grid on, axis tight square
figure(3), subplot(2,2,4)
plot(set_of_concentration,set_of_params(:,5),'LineStyle','-','LineWidth',2,'Color','k','DisplayName','')
xlabel('d (cm)'), ylabel('r-squared (a.u)'), set(gca,'fontsize',16)
title('R-squared of Inverse Distance Model'), grid on, axis tight square
clearvars i_concentration

[params_A,err_A] = fitPlaneModel(log(set_of_mua(2:end)), log(set_of_params(2:end,1)));
A_est = exp(params_A(1)).*((set_of_mua).^params_A(2));
[params_B,err_B] = fitPlaneModel(log(set_of_mua(2:end)), log(set_of_params(2:end,2)));
B_est = exp(params_B(1)).*((set_of_mua).^params_B(2));
[params_C,err_C] = fitPlaneModel(log(set_of_mua(2:end)), log(set_of_params(2:end,3)));
C_est = exp(params_C(1)).*((set_of_mua).^params_C(2));

figure(4), subplot(1,3,1)
plot(set_of_mua,set_of_params(:,1),'b','LineWidth',2,'LineStyle','-', 'DisplayName','A'), hold on
plot(set_of_mua,A_est,             'b','LineWidth',2,'LineStyle','-.','DisplayName',['estim. A, r^2 ~ ',num2str(err_A)]), hold on
plot(set_of_mua,set_of_params(:,2),'r','LineWidth',2,'LineStyle','-', 'DisplayName','B'), hold on
plot(set_of_mua,B_est,             'r','LineWidth',2,'LineStyle','-.','DisplayName',['estim. B, r^2 ~ ',num2str(err_B)]), hold on
plot(set_of_mua,set_of_params(:,3),'g','LineWidth',2,'LineStyle','-', 'DisplayName','C'), hold on
plot(set_of_mua,C_est,             'g','LineWidth',2,'LineStyle','-.','DisplayName',['estim. C, r^2 ~ ',num2str(err_C)]), hold on
xlabel('\mu_a (cm^{-1})');
ylabel('Inverse-distance parameters'); set(gca,'YScale','log')
set(gca,'FontSize',16);
title('evolution of inv. dist. parameters');
grid on; axis square; legend('show','Location','northeast')
disp(num2str(['C_A = ',num2str(exp(params_A(1)),'%.2f'),', p_A = ',num2str(params_A(2),'%.2f'),', @ r2 =  ',num2str(err_A)])), clearvars err params
disp(num2str(['C_B = ',num2str(exp(params_B(1)),'%.2f'),', p_B = ',num2str(params_B(2),'%.2f'),', @ r2 =  ',num2str(err_B)])), clearvars err params
disp(num2str(['C_C = ',num2str(exp(params_C(1)),'%.2f'),', p_C = ',num2str(params_C(2),'%.2f'),', @ r2 =  ',num2str(err_C)])), clearvars err params
clearvars A_est B_est C_est
clearvars err_A err_B err_C

figure(4), subplot(1,3,2)
for i_concentration = 2:length(set_of_concentration)
    the_color = [(i_concentration-1)./(length(set_of_concentration)-1) 0 1-(i_concentration-1)./(length(set_of_concentration)-1)];
    d  = squeeze(set_of_ODs(i_concentration,:,1));
    OD = -log(squeeze(set_of_ODs(i_concentration,:,2)));
    G = OD(1);
    the_dpf = set_of_ODs_d_mua(i_concentration,:);
    plot(d,(OD-G)./the_dpf./d,'LineStyle','-','LineWidth',1,'Color',the_color,'DisplayName',num2str(set_of_concentration(i_concentration))), hold on
    plot(d,set_of_mua(i_concentration).*ones(size(d)),'LineStyle','-.','LineWidth',2,'Color',the_color,'HandleVisibility','off'), hold on
    xlabel('d (cm)'), ylabel('\mu_a (cm^{-1})'), set(gca,'fontsize',16)
    title('Estimation of \mu_a using true model (assisted)'), grid off, axis tight square
    clearvars the_color the_dpf d OD G
end
clearvars i_concentration

figure(4), subplot(1,3,3)
for i_concentration = 2:length(set_of_concentration)
    the_color = [(i_concentration-1)./(length(set_of_concentration)-1) 0 1-(i_concentration-1)./(length(set_of_concentration)-1)];
    d  = squeeze(set_of_ODs(i_concentration,:,1));
    OD = -log(squeeze(set_of_ODs(i_concentration,:,2)));
    G = OD(1);
    A_est = exp(params_A(1)).*(set_of_mua(i_concentration)).^params_A(2);
    B_est = exp(params_B(1)).*(set_of_mua(i_concentration)).^params_B(2);
    C_est = exp(params_C(1)).*(set_of_mua(i_concentration)).^params_C(2);
    the_dpf = A_est + B_est./d + C_est.*log(d)./d;
    plot(d,(OD-G)./the_dpf./d,'LineStyle','-','LineWidth',2,'Color',the_color,'DisplayName',num2str(set_of_concentration(i_concentration))), hold on
    plot(d,set_of_mua(i_concentration).*ones(size(d)),'LineStyle','-.','LineWidth',1,'Color',the_color,'HandleVisibility','off'), hold on
    xlabel('d (cm)'), ylabel('\mu_a (cm^{-1})'), set(gca,'fontsize',16)
    title('Estimation of \mu_a using inv. dist. (assisted)'), grid off, axis tight square
    clearvars the_color A_est B_est C_est d OD G
end
clearvars i_concentration

% Estimate mu_a when the true absorption coefficient is unknown
set_of_estimated_c_using_true_non_assisted = nan(length(set_of_concentration),1);
set_of_estimated_c_using_invd_non_assisted = nan(length(set_of_concentration),1);
set_of_estimated_mu_using_true_non_assisted = nan(length(set_of_concentration),1);
set_of_estimated_mu_using_invd_non_assisted = nan(length(set_of_concentration),1);
for i_c = 1:length(set_of_concentration)
    [~,idx_dist] = min(abs(squeeze(set_of_ODs(i_c,:,1))-SDS)); % disp(['SDS = ',num2str(SDS),' ~ ',num2str(squeeze(set_of_ODs(i_c,idx_dist,1))),' = SDS'])
    OD = -log(squeeze(set_of_ODs(:,:,2)));
    G = OD(:,1);
    OD = OD - G;
    true_mu_estim_mu = abs( (OD(i_c,idx_dist)./set_of_ODs_d_mua(:,idx_dist)./squeeze(set_of_ODs(:,idx_dist,1))) - set_of_mua.' );
    % true_mu_estim_mu = abs( (OD(i_c,idx_dist)./set_of_ODs_d_mua(:,idx_dist)./squeeze(set_of_ODs(:,idx_dist,1)))./set_of_mua.' - 1 );
    [~,idx] = min(true_mu_estim_mu);
    set_of_estimated_mu_using_true_non_assisted(i_c) = set_of_mua(idx);
    set_of_estimated_c_using_true_non_assisted(i_c) = set_of_concentration(idx);
    clearvars true_mu_estim_mu idx
    A_est = exp(params_A(1)).*(set_of_mua.').^params_A(2);
    B_est = exp(params_B(1)).*(set_of_mua.').^params_B(2);
    C_est = exp(params_C(1)).*(set_of_mua.').^params_C(2);
    the_dpf = A_est + B_est./squeeze(set_of_ODs(:,idx_dist,1)) + C_est.*log(squeeze(set_of_ODs(:,idx_dist,1)))./squeeze(set_of_ODs(:,idx_dist,1));
    true_mu_estim_mu = abs( (OD(i_c,idx_dist)./the_dpf./squeeze(set_of_ODs(:,idx_dist,1))) - set_of_mua.' );
    % plot(true_mu_estim_mu), hold on
    % true_mu_estim_mu = abs( (OD(i_c,idx_dist)./the_dpf./squeeze(set_of_ODs(:,idx_dist,1)))./set_of_mua.' - 1 );
    % plot(true_mu_estim_mu), hold off
    [~,idx] = min(true_mu_estim_mu);
    set_of_estimated_mu_using_invd_non_assisted(i_c) = set_of_mua(idx);
    set_of_estimated_c_using_invd_non_assisted(i_c) = set_of_concentration(idx);
    clearvars true_mu_estim_mu idx
    clearvars idx_dist OD G the_dpf A_est B_est C_est
end
clearvars i_c

figure(5), subplot(1, 2, 1)
plot(set_of_mua,set_of_estimated_mu_using_true_non_assisted,'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',16,'DisplayName','data','Color','r'), hold on
plot([0 34],[0 34],'k--','DisplayName','1-to-1'), hold off
xlabel('true \mu_a (cm^{-1})'), ylabel('estimated \mu_a (cm^{-1})'), set(gca,'fontsize',16), hold off
title('true vs. estimated \mu_a non-assist. true DPF'), grid on, axis square
axis([0 5 0 5])
legend('show','Location','best')
figure(5), subplot(1, 2, 2)
plot(set_of_mua,set_of_estimated_mu_using_invd_non_assisted,'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',16,'DisplayName','data','Color','r'), hold on
plot([0 34],[0 34],'k--','DisplayName','1-to-1'), hold off
xlabel('true \mu_a (cm^{-1})'), ylabel('estimated \mu_a (cm^{-1})'), set(gca,'fontsize',16), hold off
title('true vs. estimated \mu_a non-assist. inv. dist. DPF'), grid on, axis square
axis([0 5 0 5])
legend('show','Location','best')

figure(6), subplot(1, 2, 1)
plot(set_of_concentration,set_of_estimated_c_using_true_non_assisted,'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',16,'DisplayName','data','Color','r'), hold on
plot([0 0.5],[0 0.5],'k--','DisplayName','1-to-1'), hold off
xlabel('true concentration (% v/v)'), ylabel('estimated concentration (% v/v)'), set(gca,'fontsize',16), hold off
title('true vs. estimated concentration non-assist. true DPF'), grid on, axis square
axis([0 0.5 0 0.5])
legend('show','Location','best')
figure(6), subplot(1, 2, 2)
plot(set_of_concentration,set_of_estimated_c_using_invd_non_assisted,'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',16,'DisplayName','data','Color','r'), hold on
plot([0 0.5],[0 0.5],'k--','DisplayName','1-to-1'), hold off
xlabel('true concentration (% v/v)'), ylabel('estimated concentration (% v/v)'), set(gca,'fontsize',16), hold off
title('true vs. estimated concentration non-assist. inv. dist. DPF'), grid on, axis square
axis([0 0.5 0 0.5])
legend('show','Location','best')
end

function [y_bind,x_bind] = get_OD (the_filename,dlta_d)
Lx = 19.1; Ly = 19.1; Lz = 6.0;  % The size of the computational domain (in cm) to prevent reflections from the boundaries.
t_db = load(the_filename); clearvars the_filename

% % removing noisy points
% u = round(-log(t_db.w)./t_db.s,4); u_unique = unique(u); freq = nan(size(u_unique));
% for i_f = 1:length(u_unique), freq(i_f) = sum(u==u_unique(i_f)); end
% [~,i_fx] = max(freq);
% t_db.x_in = t_db.x_in(u==u_unique(i_fx));
% t_db.y_in = t_db.y_in(u==u_unique(i_fx));
% t_db.z_in = t_db.z_in(u==u_unique(i_fx));
% t_db.x_ot = t_db.x_ot(u==u_unique(i_fx));
% t_db.y_ot = t_db.y_ot(u==u_unique(i_fx));
% t_db.z_ot = t_db.z_ot(u==u_unique(i_fx));
% t_db.s = t_db.s(u==u_unique(i_fx));
% t_db.w = t_db.w(u==u_unique(i_fx));
% disp(num2str([length(u) sum(u==u_unique(i_fx)) sum(u==u_unique(i_fx))/length(u)*100]))
% t_db.no_of_photons = t_db.no_of_photons-sum(u~=u_unique(i_fx));
% clearvars u u_unique freq i_f i_fx

% 1-D sorting
d_diff_edges = 0:dlta_d:sqrt((Lx/2).^2+(Ly/2).^2);
idx = t_db.z_ot<=0;
[~,~,index_in] = histcounts(sqrt( ...
    (t_db.x_in(idx)-t_db.x_ot(idx)).^2 + ...
    (t_db.y_in(idx)-t_db.y_ot(idx)).^2 + ...
    (t_db.z_in(idx)-t_db.z_ot(idx)).^2),d_diff_edges); % d_diffuse bins
clearvars d_diff_edges Lx Ly Lz

% I vs. d
fun_x = @mean; fun_y = @sum;
x_temp = sqrt( ...
    (t_db.x_in(idx)-t_db.x_ot(idx)).^2 + ...
    (t_db.y_in(idx)-t_db.y_ot(idx)).^2 + ...
    (t_db.z_in(idx)-t_db.z_ot(idx)).^2);
y_temp = t_db.w(idx);
x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
y_bind = accumarray(index_in,y_temp,[],fun_y,nan);
y_bind = ...
    (            y_bind./(2.*pi.*dlta_d.*x_bind)) ./ ...
    (t_db.no_of_photons./(   pi.*dlta_d.*dlta_d));
end
function [params,err] = fitPlaneModel(X, Z)
% Fits Z = C + a*X to gridded data
A = [ones(numel(X), 1), X(:)];  % Design matrix: [1 X]
b = Z(:);                         % Observed values
% Solve for [C; a; b]
params = A \ b;
err = rsquared(b,A*params);
end
function [r2] = rsquared(y,yhat)
y    = y(:);
yhat = yhat(:);
idx = ~isnan(y)&~isnan(yhat)&~isinf(y)&~isinf(yhat);
y    = y(idx);
yhat = yhat(idx);
SS_res = sum((y - yhat).^2);    % Residual sum of squares
SS_tot = sum((y - mean(y)).^2); % Total sum of squares
r2 = 1 - (SS_res/SS_tot);
end
