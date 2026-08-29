function [] = Photon_65_1 ()
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
% Based on Version 29, adapted for TiO2 experimental dataset.
% analysis data

clc
close all
set(0,'DefaultFigureWindowStyle','docked')

set_of_concentration = [nan, 0.00:0.01:0.5]; % set of concentrations
set_of_ODs = nan(length(set_of_concentration),250,2);
for i_concentration = 1:length(set_of_concentration)
    the_concentration = set_of_concentration(i_concentration);
    the_color = [(i_concentration-1)./(length(set_of_concentration)-1) 0 1-(i_concentration-1)./(length(set_of_concentration)-1)];

    the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_65_concentration_',num2str(the_concentration),'.mat' ];
    [y_bind,x_bind] = get_OD (the_filename,0.17);
    set_of_ODs(i_concentration,1:length(x_bind),1) = x_bind;
    set_of_ODs(i_concentration,1:length(y_bind),2) = y_bind;

    figure(1), subplot(1,2,1), plot(x_bind,-log(y_bind),'LineStyle','-','LineWidth',2,'Color',the_color,'DisplayName',num2str(the_concentration)), hold on
    xlabel('d (cm)'), ylabel('OD (a.u.)'), set(gca,'fontsize',16)
    title('Effect of Inclusion Concentration on Optical Density'), grid off, axis tight square
    clearvars the_concentration the_color the_filename y_bind x_bind
end
set_of_mua_s = 0.01 + 37.*set_of_concentration;
clearvars i_concentration

dims = axis;
set_of_ODs_at_3cm = nan(length(set_of_concentration),1);
for i_concentration = 1:length(set_of_concentration)
    [~,idx] = min(abs(set_of_ODs(i_concentration,:,1)-3));
    figure(1), subplot(1,2,1), plot([set_of_ODs(i_concentration,idx,1) set_of_ODs(i_concentration,idx,1)],[dims(3) dims(4)],'LineStyle','-.','LineWidth',2,'Color','k','HandleVisibility','off'),
    set_of_ODs_at_3cm(i_concentration) = set_of_ODs(i_concentration,idx,2);
    clearvars idx
end
figure(1), subplot(1,2,2), plot(set_of_concentration,-log(set_of_ODs_at_3cm),'LineStyle','-','LineWidth',2,'Color','k'),
xlabel('concentration (%v/v)'), ylabel('OD (a.u.)'), set(gca,'fontsize',16)
title('Effect of Inclusion Concentration on Optical Density at SDS = 3 cm'), grid on, axis square
clearvars dims set_of_ODs_at_3cm i_concentration



set_of_dpfs = nan(length(set_of_concentration),250,2);
set_of_params = nan(length(set_of_concentration),3);
for i_concentration = 2:length(set_of_concentration)
    x_bind =      squeeze(set_of_ODs(i_concentration,:,1)).';
    y_bind = -log(squeeze(set_of_ODs(i_concentration,:,2))).';
    G = y_bind(1);

    % model fitting
    x = x_bind(~isnan(x_bind) & ~isnan(y_bind) & x_bind > 0);
    y = y_bind(~isnan(x_bind) & ~isnan(y_bind) & x_bind > 0);
    Y = (y - G) ./ x ./ set_of_mua_s(i_concentration); % response #1
    Y = (y - G) ./ x;                                  % response #2
    set_of_dpfs(i_concentration,1:length(x),1) = x;    % save DPF's - x
    set_of_dpfs(i_concentration,1:length(Y),2) = Y;    % save DPF's - dpf
    invX = 1 ./ x;                                     % predictor 1
    logX_over_X = log(x) ./ x;                         % predictor 2
    T = table(Y, invX, logX_over_X, 'VariableNames', {'Y','invX','logX_over_X'});
    mdl = fitlm(T, 'Y ~ 1 + invX + logX_over_X');      % A, B, C
    clearvars T x y Y invX logX_over_X
    set_of_params(i_concentration,1) = mdl.Coefficients.Estimate(1);
    set_of_params(i_concentration,2) = mdl.Coefficients.Estimate(2);
    set_of_params(i_concentration,3) = mdl.Coefficients.Estimate(3);
    set_of_params(i_concentration,4) = mdl.Coefficients.pValue(2);
    clearvars x_bind y_bind the_mua the_color mdl G
end
clearvars i_concentration

[params_A,err_A] = fitPlaneModel(log(set_of_mua_s(2:end)), log(set_of_params(2:end,1)));
A_est = exp(params_A(1)).*((set_of_mua_s).^params_A(2));
[params_B,err_B] = fitPlaneModel(log(set_of_mua_s(2:end)), log(set_of_params(2:end,2)));
B_est = exp(params_B(1)).*((set_of_mua_s).^params_B(2));
[params_C,err_C] = fitPlaneModel(log(set_of_mua_s(2:end)), log(set_of_params(2:end,3)));
C_est = exp(params_C(1)).*((set_of_mua_s).^params_C(2));

figure(2), subplot(1,2,1)
plot(set_of_mua_s,set_of_params(:,1),'b','LineWidth',2,'LineStyle','-', 'DisplayName','A'), hold on
plot(set_of_mua_s,A_est,             'b','LineWidth',2,'LineStyle','-.','DisplayName',['estim. A, r^2 ~ ',num2str(err_A)]), hold on
plot(set_of_mua_s,set_of_params(:,2),'r','LineWidth',2,'LineStyle','-', 'DisplayName','B'), hold on
plot(set_of_mua_s,B_est,             'r','LineWidth',2,'LineStyle','-.','DisplayName',['estim. B, r^2 ~ ',num2str(err_B)]), hold on
plot(set_of_mua_s,set_of_params(:,3),'g','LineWidth',2,'LineStyle','-', 'DisplayName','C'), hold on
plot(set_of_mua_s,C_est,             'g','LineWidth',2,'LineStyle','-.','DisplayName',['estim. C, r^2 ~ ',num2str(err_C)]), hold on
xlabel('\mu_a (cm^{-1})');
ylabel('Inverse-distance parameters'); set(gca,'YScale','log')
set(gca,'FontSize',16);
title('Effect of Inclusion Absorption on Optical Density (SDS = 3 cm)');
grid on; axis square; legend('show','Location','northeast')
disp(num2str(['C_A = ',num2str(exp(params_A(1)),'%.2f'),', p_A = ',num2str(params_A(2),'%.2f'),', @ r2 =  ',num2str(err_A)])), clearvars err params
disp(num2str(['C_B = ',num2str(exp(params_B(1)),'%.2f'),', p_B = ',num2str(params_B(2),'%.2f'),', @ r2 =  ',num2str(err_B)])), clearvars err params
disp(num2str(['C_C = ',num2str(exp(params_C(1)),'%.2f'),', p_C = ',num2str(params_C(2),'%.2f'),', @ r2 =  ',num2str(err_C)])), clearvars err params
clearvars A_est B_est C_est
clearvars err_A err_B err_C

figure(2), subplot(1,2,2)
for i_concentration = 2:length(set_of_concentration)
    the_color = [(i_concentration-1)./(length(set_of_concentration)-1) 0 1-(i_concentration-1)./(length(set_of_concentration)-1)];
    A_est = exp(params_A(1)).*(set_of_mua_s(i_concentration)).^params_A(2);
    B_est = exp(params_B(1)).*(set_of_mua_s(i_concentration)).^params_B(2);
    C_est = exp(params_C(1)).*(set_of_mua_s(i_concentration)).^params_C(2);

    the_dpf = A_est + B_est./squeeze(set_of_ODs(i_concentration,:,1)) + C_est.*log(squeeze(set_of_ODs(i_concentration,:,1))./squeeze(set_of_ODs(i_concentration,:,1)));
    plot(squeeze(set_of_dpfs(i_concentration,:,1)),squeeze(set_of_dpfs(i_concentration,:,2)),'LineStyle','-','LineWidth',2,'Color',the_color,'DisplayName',['C = ',num2str(set_of_concentration(i_concentration))]), hold on
    plot(squeeze(set_of_dpfs(i_concentration,:,1)),the_dpf,'LineStyle','-.','LineWidth',2,'Color',the_color,'HandleVisibility','off'), hold on
    xlabel('d (cm)'), ylabel('DPF (a.u.)'), set(gca,'fontsize',16), set(gca,'YScale','log')
    title('Effect of Inclusion Concentration on Optical Density'), grid off, axis tight square
    clearvars the_color A_est B_est C_est the_dpf
end
% clearvars params_A params_B params_C i_concentration



figure(3), subplot(1,2,1)
for i_concentration = 2:length(set_of_concentration)
    the_color = [(i_concentration-1)./(length(set_of_concentration)-1) 0 1-(i_concentration-1)./(length(set_of_concentration)-1)];
    G = -log(squeeze(set_of_ODs(i_concentration,1,2)));
    A_est = exp(params_A(1)).*(set_of_mua_s(i_concentration)).^params_A(2);
    B_est = exp(params_B(1)).*(set_of_mua_s(i_concentration)).^params_B(2);
    C_est = exp(params_C(1)).*(set_of_mua_s(i_concentration)).^params_C(2);
    the_dpf = A_est + B_est./squeeze(set_of_dpfs(i_concentration,:,1)) + C_est.*log(squeeze(set_of_dpfs(i_concentration,:,1)))./squeeze(set_of_dpfs(i_concentration,:,1));
    plot(squeeze(set_of_ODs(i_concentration,:,1)),(-log(squeeze(set_of_ODs(i_concentration,:,2)))-G)./the_dpf./squeeze(set_of_ODs(i_concentration,:,1)),'LineStyle','-','LineWidth',2,'Color',the_color,'DisplayName',num2str(set_of_concentration(i_concentration))), hold on
    plot(squeeze(set_of_ODs(i_concentration,:,1)),set_of_mua_s(i_concentration).*ones(size(squeeze(set_of_ODs(i_concentration,:,1)))),'LineStyle','-.','LineWidth',1,'Color',the_color,'HandleVisibility','off'), hold on
    xlabel('d (cm)'), ylabel('\mu_a (cm^{-1})'), set(gca,'fontsize',16)
    title('Estimation of Concentration'), grid off, axis tight square
    clearvars the_color A_est B_est C_est the_dpf
end
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
% Flatten everything to column vectors
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
