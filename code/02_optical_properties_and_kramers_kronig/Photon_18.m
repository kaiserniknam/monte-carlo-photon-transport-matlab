function [] = Photon_18
% Repository group: 02_optical_properties_and_kramers_kronig
% Version role: primary simulation or study entry point
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: plot mua & mus vs lambda --from selected references
% ref. #1: Assessment of the size, position, and optical properties of breast tumors in vivo by noninvasive optical methods
% ref. #2: Quantitative Absorption and Scattering Spectra in Thick Tissues Using Broadband Diffuse Optical Spectroscopy

clc
close all

% formula
brst_db = readtable('DB/breast.normal.csv');
tumr_db = readtable('DB/breast.tumor.csv');
% data
brst(1).lambda = 690; brst(1).mua = 0.030; brst(1).mus = 12.0; brst(1).g = 0.95; brst(1).n = 1.4;
brst(2).lambda = 825; brst(2).mua = 0.040; brst(2).mus = 11.0; brst(2).g = 0.95; brst(2).n = 1.4;
tumr(1).lambda = 690; tumr(1).mua = 0.084; tumr(1).mus = 15.0; tumr(1).g = 0.95; tumr(1).n = 1.4;
tumr(2).lambda = 825; tumr(2).mua = 0.085; tumr(2).mus = 12.7; tumr(2).g = 0.95; tumr(2).n = 1.4;
% fit
bs_brst = -log(brst(2).mus/brst(1).mus)/log(brst(2).lambda/brst(1).lambda);
as_brst = brst(1).mus/brst(1).lambda^(-bs_brst);
bs_tumr = -log(tumr(2).mus/tumr(1).mus)/log(tumr(2).lambda/tumr(1).lambda);
as_tumr = tumr(1).mus/tumr(1).lambda^(-bs_tumr);

subplot(1,2,1)
plot(brst_db.Var1,brst_db.Var2*10,'DisplayName','breast: model','LineStyle','-','Color','b','LineWidth',1.5), hold on
plot(brst(1).lambda,brst(1).mua,'Marker','s','MarkerFaceColor','b','MarkerEdgeColor','k','DisplayName','breast: data #1','MarkerSize',16,'LineStyle','none'), hold on
plot(brst(2).lambda,brst(2).mua,'Marker','o','MarkerFaceColor','b','MarkerEdgeColor','k','DisplayName','breast: data #2','MarkerSize',16,'LineStyle','none'), hold on
plot(tumr_db.Var1,tumr_db.Var2*10,'DisplayName','tumor: model', 'LineStyle','-','Color','r','LineWidth',1.5), hold on
plot(tumr(1).lambda,tumr(1).mua,'Marker','s','MarkerFaceColor','r','MarkerEdgeColor','k','DisplayName','tumor: data #1','MarkerSize',16,'LineStyle','none'), hold on
plot(tumr(2).lambda,tumr(2).mua,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','k','DisplayName','tumor: data #2','MarkerSize',16,'LineStyle','none'), hold on
xlabel('wavelength (\lambda, nm)'), ylabel('\mu_a (cm^{-1})'), title('\mu_a vs. wavelength')
axis tight, axis square, set(gca,'fontsize',24)
legend('show','Location','best')

subplot(1,2,2)
plot(brst_db.Var1,as_brst.*(brst_db.Var1.^(-bs_brst)),'DisplayName','breast: model','LineStyle','-','Color','b','LineWidth',1.5), hold on
plot(brst(1).lambda,brst(1).mus,'Marker','s','MarkerFaceColor','b','MarkerEdgeColor','k','DisplayName','breast: data #1','MarkerSize',16,'LineStyle','none'), hold on
plot(brst(2).lambda,brst(2).mus,'Marker','o','MarkerFaceColor','b','MarkerEdgeColor','k','DisplayName','breast: data #2','MarkerSize',16,'LineStyle','none'), hold on
plot(tumr_db.Var1,as_tumr.*(tumr_db.Var1.^(-bs_tumr)),'DisplayName','tumor: model', 'LineStyle','-','Color','r','LineWidth',1.5), hold on
plot(tumr(1).lambda,tumr(1).mus,'Marker','s','MarkerFaceColor','r','MarkerEdgeColor','k','DisplayName','tumor: data #1','MarkerSize',16,'LineStyle','none'), hold on
plot(tumr(2).lambda,tumr(2).mus,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','k','DisplayName','tumor: data #2','MarkerSize',16,'LineStyle','none'), hold on
xlabel('wavelength (\lambda, nm)'), ylabel('\mu''_s (cm^{-1})'), title('\mu''_s vs. wavelength')
axis tight, axis square, set(gca,'fontsize',24)
legend('show','Location','best')
end
