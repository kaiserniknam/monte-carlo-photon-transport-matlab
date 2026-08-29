function [] = Photon_26 ()
% Repository group: 02_optical_properties_and_kramers_kronig
% Version role: primary simulation or study entry point
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Arnab's paper
% Designing refractive index fluids using the Kramers–Kronig relations
% trying to implement Kramers-Kronig, a sanity check for Fig. 1

clc
close all
set(0,'DefaultFigureWindowStyle','docked')

% load data
d_lmbd = 1; lmbd_s = 1:d_lmbd:1000;
kapa_s = ones(size(lmbd_s)).*eps; kapa_s(100<=lmbd_s&lmbd_s<=400) = 2.5;
RI___n = nan(size(lmbd_s));
for i_lmbda = 1:length(lmbd_s)
    s = (2/pi).*(kapa_s./lmbd_s).*(1./(1-(lmbd_s./lmbd_s(i_lmbda)).^2));
    s = s(~isinf(s));
    RI___n(i_lmbda) = sum(s).*d_lmbd; clearvars s
end
clearvars i_lmbda db
% Plotting absorbed photon positions (optional, for visualization)
figure(1), subplot(1,2,1)
plot(lmbd_s,kapa_s,'LineWidth',2,'DisplayName','\Kappa','Color','b')
xlabel('wavelength (nm)'), ylabel('\kappa'), title('Given')
set(gca,'fontsize',24), axis tight, axis square, grid on, hold on
set(gca,'xtick',[0 300 600 900]), xlim([0 1000])
set(gca,'ytick',[0 1 2 3]), ylim([0 2.75])
figure(1), subplot(1,2,2)
plot(lmbd_s,RI___n,'LineWidth',2,'DisplayName','n','Color','b'), hold on
plot(lmbd_s,zeros(size(lmbd_s)),'LineWidth',2,'DisplayName','zero','Color','k'), hold on
xlabel('wavelength (nm)'), ylabel('n'), title('Calculated')
set(gca,'fontsize',24), axis tight, axis square, grid on, hold on
set(gca,'xtick',[0 300 600 900]), xlim([0 1000])
set(gca,'ytick',[-4 0 +4 +8]), ylim([-4 +8])

clear
% load data
db = readtable('DB/Tianqi_Sai_Fig2.csv');
lmbd_s = db.Var1;
kapa_s = db.Var2;
d_lmbd = (mean(diff(lmbd_s)));
RI___n = nan(size(lmbd_s));
for i_lmbda = 1:length(lmbd_s)
    s = (2/pi).*(kapa_s./lmbd_s).*(1./(1-(lmbd_s./lmbd_s(i_lmbda)).^2));
    s = s(~isinf(s));
    RI___n(i_lmbda) = sum(s).*d_lmbd; clearvars s
end
clearvars i_lmbda db
% Plotting absorbed photon positions (optional, for visualization)
figure(2), subplot(1,2,1)
plot(lmbd_s,kapa_s,'LineWidth',2,'DisplayName','\Kappa','Color','b')
xlabel('wavelength (nm)'), ylabel('\kappa'), title('Given')
set(gca,'fontsize',24), axis tight, axis square, grid on, hold on
figure(2), subplot(1,2,2)
plot(lmbd_s,RI___n,'LineWidth',2,'DisplayName','n','Color','b')
xlabel('wavelength (nm)'), ylabel('n'), title('Calculated')
set(gca,'fontsize',24), axis tight, axis square, grid on, hold on
end
