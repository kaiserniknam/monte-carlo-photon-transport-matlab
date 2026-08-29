function [] = Photon_25 ()
% Repository group: 02_optical_properties_and_kramers_kronig
% Version role: primary simulation or study entry point
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Arnab's paper
% Designing refractive index fluids using the Kramers–Kronig relations
% trying to implement Kramers-Kronig

clc
close all
set(0,'DefaultFigureWindowStyle','docked')

% load data
db = readtable('DB/refrective index_example.csv');
lmbd_s = db.Var1;
% alfa_s = db.Var2; kapa_s = alfa_s.*lmbd_s/(4*pi);
kapa_s = db.Var2;
d_lmbd = (mean(diff(lmbd_s(11:end))));
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
figure(1), subplot(1,2,2)
plot(lmbd_s,RI___n,'LineWidth',2,'DisplayName','n','Color','b')
xlabel('wavelength (nm)'), ylabel('n'), title('Calculated')
set(gca,'fontsize',24), axis tight, axis square, grid on, hold on
end
