function [] = Photon_81_3 ()
% Repository group: 07_dpf_parameter_studies
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% This version computes w and s as functions of d across varying mu_s and mu_a,
% based on Case 33 but using a different set of optical properties.
% analysis data: comapre with 58
% the similarity is due to g!

clc
close all


db33 = load('/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_58.mat');
db81 = load('/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_81.mat');

mua = 0.0275; mus = 45;
plot(squeeze(db33.set_of__ds(db33.set_of_mua==      mua   ,db33.set_of_mus==mus,:)),squeeze(db33.set_of_ODs(db33.set_of_mua==      mua   ,db33.set_of_mus==mus,:)),'DisplayName','58',LineWidth=2,Color='b'), hold on
plot(squeeze(db81.set_of__ds(db81.set_of_mua==round(mua,2),db81.set_of_mus==mus,:)),squeeze(db81.set_of_ODs(db81.set_of_mua==round(mua,2),db81.set_of_mus==mus,:)),'DisplayName','81',LineWidth=2,Color='r'), hold on
xlabel('SDS (cm)')
ylabel('OD (a.u.)')
set(gca,'fontsize',24')
grid on, legend('show','Location','southeast')
axis square
end
