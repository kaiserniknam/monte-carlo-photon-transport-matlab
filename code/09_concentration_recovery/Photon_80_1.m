function [] = Photon_80_1 ()
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
% Based on Version 79, a better absorption coefficient (mu_a)
% data collection

clc
close all

set_of_verse = 0:9;
set_of_concentration = [nan, 0.00:0.01:0.15]; % set of concentrations
for i_concentration = 1:length(set_of_concentration)

    x_in = []; y_in = []; z_in = [];
    x_ot = []; y_ot = []; z_ot = [];
    s = []; w = [];
    no_of_photons = 0;
    M_raw = [];
    the_concentration = set_of_concentration(i_concentration);
    mua_phnt = 0; mus_phnt = 0; g_phnt = 0; n_phnt = 0;
    mua_incl = 0; mus_incl = 0; g_incl = 0; n_incl = 0;

    for i_rep = set_of_verse
        the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_80_concentration_',num2str(the_concentration),'_',num2str(i_rep),'.mat'];
        db = load(the_filename);

        x_in = [x_in ; db.x_in]; y_in = [y_in ; db.y_in]; z_in = [z_in ; db.z_in];
        x_ot = [x_ot ; db.x_ot]; y_ot = [y_ot ; db.y_ot]; z_ot = [z_ot ; db.z_ot];
        s = [s ; db.s];
        w = [w ; db.w];
        no_of_photons = no_of_photons + db.no_of_photons;

        M_raw = db.M_raw;
        mua_phnt = db.mua_phnt; mus_phnt = db.mus_phnt; g_phnt = db.g_phnt; n_phnt = db.n_phnt;
        mua_incl = db.mua_incl; mus_incl = db.mus_incl; g_incl = db.g_incl; n_incl = db.n_incl;

        clearvars the_filename db
    end
    clearvars i_rep

    save  (['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_80_concentration_',num2str(the_concentration),'.mat'], ...
        'x_in', 'y_in', 'z_in', 'x_ot', 'y_ot', 'z_ot','s', 'w', 'no_of_photons', 'M_raw', ...
        'mua_phnt', 'mus_phnt', 'g_phnt', 'n_phnt', ...
        'mua_incl', 'mus_incl', 'g_incl', 'n_incl')

        clearvars x_in y_in z_in x_ot y_ot z_ot s w no_of_photons M_raw
        clearvars mua_phnt mus_phnt g_phnt n_phnt
        clearvars mua_incl mus_incl g_incl n_incl
        clearvars the_concentration
end
end
