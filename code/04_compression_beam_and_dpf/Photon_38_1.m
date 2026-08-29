function [] = Photon_38_1 ()
% Repository group: 04_compression_beam_and_dpf
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: s vs. d for breast tissue (Breast is emitting / no tumor)
% For Arnab research / SPIE WEST 2025; analysis

clc
close all

z_air = 0.0; % the thickness of air layer
Lx = 19.1; Ly = 19.1; Lz = 5.0; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
beam_phi = 0;     % polar angle of beam in cm
beam_tht = 0;     % azimuthal angle of beam in cm
beam_X = 0.0;     % X deviation of beam in cm
beam_Y = 0.0;     % Y deviation of beam in cm

db = load (['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_38.mat' ]);


% scematics
subplot(1,2,1),
cubeHalfSize_x = Lx/2;
cubeHalfSize_y = Ly/2;
cubeHalfSize_z = Lz/2;
% Cube vertices
vertices = [
    -cubeHalfSize_x -cubeHalfSize_y  0;
     cubeHalfSize_x -cubeHalfSize_y  0;
     cubeHalfSize_x  cubeHalfSize_y  0;
    -cubeHalfSize_x  cubeHalfSize_y  0;
    -cubeHalfSize_x -cubeHalfSize_y  -2*cubeHalfSize_z;
     cubeHalfSize_x -cubeHalfSize_y  -2*cubeHalfSize_z;
     cubeHalfSize_x  cubeHalfSize_y  -2*cubeHalfSize_z;
    -cubeHalfSize_x  cubeHalfSize_y  -2*cubeHalfSize_z;
];
% Cube faces
faces = [
    1 2 3 4; % Bottom
    5 6 7 8; % Top
    1 2 6 5; % Front
    2 3 7 6; % Right
    3 4 8 7; % Back
    4 1 5 8; % Left
];
patch('Vertices', vertices, 'Faces', faces, ...
      'FaceColor', 'b', 'EdgeColor', 'black', 'LineWidth', 1,'FaceAlpha', 0.1); hold on
% Draw sources
plot3(db.st(:,1),db.st(:,2),-db.st(:,3),'r.','MarkerSize',16), hold on
set(gca,'xtick',[round(min(db.x_c),1) 0 round(max(db.x_c),1)]), xlabel(['x-axis (cm)'])
set(gca,'ytick',[round(min(db.x_c),1) 0 round(max(db.x_c),1)]), ylabel(['y-axis (cm)'])
set(gca,'ztick',[-Lz 0]), zlabel(['z-axis (cm)']), title('The Schematics')
set(gca,'fontsize',12), axis equal, axis tight, grid on; view([35 5])
clearvars faces vertices cubeHalfSize_x cubeHalfSize_y cubeHalfSize_z

% scematics
set_of_N = 10000:1000:1000000;
set_of_I = nan(size(set_of_N));
for i_N = 1:length(set_of_N)
    idx = randperm(db.no_of_surf_photons,set_of_N(i_N));
    % set_of_I(i_N) = -log(sum(db.w(idx))/((db.no_of_surf_photons+db.no_of_dpth_photons)*set_of_N(i_N)/(db.no_of_surf_photons)));
    % set_of_I(i_N) = -log(sum(db.w(idx))/set_of_N(i_N));
    % set_of_I(i_N) = log(sum(db.w(idx)));
    set_of_I(i_N) = (sum(db.w(idx)))./(sum(db.w(:)));
end
subplot(1,2,2)
plot(set_of_N,set_of_I,'Color','k','LineStyle','-','LineWidth',2),
xlim([min(set_of_N)-5000 max(set_of_N)+5000]), xlabel('# of point sources (activated np)'),
ylim([0 1] )                                 , ylabel('\Sigma I/I_0 (a.u.)'),
title('Intemsity at surface vs. activated np'),
set(gca,'fontsize',16), axis square, grid on, hold on
end
