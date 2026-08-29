function [] = Photon_88()
% Repository group: 11_validation_and_calibration
% Version role: primary simulation or study entry point
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Monte Carlo simulations for Arnab's 1% Intralipid + India Ink phantom experiment.
% Purpose: estimate mu_s by comparing simulated OD vs. source-detector separation with experiment at known mu_a and g.
% This version sweeps mu_s while keeping mu_a, g, refractive index, geometry, and source configuration fixed.

clc
close all

%% Simulation parameters
z_air = 0.0;                       % Air-layer thickness (cm)
nPhotonsReq = 1e5;                % Number of photons requested per Monte Carlo run
nExamplePaths = 1000000;          % Total number of detected photon paths to collect

%% Computational domain and source geometry
Lx = 29.1; Ly = 29.1; Lz = 06.0; % Computational-domain dimensions (cm)
beam_phi = 0;                     % Beam azimuthal angle (rad)
beam_tht = 0;                     % Beam polar angle (rad)
beam_X = 0.0;                     % Beam x-position (cm)
beam_Y = 0.0;                     % Beam y-position (cm)

%% Optical properties
n = 1.33; g = 0.90;               % Refractive index and anisotropy factor
set_of_mua = [0.4];               % Absorption coefficients to simulate (cm^-1)
set_of_mus = 1:100;               % Scattering coefficients to simulate (cm^-1)

%% Run simulations over optical-property combinations
for i_a = 1:length(set_of_mua)
    for i_s = 1:length(set_of_mus)

        mua = set_of_mua(i_a);
        mus = set_of_mus(i_s);

        the_filename = ['Photon_88_mua_',sprintf('%.4f',set_of_mua(i_a)),'_mus_',sprintf('%.2f',set_of_mus(i_s)),'.mat' ];

        % Run only if the corresponding simulation file does not already exist
        if ~exist(the_filename,'file')

            [x, y, z, d, s, w, c, a, no_of_photons, M_raw] = do_simulation ( ...
                [mua,mus,g   ,n  ], ...
                [nan,nan,nan ,nan], ...
                [Lx,Ly,Lz],[nan,nan,nan],[nan,nan,nan],z_air,001,nPhotonsReq,nExamplePaths,beam_X,beam_Y,beam_phi,beam_tht);

            % Save detected-photon information and the raw Monte Carlo result
            save(the_filename,'x','y','z','d','s','w','c','a','no_of_photons','M_raw')

            clearvars x y z d s w c a no_of_photons M_raw
            clearvars mua mus
        end

        clearvars the_filename
    end
end
end


function [x, y, z, d, s, w, c, a, no_of_photons, M_raw] = do_simulation(opt_bckg,opt_sgnl,cmd_size,sig_pos,sig_size,zSurface,wvlngth,nPhotonsReq,nExamplePaths,beam_X,beam_Y,beam_phi,beam_tht)

%% Computational geometry
model = MCmatlab.model;
dl = 0.1;                                            % Spatial resolution (cm)
model.G.nx = round(cmd_size(1)./dl);                 % Number of voxels along x
model.G.ny = round(cmd_size(2)./dl);                 % Number of voxels along y
model.G.nz = round(cmd_size(3)./dl);                 % Number of voxels along z
model.G.Lx = cmd_size(1);                            % Domain size along x (cm)
model.G.Ly = cmd_size(2);                            % Domain size along y (cm)
model.G.Lz = cmd_size(3);                            % Domain size along z (cm)

model.G.mediaPropertiesFunc = @mediaPropertiesFunc;  % Function defining optical properties
model.G.mediaPropParams = {opt_bckg(1),opt_bckg(2),opt_bckg(3),opt_bckg(4),opt_sgnl(1),opt_sgnl(2),opt_sgnl(3),opt_sgnl(4)};
model.G.geomFunc = @geometryDefinition;              % Function defining the spatial distribution of media
model.G.geomFuncParams = {zSurface,2,3,sig_pos(1),sig_pos(2),sig_pos(3),sig_size(1),sig_size(2),sig_size(3)};
% model = plot(model,'G');

clearvars opt_bckg opt_sgnl cmd_size sig_pos sig_size dl

%% Monte Carlo settings
model.MC.simulationTimeRequested = 0.25;             % Requested simulation time (min)
model.MC.matchedInterfaces = false;                  % Account for refractive-index mismatch at interfaces
model.MC.boundaryType = 1;                           % All six boundaries are escaping
model.MC.nPhotonsRequested = nPhotonsReq;            % Number of photons requested per run
model.MC.nExamplePaths = nPhotonsReq;                % Number of photon paths stored by MCmatlab
model.MC.wavelength = wvlngth;                       % Wavelength (nm)

clearvars wvlngth nPhotonsReq

%% Pencil-beam source
model.MC.lightSource.sourceType = 0;                  % Pencil beam
model.MC.silentMode = false;
model.MC.lightSource.xFocus = beam_X;                 % Beam focus x-position (cm)
model.MC.lightSource.yFocus = beam_Y;                 % Beam focus y-position (cm)
model.MC.lightSource.zFocus = zSurface;               % Beam focus z-position (cm)
model.MC.lightSource.theta = beam_tht;                % Beam polar angle (rad)
model.MC.lightSource.phi = beam_phi;                  % Beam azimuthal angle (rad)

clearvars zSurface beam_phi beam_tht

%% Run Monte Carlo simulations until the requested number of photon paths is collected
scatter_stat = [];
no_of_photons = 0;

while size(scatter_stat,1) < nExamplePaths

    t_model = runMonteCarlo(model);
    M_raw = t_model.G.M_raw;

    % Identify the beginning and end of each stored photon trajectory
    spratrs = find(isnan(t_model.MC.examplePaths(1,:)));
    strt_pnts = spratrs(2:2:end); strt_pnts = strt_pnts(1:end-1);
    fnsh_pnts = spratrs(1:2:end); fnsh_pnts = fnsh_pnts(2:end-0);

    % Columns: x, y, z, pathlength, weight, dx, dy, dz, exit-surface code
    t_scatter_stat = nan(length(fnsh_pnts),3+1+1+3+1);
    clearvars spratrs

    for idx = 1:length(strt_pnts)

        Photon_Path = t_model.MC.examplePaths(:,strt_pnts(idx)+1:fnsh_pnts(idx)-1);

        % Photon exits through the top surface: diffuse reflectance
        if Photon_Path(3,end) <= 0.0

            t_scatter_stat(idx,:) = [...
                Photon_Path(1,end), ...
                Photon_Path(2,end), ...
                Photon_Path(3,end), ...
                sum(sqrt(sum((Photon_Path(1:3,2:end-0)-Photon_Path(1:3,1:end-1)).^2,1))), ...
                Photon_Path(4,end), ...
                (Photon_Path(1,end)-Photon_Path(1,end-1)), ...
                (Photon_Path(2,end)-Photon_Path(2,end-1)), ...
                (Photon_Path(3,end)-Photon_Path(3,end-1)), ...
                (0), ...
                ];

        % Photon exits through the bottom surface: transmission
        elseif Photon_Path(3,end) >= t_model.G.Lz

            t_scatter_stat(idx,:) = [...
                Photon_Path(1,end), ...
                Photon_Path(2,end), ...
                Photon_Path(3,end), ...
                sum(sqrt(sum((Photon_Path(1:3,2:end-0)-Photon_Path(1:3,1:end-1)).^2,1))), ...
                Photon_Path(4,end), ...
                (Photon_Path(1,end)-Photon_Path(1,end-1)), ...
                (Photon_Path(2,end)-Photon_Path(2,end-1)), ...
                (Photon_Path(3,end)-Photon_Path(3,end-1)), ...
                (1), ...
                ];

        % Ignore photons exiting through the lateral boundaries
        else
        end

        clearvars Photon_Path
    end

    clearvars strt_pnts fnsh_pnts idx t_model

    % Keep only photons exiting through the top or bottom surfaces
    surf_photons = ~isnan(t_scatter_stat(:,1));
    no_of_photons = no_of_photons + size(t_scatter_stat,1);

    disp('****************************************************')
    disp(['progress ... ',num2str(size(scatter_stat,1)/nExamplePaths*100),'%'])
    disp('****************************************************')

    t_scatter_stat = t_scatter_stat(surf_photons,:);
    scatter_stat = [scatter_stat ; t_scatter_stat];

    clearvars surf_photons t_scatter_stat
end

%% Extract detected-photon quantities
x = scatter_stat(:,1);                                                 % Photon exit x-coordinate (cm)
y = scatter_stat(:,2);                                                 % Photon exit y-coordinate (cm)
z = scatter_stat(:,3);                                                 % Photon exit z-coordinate (cm)
d = sqrt((scatter_stat(:,1)-beam_X).^2+(scatter_stat(:,2)-beam_Y).^2); % Radial source-detector separation (cm)
s = scatter_stat(:,4);                                                 % Total photon pathlength (cm)
w = scatter_stat(:,5);                                                 % Photon weight
a = scatter_stat(:,8)./sqrt(scatter_stat(:,6).^2 + scatter_stat(:,7).^2 + scatter_stat(:,8).^2); % Cosine of exit angle relative to z-axis
c = scatter_stat(:,9);                                                 % Exit code: 0 = diffuse reflectance, 1 = transmission
end
function M = geometryDefinition(X,Y,Z,parameters)
% Define the spatial distribution of the simulation media.

M = ones(size(X));                                                     % Initialize all voxels as air
zSurface = double(parameters{1});                                      % Top surface of the phantom (cm)
M(Z>=zSurface) = double(parameters{2});                                % Assign homogeneous phantom medium

M(sqrt( ...
    ((X-double(parameters{4})         ).^2)./(double(parameters{7}).^2) + ...
    ((Y-double(parameters{5})         ).^2)./(double(parameters{8}).^2) + ...
    ((Z-double(parameters{6})-zSurface).^2)./(double(parameters{9}).^2)) ...
    <=1) = double(parameters{3});                                      % Assign inclusion medium when present

M = uint8(M);
end
function mediaProperties = mediaPropertiesFunc(var)
% Define the optical properties of each simulation medium.

mediaProperties = MCmatlab.mediumProperties;

j = 1;
mediaProperties(j).name = 'air';
mediaProperties(j).mua = 1e-4;                                        % Absorption coefficient (cm^-1)
mediaProperties(j).mus = 1e-5;                                        % Scattering coefficient (cm^-1)
mediaProperties(j).g = 1;                                             % Scattering anisotropy factor
mediaProperties(j).n = 1;                                             % Refractive index

j = 2;
mediaProperties(j).name = 'breast tissue';
mediaProperties(j).mua = double(var{1});                              % Absorption coefficient (cm^-1)
mediaProperties(j).mus = double(var{2});                              % Scattering coefficient (cm^-1)
mediaProperties(j).g = double(var{3});                                % Scattering anisotropy factor
mediaProperties(j).n = double(var{4});                                % Refractive index

j = 3;
mediaProperties(j).name = 'tumor tissue';
mediaProperties(j).mua = double(var{5});                              % Absorption coefficient (cm^-1)
mediaProperties(j).mus = double(var{6});                              % Scattering coefficient (cm^-1)
mediaProperties(j).g = double(var{7});                                % Scattering anisotropy factor
mediaProperties(j).n = double(var{8});                                % Refractive index
end
