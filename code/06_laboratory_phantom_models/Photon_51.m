function [] = Photon_51 ()
% Repository group: 06_laboratory_phantom_models
% Version role: primary simulation or study entry point
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: simulating phantoms prepared by Arnab
% the same as Photon_47, but different optical properties for India Ink
% (constant mus, lower mua, higher g)(taken from Photon_50)

clc
close all

wvlnt = 760; % sample wavelength (nm)
nPhotonsReq = 1e5;  % number of requested photon
nExamplePaths = 1000000; % number of example photon paths
% Optical properties of TiO2 in agar
mua_TiO2 = .003; mus_TiO2 = 16; g_TiO2 = 0.88; n_TiO2 = 1.36;
% Optical properties of India Ink
% Spinelli, et al. Determination of reference values for optical properties of liquid phantoms based on Intralipid and India ink
mua_indk = 3.24; mus_indk = 0.1; g_indk = 0.9; n_indk = 1.34;
set_of_TiO2 = [0.50, 2.50, 5.0, 7.5]; % percents of TiO2
set_of_indk = [0.01, 0.05, 0.25, 1.0, 5.0, 10.0, 20.0]; % percents of India ink
dl = 0.02; Nz = 10/dl; Ny = 10/dl; Nx = 5/dl;
the_Tumr_X = 0; the_Tumr_Y = 0;   the_Tumr_Z = 1;   % Tumor center location all in cm
the_Tumr_L = 4; the_Tumr_W =.3/2; the_Tumr_H =.3/2; % Tumor size all in cm
the_beam_X = 3; the_beam_Y = 0;   the_beam_Z = 0;   % Source location in cm

for i_TiO2 = 1:length(set_of_TiO2)
    for i_indk = 1:length(set_of_indk)
        the_TiO2 = set_of_TiO2(i_TiO2);
        the_indk = set_of_indk(i_indk);

        T_mua_TiO2 = mua_TiO2;
        T_mus_TiO2 = mus_TiO2*the_TiO2;
        T_g_TiO2   = g_TiO2;
        T_n_TiO2   = n_TiO2;

        T_mua_indk = mua_indk*the_indk;
        T_mus_indk = mus_indk;
        T_g_indk   = g_indk;
        T_n_indk   = n_indk;

        % calculate optical properties of tumor
        [p_in, w_in, p_ot, w_ot, s, dl, no_of_photons] = do_simulation ( ...
            [Nz, Ny, Nx], ...
            [T_mua_TiO2,T_mus_TiO2,T_g_TiO2,T_n_TiO2], ...
            [T_mua_indk,T_mus_indk,T_g_indk,T_n_indk], ...
            dl,wvlnt,nPhotonsReq,nExamplePaths,...
            the_Tumr_X,the_Tumr_Y,the_Tumr_Z, ...
            the_Tumr_L,the_Tumr_W,the_Tumr_H, ...
            the_beam_X,the_beam_Y,the_beam_Z);

        save  (['Photon_51_TiO2_',num2str(the_TiO2),'_Indik_',num2str(the_indk),'.mat' ], ...
            'p_in','w_in', ...
            'p_ot','w_ot', ...
            'T_mua_TiO2', 'T_mus_TiO2', 'T_g_TiO2', 'T_n_TiO2', ...
            'T_mua_indk', 'T_mus_indk', 'T_g_indk', 'T_n_indk', ...
            'the_TiO2', 'the_indk', 'the_beam_X', 'the_beam_Y', 'the_beam_Z', 'the_Tumr_X', 'the_Tumr_Y', 'the_Tumr_Z', 'the_Tumr_L', 'the_Tumr_W', 'the_Tumr_H', ...
            's','dl','no_of_photons')

        clearvars p_in w_in
        clearvars p_ot w_ot
        clearvars s no_of_photons
        clearvars the_indk  the_TiO2
        clearvars T_mua_indk T_mus_indk T_g_indk T_n_indk
        clearvars T_mua_TiO2 T_mus_TiO2 T_g_TiO2 T_n_TiO2
    end
end
end

function [p_in, w_in, p_ot, w_ot, s, dl, no_of_photons] = do_simulation (dom_size , opt_TiO2, opt_Indk, dl, wvlngth, nPhotonsReq, nExamplePaths, the_Tumr_X, the_Tumr_Y, the_Tumr_Z, the_Tumr_L, the_Tumr_W, the_Tumr_H, the_beam_X, the_beam_Y, the_beam_Z)
%% Geometry definition
model = MCmatlab.model;
model.G.nx                = dom_size(1);  % Number of bins in the x direction
model.G.ny                = dom_size(2);  % Number of bins in the y direction
model.G.nz                = dom_size(3);  % Number of bins in the z direction
model.G.Lx                = dom_size(1)*dl; % [cm] x size of simulation cuboid
model.G.Ly                = dom_size(2)*dl; % [cm] y size of simulation cuboid
model.G.Lz                = dom_size(3)*dl; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc;  % Media properties defined as a function at the end of this file
model.G.mediaPropParams     = {opt_TiO2(1),opt_TiO2(2),opt_TiO2(3),opt_TiO2(4),opt_Indk(1),opt_Indk(2),opt_Indk(3),opt_Indk(4)}; % A cell array that you can use to contain all sorts of inputs you would like to use inside the mediaPropertiesFunc
model.G.geomFunc            = @geometryDefinition;   % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.
model.G.geomFuncParams      = {the_Tumr_X, the_Tumr_Y, the_Tumr_Z, the_Tumr_L, the_Tumr_W, the_Tumr_H}; % A cell array that you can use to contain all sorts of inputs you would like to use inside the geomFunc
% model = plot(model,'G');
clearvars dom_size
clearvars opt_TiO2 opt_Indk
clearvars the_Tumr_X the_Tumr_Y the_Tumr_Z
clearvars the_Tumr_L the_Tumr_W the_Tumr_H

%% Monte Carlo simulation
model.MC.simulationTimeRequested  = 0.25;    % [min] Time duration of the simulation
model.MC.matchedInterfaces        = false;   % Assumes all refractive indices are the same
model.MC.boundaryType             = 1;       % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.nPhotonsRequested = nPhotonsReq;    % The time to run the MC simualtion for, in minutes. The number of photos launched will vary from run to run.
model.MC.nExamplePaths = nExamplePaths;      % The code stores the paths of the first N photons for subsequent visualization during the plotting.
model.MC.wavelength = wvlngth;               % [nm] Excitation wavelength, used for determination of optical properties for excitation light
clearvars wvlngth nPhotonsReq

%% For a pencil beam, the "focus" is just a point that the beam goes through, here set to be the center of the cuboid:
model.MC.lightSource.sourceType = 0; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.silentMode = false;
beam_X = the_beam_X; beam_Y = the_beam_Y; beam_Z = the_beam_Z; beam_tht = 0; beam_phi = 0;
model.MC.lightSource.xFocus       = beam_X;   % [cm] x position of focus
model.MC.lightSource.yFocus       = beam_Y;   % [cm] y position of focus
model.MC.lightSource.zFocus       = beam_Z;   % [cm] z position of focus
model.MC.lightSource.theta        = beam_tht; % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = beam_phi; % [rad] Azimuthal angle of beam center axis
clearvars the_beam_X the_beam_Y the_beam_Z beam_X beam_Y beam_Z beam_phi beam_tht

%% These lines will run the Monte Carlo simulation with the provided parameters and subsequently plot the results:
% figure(17)
scatter_stat = [];
no_of_photons = 0;
while size(scatter_stat,1) < nExamplePaths
    t_model = runMonteCarlo(model);
    M_raw = t_model.G.M_raw;
    % t_model = plot(t_model,'MC');
    % do calc
    spratrs = find(isnan(t_model.MC.examplePaths(1,:)));
    strt_pnts = spratrs(2:2:end); strt_pnts = strt_pnts(1:end-1);
    fnsh_pnts = spratrs(1:2:end); fnsh_pnts = fnsh_pnts(2:end-0);
    t_scatter_stat = nan(length(fnsh_pnts),3+1+3+1+1); clearvars spratrs
    for idx = 1:length(strt_pnts)
        Photon_Path = t_model.MC.examplePaths(:,strt_pnts(idx)+1:fnsh_pnts(idx)-1);
        if Photon_Path(3,end) <= 0.0
            % plot3(Photon_Path(1,:),Photon_Path(2,:),-Photon_Path(3,:),'r'), hold on
            t_scatter_stat(idx,:) = [...
                Photon_Path(1,1), ...
                Photon_Path(2,1), ...
                Photon_Path(3,1), ...
                Photon_Path(4,1), ...
                Photon_Path(1,end), ...
                Photon_Path(2,end), ...
                Photon_Path(3,end), ...
                Photon_Path(4,end), ...
                sum(sqrt(sum((Photon_Path(1:3,2:end-0)-Photon_Path(1:3,1:end-1)).^2,1))), ...
                ];
        end
        clearvars Photon_Path
    end
    clearvars strt_pnts fnsh_pnts idx t_model

    idx_photons = ~isnan(t_scatter_stat(:,1));
    t_scatter_stat = t_scatter_stat(idx_photons,:); clearvars idx_photons

    no_of_photons = no_of_photons + size(t_scatter_stat,1);
    scatter_stat = [scatter_stat ; t_scatter_stat];
    disp('****************************************************')
    disp(['progress ... ',num2str(size(scatter_stat,1)/nExamplePaths*100),'%'])
    disp('****************************************************')
    clearvars surf_photons t_scatter_stat
end
% figure(17), hold off, axis equal, view([0 0]), grid on
% axis([min(model.G.x) max(model.G.x) min(model.G.y) max(model.G.y) -max(model.G.z) -min(model.G.z)])
% xlabel('x (cm)'), ylabel('y (cm)'), zlabel('z (cm)'), set(gca,'fontsize',24)

clearvars model
% histcount
p_in = scatter_stat(:,1:3);                                                % source positions
w_in = scatter_stat(:,4);                                                  % source weights
p_ot = scatter_stat(:,5:7);                                                % detedtor positions
w_ot = scatter_stat(:,8);                                                  % detedtor weights
s    = scatter_stat(:,9);                                                  % true photon pathlength
end
function M = geometryDefinition (X,Y,Z,parameters)
% Geometry function(s) (see readme for details)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% mediaPropertiesFunc) at each voxel location.
M = 1.*ones(size(X)); % all Air
M = 2.*ones(size(X)); % all TiO2 in Agar
% Logical mask for cylinder aligned along X axis
mask = ( ((Y-parameters{2})./parameters{5}).^2 + ((Z-parameters{3})./parameters{6}).^2 <= 1 ) & ...
    (X >= (parameters{1}-(parameters{4}/2))) & ...
    (X <= (parameters{1}+(parameters{4}/2)));
% Set values inside the cylinder to 3
M(mask) = 3;
M = uint8(M);
end
function mediaProperties = mediaPropertiesFunc (var)
mediaProperties = MCmatlab.mediumProperties;
% Put in your own media property definitions below at 800 nm
j = 1;
mediaProperties(j).name  = 'air';
mediaProperties(j).mua   = 1e-4; % Absorption coefficient [cm^-1]
mediaProperties(j).mus   = 1e-5; % Scattering coefficient [cm^-1]
mediaProperties(j).g     = 1;    % Henyey-Greenstein scattering anisotropy
mediaProperties(j).n     = 1;    % The refractive index

j = 2;
mediaProperties(j).name  = 'TiO2 in Agar';
mediaProperties(j).mua   = double(var{1}); % Absorption coefficient [cm^-1]
mediaProperties(j).mus   = double(var{2}); % Scattering coefficient [cm^-1]
mediaProperties(j).g     = double(var{3}); % Henyey-Greenstein scattering anisotropy
mediaProperties(j).n     = double(var{4}); % The refractive index

j = 3;
mediaProperties(j).name  = 'India Ink';
mediaProperties(j).mua   = double(var{5}); % Absorption coefficient [cm^-1]
mediaProperties(j).mus   = double(var{6}); % Scattering coefficient [cm^-1]
mediaProperties(j).g     = double(var{7}); % Henyey-Greenstein scattering anisotropy
mediaProperties(j).n     = double(var{8}); % The refractive index
end
