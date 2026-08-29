function [] = Photon_62 ()
% Repository group: 07_dpf_parameter_studies
% Version role: primary simulation or study entry point
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: studies the effect of a gradient versus a uniform dye distribution

clc
close all

% 1D tartrazine diffusion into chicken (one-side boundary, semi-infinite approx)
% Model: C(z,t) = C0 * erfc( z / (2*sqrt(D*t)) )
% Assumptions:
%   - Chicken initially has no tartrazine: C(z,0)=0 for z>0
%   - Surface at z=0 is held at constant concentration C0 (large reservoir, always wet)
%   - One-dimensional diffusion along depth z (no lateral gradients)
%   - Tissue is effectively semi-infinite for the time window of interest
% Notes:
%   - erfc() is the complementary error function.
%   - At t=0, the formula is singular; we handle t=0 explicitly.

% Geometry (cm)
z_air = 0.0;                 % Thickness of an air layer above tissue [cm] (not used in 1D diffusion eqn)
z_bot = 3.0;                 % Dye penetration depth [cm]
Lx = 19.1; Ly = 19.1;        % Lateral domain sizes [cm] (not used in 1D diffusion eqn)
Lz = 6.0;                    % Tissue thickness in z [cm] (only used to set max depth for plotting)

% % Diffusion / boundary condition
% D   = 6.6e-7;                % Diffusion coefficient in chicken breast [cm^2/s]
% C0  = 0.6;                   % Surface concentration at z=0 (reservoir) [M]
%
% % Time and depth grids
% t_h = 0:60*60:30*24*3500;                % Time vector [s] (0 to 1 hour in 1-second steps)
% z_h = 0:0.1:Lz;             % Depth vector [cm] (0 to Lz in 0.01 cm = 0.1 mm steps)
%
% % Concentration field C(z,t)
% [Z, T] = ndgrid(z_h, t_h);
% C = zeros(size(Z));          % Preallocate [M]
% % Avoid division-by-zero at t=0; for t>0 use erfc solution
% mask = (T > 0);
% C(mask) = C0 .* erfc( Z(mask) ./ (2 .* sqrt(D .* T(mask))) );
%
% % Enforce boundary/initial conditions explicitly (optional clarity)
% C(z_h==0, :) = C0;            % Surface held at C0 for all t
% C(:, t_h==0) = 0;             % Initially zero inside (including z=0 gets overwritten above)
% pcolor(T./3600/24,Z,C), shading interp, xlabel('time (days)'), ylabel('Depth z (cm)'), h = colorbar; ylabel(h,'C (M)'), colormap jet, set(gca,'fontsize',18), clim([0 C0])



% Optical properties
nPhotonsReq = 1e5; % number of requested photon
nExamplePaths = 1000000; % number of example photon paths
beam_phi = 0;     % polar angle of beam in cm
beam_tht = 0;     % azimuthal angle of beam in cm
beam_X = 0.0;     % X deviation of beam in cm
beam_Y = 0.0;     % Y deviation of beam in cm

g = 0.9;
no_of_layer = 7;
z_edges = linspace(z_air,z_bot,no_of_layer+1);
n_centr = linspace(1.44,1.35,no_of_layer);
mus_lyr = 300.*((1.46./n_centr - 1)./(1.46./1.35 - 1)+eps).^2; % 1/cm
mua = 0.3341; % 1/cm
% Dye distribution modes:
% 0   → no dye (high scattering)
% 1/2 → layered (gradient) distribution (low → high scattering with depth)
% 1   → uniform distribution (low scattering)
% -1  → randomly layered distribution (layers randomly distributed)
set_of_status = [0, 1/2, 1, -1];

for i_stat = 1:length(set_of_status)
    the_status = set_of_status(i_stat);
    the_filename = ['Photon_62_status_',sprintf('%.2f',set_of_status(i_stat)),'.mat'];
    if ~exist(the_filename,'file')
        [x_in, y_in, z_in, x_ot, y_ot, z_ot, s, w, no_of_photons, M_raw] = do_simulation ( ...
            mua,mus_lyr,g,n_centr,z_edges, ...
            [Lx,Ly,Lz],the_status,001,nPhotonsReq,nExamplePaths,beam_X,beam_Y,beam_phi,beam_tht);
        save  (the_filename, ...
            'x_in', 'y_in', 'z_in', 'x_ot', 'y_ot', 'z_ot','s', 'w', 'no_of_photons', 'M_raw', 'mua', 'mus_lyr', 'n_centr', 'g', 'z_edges')
        clearvars x_in y_in z_in x_ot y_ot z_ot s w no_of_photons M_raw
    end
    clearvars the_filename the_status
    clearvars the_filename
end
end

function [x_in, y_in, z_in, x_ot, y_ot, z_ot, s, w, no_of_photons, M_raw] = do_simulation (mua,mus_lyr,g,n_lyr,z_edges,cmd_size,the_status,wvlngth,nPhotonsReq,nExamplePaths,beam_X,beam_Y,beam_phi,beam_tht)
%% Geometry definition
model = MCmatlab.model;
dl = 0.05;                                           % spatial resolution (cm)
model.G.nx                = round(cmd_size(1)./dl);  % Number of bins in the x direction
model.G.ny                = round(cmd_size(2)./dl);  % Number of bins in the y direction
model.G.nz                = round(cmd_size(3)./dl);  % Number of bins in the z direction
model.G.Lx                = cmd_size(1); % [cm] x size of simulation cuboid
model.G.Ly                = cmd_size(2); % [cm] y size of simulation cuboid
model.G.Lz                = cmd_size(3); % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc;  % Media properties defined as a function at the end of this file
model.G.mediaPropParams     = {mua,mus_lyr,g,n_lyr,z_edges,the_status}; % A cell array that you can use to contain all sorts of inputs you would like to use inside the mediaPropertiesFunc
model.G.geomFunc            = @geometryDefinition;   % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.
model.G.geomFuncParams      = {z_edges};    % A cell array that you can use to contain all sorts of inputs you would like to use inside the geomFunc
% model = plot(model,'G');
clearvars cmd_size dl mua mus_lyr g n_lyr z_edges the_status z_edges zSurface

%% Monte Carlo simulation
model.MC.simulationTimeRequested  = 0.25;    % [min] Time duration of the simulation
model.MC.matchedInterfaces        = false;   % Assumes all refractive indices are the same
model.MC.boundaryType             = 1;       % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.nPhotonsRequested = nPhotonsReq;    % The time to run the MC simualtion for, in minutes. The number of photos launched will vary from run to run.
model.MC.nExamplePaths = nPhotonsReq;        % The code stores the paths of the first N photons for subsequent visualization during the plotting.
model.MC.wavelength = wvlngth;               % [nm] Excitation wavelength, used for determination of optical properties for excitation light
clearvars wvlngth nPhotonsReq

%% For a pencil beam, the "focus" is just a point that the beam goes through, here set to be the center of the cuboid:
model.MC.lightSource.sourceType   = 0;    % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.silentMode = false;
model.MC.lightSource.xFocus       = beam_X;   % [cm] x position of focus
model.MC.lightSource.yFocus       = beam_Y;   % [cm] y position of focus
model.MC.lightSource.zFocus       = 0;        % [cm] z position of focus
model.MC.lightSource.theta        = beam_tht; % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = beam_phi; % [rad] Azimuthal angle of beam center axis
clearvars zSurface beam_phi beam_tht beam_X beam_Y

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
    t_scatter_stat = nan(length(fnsh_pnts),3+3+1+1); clearvars spratrs
    for idx = 1:length(strt_pnts)
        Photon_Path = t_model.MC.examplePaths(:,strt_pnts(idx)+1:fnsh_pnts(idx)-1);
        % plot3(Photon_Path(1,:),Photon_Path(2,:),-Photon_Path(3,:)), hold on
        t_scatter_stat(idx,:) = [...
            (Photon_Path(1,1)), ...
            (Photon_Path(2,1)), ...
            (Photon_Path(3,1)), ...
            Photon_Path(1,end), ...
            Photon_Path(2,end), ...
            Photon_Path(3,end), ...
            sum(sqrt(sum((Photon_Path(1:3,2:end-0)-Photon_Path(1:3,1:end-1)).^2,1))), ...
            Photon_Path(4,end), ...
            ];
        clearvars Photon_Path
    end
    clearvars strt_pnts fnsh_pnts idx t_model

    no_of_photons = no_of_photons + size(t_scatter_stat,1);
    scatter_stat = [scatter_stat ; t_scatter_stat];
    disp('****************************************************')
    disp(['fluence progress ... ',num2str(size(scatter_stat,1)/nExamplePaths*100),'%'])
    disp('****************************************************')
    clearvars surf_photons t_scatter_stat
end
% figure(17), hold off, axis equal, view([-90 0]), grid on
% axis([min(model.G.x) max(model.G.x) min(model.G.y) max(model.G.y) -max(model.G.z) -min(model.G.z)])
% xlabel('x (cm)'), ylabel('y (cm)'), zlabel('z (cm)'), set(gca,'fontsize',24)
% figure(17), hold off, axis equal, view([-0 +0]), grid on
% axis([-7.5 +7.5 -7.5 +7.5 -max(model.G.z) -min(model.G.z)])
% xlabel('x (cm)'), ylabel('y (cm)'), zlabel('z (cm)'), set(gca,'fontsize',24)

% histcount
x_in = scatter_stat(:,1);                                                 % x -> input
y_in = scatter_stat(:,2);                                                 % y -> input
z_in = scatter_stat(:,3);                                                 % y -> input
x_ot = scatter_stat(:,4);                                                 % x -> output
y_ot = scatter_stat(:,5);                                                 % y -> output
z_ot = scatter_stat(:,6);                                                 % y -> output
s = scatter_stat(:,7);                                                    % s -> true distance
w = scatter_stat(:,8);                                                    % w -> intensity
end

function M = geometryDefinition (X,Y,Z,parameters)
% Geometry function(s) (see readme for details)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% mediaPropertiesFunc) at each voxel location.
M = ones(size(X)); % all Air
z_edges  = double(parameters{1});
M(Z<=z_edges(1)) = 1; % air
for j = 2:length(z_edges)
    M(z_edges(j-1)<=Z&Z<=z_edges(j)) = j; % tissue j+1
end
M(z_edges(j)<=Z) = j;
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

% extract parameters
mua        = double(var{1});
mus_lyr    = double(var{2});
g          = double(var{3});
n_lyr      = double(var{4});
z_edges    = double(var{5});
the_status = double(var{6});

% set order of layers
if     the_status == 0 % no dye
    set_of_mua = mua.*ones(size(n_lyr));
    set_of_mus = mus_lyr(end).*ones(size(n_lyr));
    set_of___n = n_lyr(end).*ones(size(n_lyr));
    set_of___g = g.*ones(size(n_lyr));
elseif the_status == 1/2 % gradient dye
    set_of_mua = mua.*ones(size(n_lyr));
    set_of_mus = mus_lyr;
    set_of___n = n_lyr;
    set_of___g = g.*ones(size(n_lyr));
elseif the_status == 1 % uniform dye
    set_of_mua = mua.*ones(size(n_lyr));
    set_of_mus = mus_lyr(1).*ones(size(n_lyr));
    set_of___n = n_lyr(1).*ones(size(n_lyr));
    set_of___g = g.*ones(size(n_lyr));
else                   % random dye layers
    set_of_mua = mua.*ones(size(n_lyr));
    set_of_mus = mus_lyr(randperm(length(mus_lyr)));
    set_of___n = n_lyr(randperm(length(n_lyr)));
    set_of___g = g.*ones(size(n_lyr));
end
clearvars mua mus_lyr n_lyr g the_status j

for j = 2:length(z_edges)
    mediaProperties(j).name  = ['layer #',num2str(j)];
    mediaProperties(j).mua   = set_of_mua(j-1); % Absorption coefficient [cm^-1]
    mediaProperties(j).mus   = set_of_mus(j-1); % Scattering coefficient [cm^-1]
    mediaProperties(j).g     = set_of___g(j-1); % Henyey-Greenstein scattering anisotropy
    mediaProperties(j).n     = set_of___n(j-1); % The refractive index
end
end
