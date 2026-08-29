function [] = Photon_59 ()
% Repository group: 07_dpf_parameter_studies
% Version role: primary simulation or study entry point
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: w & s vs. d over mu_s & mu_a: similar to 58 but with a non-laser source

clc
close all

z_air = 0.0; % the thickness of air layer
nPhotonsReq = 1e5; % number of requested photon
nExamplePaths = 1000000; % number of example photon paths
% Optical & size properties
Lx = 29.1; Ly = 29.1; Lz = 06.0; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
beam_phi = 0;     % polar angle of beam in cm
beam_tht = 0;     % azimuthal angle of beam in cm
beam_X = 0.0;     % X deviation of beam in cm
beam_Y = 0.0;     % Y deviation of beam in cm

n = 1.33; g = 0.93;
set_of_mua = [0.0275];
set_of_mus = [35/4,35,35*4];

for i_a = 1:length(set_of_mua)
    for i_s = 1:length(set_of_mus)
        mua = set_of_mua(i_a);
        mus = set_of_mus(i_s);
        the_filename = ['Photon_59_mua_',sprintf('%.4f',set_of_mua(i_a)),'_mus_',sprintf('%.2f',set_of_mus(i_s)),'.mat' ];
        if ~exist(the_filename,'file')
            [x, y, z, d, s, w, c, no_of_photons, M_raw, p_ot, p_in] = do_simulation ( ...
                [mua,mus,g   ,n  ], ...
                [nan,nan,nan ,nan], ...
                [Lx,Ly,Lz],[nan,nan,nan],[nan,nan,nan],z_air,001,nPhotonsReq,nExamplePaths,beam_X,beam_Y,beam_phi,beam_tht);

            save  (the_filename, ...
                'x', 'y', 'z', 'd', 's', 'w', 'c', 'no_of_photons', 'M_raw', 'p_in', 'p_ot')

            clearvars x y z d s w c no_of_photons M_raw p_in p_ot
            clearvars mua mus
        end
        clearvars the_filename
    end
end
end

function [x, y, z, d, s, w, c, no_of_photons, M_raw, p_ot, p_in] = do_simulation (opt_bckg,opt_sgnl,cmd_size,sig_pos,sig_size,zSurface,wvlngth,nPhotonsReq,nExamplePaths,beam_X,beam_Y,beam_phi,beam_tht)
%% Geometry definition
model = MCmatlab.model;
dl = 0.1; % spatial resolution (cm)
model.G.nx                = round(cmd_size(1)./dl);  % Number of bins in the x direction
model.G.ny                = round(cmd_size(2)./dl);  % Number of bins in the y direction
model.G.nz                = round(cmd_size(3)./dl);  % Number of bins in the z direction
model.G.Lx                = cmd_size(1); % [cm] x size of simulation cuboid
model.G.Ly                = cmd_size(2); % [cm] y size of simulation cuboid
model.G.Lz                = cmd_size(3); % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc;  % Media properties defined as a function at the end of this file
model.G.mediaPropParams     = {opt_bckg(1),opt_bckg(2),opt_bckg(3),opt_bckg(4),opt_sgnl(1),opt_sgnl(2),opt_sgnl(3),opt_sgnl(4)}; % A cell array that you can use to contain all sorts of inputs you would like to use inside the mediaPropertiesFunc
model.G.geomFunc            = @geometryDefinition;   % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.
model.G.geomFuncParams      = {zSurface,2,3,sig_pos(1),sig_pos(2),sig_pos(3),sig_size(1),sig_size(2),sig_size(3)}; % A cell array that you can use to contain all sorts of inputs you would like to use inside the geomFunc
% model = plot(model,'G');
clearvars opt_bckg opt_sgnl cmd_size sig_pos sig_size dl

%% Monte Carlo simulation
model.MC.simulationTimeRequested  = 0.25;    % [min] Time duration of the simulation
model.MC.matchedInterfaces        = false;   % Assumes all refractive indices are the same
model.MC.boundaryType             = 1;       % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.nPhotonsRequested = nPhotonsReq;    % The time to run the MC simualtion for, in minutes. The number of photos launched will vary from run to run.
model.MC.nExamplePaths = nPhotonsReq;        % The code stores the paths of the first N photons for subsequent visualization during the plotting.
model.MC.wavelength = wvlngth;               % [nm] Excitation wavelength, used for determination of optical properties for excitation light
clearvars wvlngth nPhotonsReq

%% For a pencil beam, the "focus" is just a point that the beam goes through, here set to be the center of the cuboid:
% LED750L: narrow divergence due to ball lens (viewing half-angle ~11°)
model.MC.lightSource.sourceType = 4;   % radial-factorizable beam
% Spatial: top-hat disk (use the ball-lens clear aperture as a reasonable approximation)
% Datasheet drawing shows ~4.7 mm lens diameter => radius = 2.35 mm = 0.235 cm  :contentReference[oaicite:1]{index=1}
model.MC.lightSource.focalPlaneIntensityDistribution.radialDistr = 0;      % 0 = Top-hat
model.MC.lightSource.focalPlaneIntensityDistribution.radialWidth = 0.235;  % [cm] radius
% Angular: top-hat cone using viewing half-angle (narrow)
model.MC.lightSource.angularIntensityDistribution.radialDistr = 0;         % 0 = Top-hat angular
model.MC.lightSource.angularIntensityDistribution.radialWidth = deg2rad(11); % [rad] viewing half-angle :contentReference[oaicite:2]{index=2}

model.MC.silentMode = false;
model.MC.lightSource.xFocus       = beam_X;   % [cm] x position of focus
model.MC.lightSource.yFocus       = beam_Y;   % [cm] y position of focus
model.MC.lightSource.zFocus       = zSurface; % [cm] z position of focus
model.MC.lightSource.theta        = beam_tht; % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = beam_phi; % [rad] Azimuthal angle of beam center axis
clearvars zSurface beam_phi beam_tht

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
    t_scatter_stat = nan(length(fnsh_pnts),3+1+1+3+1); clearvars spratrs
    for idx = 1:length(strt_pnts)
        Photon_Path = t_model.MC.examplePaths(:,strt_pnts(idx)+1:fnsh_pnts(idx)-1);
        if Photon_Path(3,end) <= 0.0 % diffuse reflectance
            % plot3(Photon_Path(1,:),Photon_Path(2,:),-Photon_Path(3,:),'r'), hold on
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
        elseif Photon_Path(3,end) >= t_model.G.Lz % transmittance
            % plot3(Photon_Path(1,:),Photon_Path(2,:),-Photon_Path(3,:),'g'), hold on
            t_scatter_stat(idx,:) = [...
                Photon_Path(1,end), ...
                Photon_Path(2,end), ...
                Photon_Path(3,end), ...
                sum(sqrt(sum((Photon_Path(1:3,2:end-0)-Photon_Path(1:3,1:end-1)).^2,1))), ...
                Photon_Path(4,end), ...
                (Photon_Path(1,1)), ...
                (Photon_Path(2,1)), ...
                (Photon_Path(3,1)), ...
                (1), ...
                ];
        else % absorbance
            % plot3(Photon_Path(1,:),Photon_Path(2,:),-Photon_Path(3,:),'k'), hold on
        end
        clearvars Photon_Path
    end
    clearvars strt_pnts fnsh_pnts idx t_model

    surf_photons = ~isnan(t_scatter_stat(:,1));
    no_of_photons = no_of_photons + size(t_scatter_stat,1);
    disp('****************************************************')
    disp(['progress ... ',num2str(size(scatter_stat,1)/nExamplePaths*100),'%'])
    disp('****************************************************')
    t_scatter_stat = t_scatter_stat(surf_photons,:);

    scatter_stat = [scatter_stat ; t_scatter_stat];
    clearvars surf_photons t_scatter_stat
end
% figure(17), hold off, axis equal, view([-90 0]), grid on
% axis([min(model.G.x) max(model.G.x) min(model.G.y) max(model.G.y) -max(model.G.z) -min(model.G.z)])
% xlabel('x (cm)'), ylabel('y (cm)'), zlabel('z (cm)'), set(gca,'fontsize',24)
% figure(17), hold off, axis equal, view([-0 +0]), grid on
% axis([-7.5 +7.5 -7.5 +7.5 -max(model.G.z) -min(model.G.z)])
% xlabel('x (cm)'), ylabel('y (cm)'), zlabel('z (cm)'), set(gca,'fontsize',24)

% histcount
p_ot = scatter_stat(:,1:3);                                            % detedtor positions
p_in = scatter_stat(:,6:8);                                            % source positions
x = scatter_stat(:,1);                                                 % x -> source to detector distance
y = scatter_stat(:,2);                                                 % y -> source to detector distance
z = scatter_stat(:,3);                                                 % y -> source to detector distance
d = sqrt((p_in(:,1)-p_ot(:,1)).^2+(p_in(:,2)-p_ot(:,2)).^2);           % d -> source to detector distance
s = scatter_stat(:,4);                                                 % s -> true distance
w = scatter_stat(:,5);                                                 % w -> intensity
c = scatter_stat(:,9);                                                 % c -> code: diffuse reflectance, 0; transmittive, 1
end
function M = geometryDefinition (X,Y,Z,parameters)
% Geometry function(s) (see readme for details)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% mediaPropertiesFunc) at each voxel location.
M = ones(size(X)); % all Air
zSurface       = double(parameters{1}); % Air thickness
M(Z>=zSurface) = double(parameters{2}); % standard (breast) tissue
M(sqrt( ...
    ((X-double(parameters{4})         ).^2)./(double(parameters{7}).^2) + ...
    ((Y-double(parameters{5})         ).^2)./(double(parameters{8}).^2) + ...
    ((Z-double(parameters{6})-zSurface).^2)./(double(parameters{9}).^2))  ...
    <=1) = double(parameters{3}); % non-standard (tumor) tissue
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
mediaProperties(j).name  = 'breast tissue';
mediaProperties(j).mua   = double(var{1}); % Absorption coefficient [cm^-1]
mediaProperties(j).mus   = double(var{2}); % Scattering coefficient [cm^-1]
mediaProperties(j).g     = double(var{3}); % Henyey-Greenstein scattering anisotropy
mediaProperties(j).n     = double(var{4}); % The refractive index

j = 3;
mediaProperties(j).name  = 'tumor tissue';
mediaProperties(j).mua   = double(var{5}); % Absorption coefficient [cm^-1]
mediaProperties(j).mus   = double(var{6}); % Scattering coefficient [cm^-1]
mediaProperties(j).g     = double(var{7}); % Henyey-Greenstein scattering anisotropy
mediaProperties(j).n     = double(var{8}); % The refractive index
end
