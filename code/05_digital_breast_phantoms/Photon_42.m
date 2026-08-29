function [] = Photon_42 ()
% Repository group: 05_digital_breast_phantoms
% Version role: primary simulation or study entry point
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: simulating phantoms prepared by Diego

clc
close all

wvlnt = 800; % sample wavelength (nm)
nPhotonsReq = 1e5;  % number of requested photon
nExamplePaths = 1000; % number of example photon paths
% Optical properties of breast and tumor [1-4]
mua_glnd = 0.55; mus_glnd = 332.7; g_glnd = 0.965; n_glnd = 1.4;
mua_adps = 0.82; mus_adps = 313.8; g_adps = 0.976; n_adps = 1.4;
set_of_percnt = [25,50,75]; % percent of adipose
set_of_versns = [0,1,2]; % version of each density
Nz = 1020; Ny = 257; Nx = 323; dl = 0.2/10; r_lowres = 0.5; dl = dl./r_lowres;

for i_percent = 1:length(set_of_percnt)
    for i_versns = 1:length(set_of_versns)
        the_version = set_of_versns(i_versns);
        the_percent = set_of_percnt(i_percent);
        the_path = ['/home/kaiser/Phantoms/Diego/Phantoms/CupB_Pd',num2str(the_percent),'_',num2str(the_version),'_y_-z_x.obj'];
        fid = fopen(the_path, 'r');
        data = fread(fid,'float32');
        fclose(fid);
        data(        data<=0.0) = 0.0; % Air
        data(0.0<data&data<0.9) = 0.1; % Fat (Adipose Tissue)
        data(0.9<=data        ) = 0.2; % Fibroglandular Tissue
        data = round(data*10);
        % Reshape the data into a 3D array
        TheImage = reshape(data, [Nz, Ny, Nx]);
        TheImage = round(imresize3(TheImage,r_lowres));
        TheImage = permute(TheImage,[1,3,2]); % TheImage = TheImage(:,:,end:-1:1);
        clearvars data the_path

        [p_in, i_in, m_in, w_in, p_ot, i_ot, m_ot, w_ot, s, dl, no_of_photons, M_raw] = do_simulation ( ...
            [mua_glnd,mus_glnd,g_glnd,n_glnd], ...
            [mua_adps,mus_adps,g_adps,n_adps], ...
            TheImage,dl,wvlnt,nPhotonsReq,nExamplePaths);

        save  (['Photon_42_Pd',num2str(the_percent),'_',num2str(the_version),'_y_-z_x.mat' ], ...
            'p_in', 'i_in', 'm_in', 'w_in', ...
            'p_ot', 'i_ot', 'm_ot', 'w_ot', ...
            's','dl','no_of_photons','M_raw')

        clearvars p_in i_in m_in w_in
        clearvars p_ot i_ot m_ot w_ot
        clearvars s no_of_photons M_raw
        clearvars fid the_percent the_version TheImage
    end
end
end

function [p_in, i_in, m_in, w_in, p_ot, i_ot, m_ot, w_ot, s, dl, no_of_photons, M_raw] = do_simulation (opt_glnd, opt_adps, M_in, dl, wvlngth, nPhotonsReq, nExamplePaths)
%% Geometry definition
model = MCmatlab.model;
model.G.nx                = size(M_in,1);  % Number of bins in the x direction
model.G.ny                = size(M_in,2);  % Number of bins in the y direction
model.G.nz                = size(M_in,3);  % Number of bins in the z direction
model.G.Lx                = size(M_in,1)*dl; % [cm] x size of simulation cuboid
model.G.Ly                = size(M_in,2)*dl; % [cm] y size of simulation cuboid
model.G.Lz                = size(M_in,3)*dl; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc;  % Media properties defined as a function at the end of this file
model.G.mediaPropParams     = {opt_adps(1),opt_adps(2),opt_adps(3),opt_adps(4),opt_glnd(1),opt_glnd(2),opt_glnd(3),opt_glnd(4)}; % A cell array that you can use to contain all sorts of inputs you would like to use inside the mediaPropertiesFunc
model.G.geomFunc            = @geometryDefinition;   % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.
model.G.geomFuncParams      = {M_in}; % A cell array that you can use to contain all sorts of inputs you would like to use inside the geomFunc
model = plot(model,'G');
clearvars opt_glnd opt_adps M_in

%% Monte Carlo simulation
model.MC.simulationTimeRequested  = 0.25;    % [min] Time duration of the simulation
model.MC.matchedInterfaces        = false;   % Assumes all refractive indices are the same
model.MC.boundaryType             = 1;       % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.nPhotonsRequested = nPhotonsReq;    % The time to run the MC simualtion for, in minutes. The number of photos launched will vary from run to run.
model.MC.nExamplePaths = nExamplePaths;      % The code stores the paths of the first N photons for subsequent visualization during the plotting.
model.MC.wavelength = wvlngth;               % [nm] Excitation wavelength, used for determination of optical properties for excitation light
clearvars wvlngth nPhotonsReq

%% For a pencil beam, the "focus" is just a point that the beam goes through, here set to be the center of the cuboid:
model.MC.lightSource.sourceType   = 0;    % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.silentMode = false;
beam_X = 0; beam_Y = 0; beam_Z = 0; beam_tht = 0; beam_phi = 0;
[ix_beam,iy_beam,iz_beam] = findIndex(beam_X,beam_Y,beam_Z,model.G.x,model.G.y,model.G.z);
while model.G.M_raw(ix_beam,iy_beam,iz_beam)==1
    beam_Z = iz_beam*dl-dl;
    iz_beam = iz_beam+1;
end
disp([ix_beam,iy_beam,iz_beam-1])
model.MC.lightSource.xFocus       = beam_X;   % [cm] x position of focus
model.MC.lightSource.yFocus       = beam_Y;   % [cm] y position of focus
model.MC.lightSource.zFocus       = beam_Z;   % [cm] z position of focus
model.MC.lightSource.theta        = beam_tht; % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = beam_phi; % [rad] Azimuthal angle of beam center axis
clearvars beam_X beam_Y beam_Z beam_phi beam_tht ix_beam iy_beam iz_beam

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
    t_scatter_stat = nan(length(fnsh_pnts),3+3+2+3+3+2+1); clearvars spratrs
    for idx = 1:length(strt_pnts)
        Photon_Path = t_model.MC.examplePaths(:,strt_pnts(idx)+1:fnsh_pnts(idx)-1);

        [ix_e,iy_e,iz_e] = findIndex(Photon_Path(1,end),  Photon_Path(2,end),  Photon_Path(3,end),  t_model.G.x,t_model.G.y,t_model.G.z); M_e = double(t_model.G.M_raw(ix_e,iy_e,iz_e)); clearvars ix_e iy_e iz_e
        [ix_1,iy_1,iz_1] = findIndex(Photon_Path(1,end-1),Photon_Path(2,end-1),Photon_Path(3,end-1),t_model.G.x,t_model.G.y,t_model.G.z); M_1 = double(t_model.G.M_raw(ix_1,iy_1,iz_1)); clearvars ix_1 iy_1 iz_1
        while (M_e==1&&M_1==1&&size(Photon_Path,2)>2)
            clearvars M_e M_1
            Photon_Path = Photon_Path(:,1:end-1);
            [ix_e,iy_e,iz_e] = findIndex(Photon_Path(1,end),  Photon_Path(2,end),  Photon_Path(3,end),  t_model.G.x,t_model.G.y,t_model.G.z); M_e = double(t_model.G.M_raw(ix_e,iy_e,iz_e)); clearvars ix_e iy_e iz_e
            [ix_1,iy_1,iz_1] = findIndex(Photon_Path(1,end-1),Photon_Path(2,end-1),Photon_Path(3,end-1),t_model.G.x,t_model.G.y,t_model.G.z); M_1 = double(t_model.G.M_raw(ix_1,iy_1,iz_1)); clearvars ix_1 iy_1 iz_1
        end
        clearvars M_e M_1

        [ix_f,iy_f,iz_f] = findIndex(Photon_Path(1,1),  Photon_Path(2,1),  Photon_Path(3,1),  t_model.G.x,t_model.G.y,t_model.G.z); M_f = double(t_model.G.M_raw(ix_f,iy_f,iz_f));
        [ix_e,iy_e,iz_e] = findIndex(Photon_Path(1,end),Photon_Path(2,end),Photon_Path(3,end),t_model.G.x,t_model.G.y,t_model.G.z); M_e = double(t_model.G.M_raw(ix_e,iy_e,iz_e));
        if size(Photon_Path,2)>2
            t_scatter_stat(idx,:) = [...
                Photon_Path(1,1), ...
                Photon_Path(2,1), ...
                Photon_Path(3,1), ...
                ix_f, ...
                iy_f, ...
                iz_f, ...
                (M_f), ...
                Photon_Path(4,1), ...
                Photon_Path(1,end), ...
                Photon_Path(2,end), ...
                Photon_Path(3,end), ...
                ix_e, ...
                iy_e, ...
                iz_e, ...
                (M_e), ...
                Photon_Path(4,end), ...
                sum(sqrt(sum((Photon_Path(1:3,2:end-0)-Photon_Path(1:3,1:end-1)).^2,1))), ...
                ];
        end
        clearvars Photon_Path M_f M_e ix_f iy_f iz_f ix_e iy_e iz_e
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

clearvars model
% histcount
p_in = scatter_stat(:,1:3);                                                % source positions
i_in = scatter_stat(:,4:6);                                                % source indices
m_in = scatter_stat(:,7);                                                  % source media
w_in = scatter_stat(:,8);                                                  % source weights
p_ot = scatter_stat(:, 9:11);                                              % detedtor positions
i_ot = scatter_stat(:,12:14);                                              % detedtor indices
m_ot = scatter_stat(:,15);                                                 % detedtor media
w_ot = scatter_stat(:,16);                                                 % detedtor weights
s    = scatter_stat(:,17);                                                 % true photon pathlength
end
function M = geometryDefinition (X,Y,Z,parameters)
% Geometry function(s) (see readme for details)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% mediaPropertiesFunc) at each voxel location.
M = parameters{1} + 1;
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
mediaProperties(j).name  = 'adipose tissue';
mediaProperties(j).mua   = double(var{1}); % Absorption coefficient [cm^-1]
mediaProperties(j).mus   = double(var{2}); % Scattering coefficient [cm^-1]
mediaProperties(j).g     = double(var{3}); % Henyey-Greenstein scattering anisotropy
mediaProperties(j).n     = double(var{4}); % The refractive index

j = 3;
mediaProperties(j).name  = 'fibroglandular tissue';
mediaProperties(j).mua   = double(var{5}); % Absorption coefficient [cm^-1]
mediaProperties(j).mus   = double(var{6}); % Scattering coefficient [cm^-1]
mediaProperties(j).g     = double(var{7}); % Henyey-Greenstein scattering anisotropy
mediaProperties(j).n     = double(var{8}); % The refractive index
end
function [out] = get_color(idx)
if     idx==1
    out = [0.0000 0.4470 0.7410];
elseif idx==2
    out = [0.8500 0.3250 0.0980];
elseif idx==3
    out = [0.9290 0.6940 0.1250];
elseif idx==4
    out = [0.4940 0.1840 0.5560];
elseif idx==5
    out = [0.4660 0.6740 0.1880];
elseif idx==6
    out = [0.3010 0.7450 0.9330];
elseif idx==7
    out = [0.6350 0.0780 0.1840];
else
    out = [0.0000 0.0000 0.0000];
end
end
function [ix,iy,iz] = findIndex(x,y,z,X,Y,Z)
[~,ix] = min(abs(X-x));
[~,iy] = min(abs(Y-y));
[~,iz] = min(abs(Z-z));
end
