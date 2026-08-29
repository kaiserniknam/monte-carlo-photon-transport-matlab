function [] = Photon_56_q ()
% Repository group: 08_oxygenation_and_proposal_studies
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor and Signal on DPF/PPF Simulations
% This version: Comparison of source–detector separation (s vs. d) in breast tissue
% Context: Includes no embedded tumor
% Dr. Das’s Grant Proposal - The same as 55 but simple plot
% Notes, simulations, and related analysis for proposal development
% ref: Monte Carlo investigation of the effect of blood volume and oxygen saturation on optical path in reflectance pulse oximetry
% S Chatterjee et. al. 2016 Biomed. Phys. Eng. Express 2 065018

clc
close all

z_air = 0.0; % the thickness of air layer
nPhotonsReq = 1e5;  % number of requested photon
nExamplePaths = 500000; % number of example photon paths
% Optical properties of breast and tumor [1-4]
Lx = 19.1; Ly = 19.1; Lz = 10.0; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
dx = 2.00; cy = 0.00;            % Tumor depth from the skin (in cm)
rx = 0.65; ry = 0.65; rz = 0.65; % Tumor radius within the breast (in cm)
beam_phi = 0;     % polar angle of beam in cm
beam_tht = 0;     % azimuthal angle of beam in cm

set_of_beam_X = -2:2:+2;               % X deviation of beam in cm
set_of_beam_Y = -2:2:+2;               % Y deviation of beam in cm
set_of_czs = [1, 2];                   % tumor depth in cm
set_of_wvl = [660 940];                % wavelength

for i_BeamX = 1:length(set_of_beam_X)
    for i_BeamY = 1:length(set_of_beam_Y)
        for i_czs = 1:length(set_of_czs)
            for i_wvl = 1:length(set_of_wvl)
                beam_X = set_of_beam_X(i_BeamX);
                beam_Y = set_of_beam_Y(i_BeamY);
                cz = set_of_czs(i_czs);
                wvlngth = set_of_wvl(i_wvl);

                the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_56_null_X_',num2str(beam_X),'_Y_',num2str(beam_Y),'_cz_',num2str(cz),'_wvlngth_',num2str(wvlngth),'.mat' ];
                [opt_brst,opt_HbO2] = give_breast_mua_mus(wvlngth,1);
                [opt_brst,opt_dxHb] = give_breast_mua_mus(wvlngth,0);

                if ~exist(the_filename,'file')
                    [p_in, p_ot, s, d, w, no_of_photons, M_raw] = do_simulation ( ...
                        [opt_brst.mua, opt_brst.mus, opt_brst.g, opt_brst.n], ...
                        [opt_HbO2.mua, opt_HbO2.mus, opt_HbO2.g, opt_HbO2.n], ...
                        [opt_dxHb.mua, opt_dxHb.mus, opt_dxHb.g, opt_dxHb.n], ...
                        [Lx,Ly,Lz],[-dx/2-rx,cy,cz],[+dx/2+rx,cy,cz],[rx,ry,rz],z_air,wvlngth,nPhotonsReq,nExamplePaths,beam_X,beam_Y,beam_phi,beam_tht);
                    save (the_filename, ...
                        'p_in','p_ot','s','d','w','no_of_photons','M_raw','opt_brst','opt_HbO2','opt_dxHb','Lx','Ly','Lz','dx','cy','cz','rx','ry','rz','wvlngth','beam_X','beam_Y','beam_phi','beam_tht')

                    clearvars p_in p_ot s d w no_of_photons M_raw
                else
                    db = load(the_filename);
                    disp(['NUll: Beam X = ',num2str(db.beam_X),', beam Y = ',num2str(db.beam_Y),', wvlngth = ',num2str(db.wvlngth),', dpth = ',num2str(db.cz),': # of photons = ',num2str(db.no_of_photons)])
                end
                clearvars beam_X beam_Y cz wvlngth opt_brst opt_HbO2 opt_dxHb
            end
        end
    end
end
end

function [p_in, p_ot, s, d, w, no_of_photons, M_raw] = do_simulation (opt_bckg,opt_HbO2,opt_dxHb,cmd_size,sgn_HbO2,sgn_dxHb,sig_size,zSurface,wvlngth,nPhotonsReq,nExamplePaths,beam_X,beam_Y,beam_phi,beam_tht)
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
model.G.mediaPropParams     = {opt_bckg(1),opt_bckg(2),opt_bckg(3),opt_bckg(4),opt_HbO2(1),opt_HbO2(2),opt_HbO2(3),opt_HbO2(4),opt_dxHb(1),opt_dxHb(2),opt_dxHb(3),opt_dxHb(4)}; % A cell array that you can use to contain all sorts of inputs you would like to use inside the mediaPropertiesFunc
model.G.geomFunc            = @geometryDefinition;   % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.
model.G.geomFuncParams      = {zSurface,2,2,2,sgn_HbO2(1),sgn_HbO2(2),sgn_HbO2(3),sgn_dxHb(1),sgn_dxHb(2),sgn_dxHb(3),sig_size(1),sig_size(2),sig_size(3)}; % A cell array that you can use to contain all sorts of inputs you would like to use inside the geomFunc
% model = plot(model,'G');
clearvars opt_bckg sgn_HbO2 sgn_dxHb cmd_size sgn_HbO2 sgn_dxHb sig_size dl

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
model.MC.lightSource.zFocus       = zSurface; % [cm] z position of focus
model.MC.lightSource.theta        = beam_tht; % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = beam_phi; % [rad] Azimuthal angle of beam center axis
clearvars zSurface beam_phi beam_tht beam_X beam_Y

%% These lines will run the Monte Carlo simulation with the provided parameters and subsequently plot the results:
% figure(17)
scatter_stat = [];
no_of_photons = 0;
while no_of_photons < nExamplePaths
    t_model = runMonteCarlo(model);
    M_raw = t_model.G.M_raw;
    % t_model = plot(t_model,'MC');
    % do calc
    spratrs = find(isnan(t_model.MC.examplePaths(1,:)));
    strt_pnts = spratrs(2:2:end); strt_pnts = strt_pnts(1:end-1);
    fnsh_pnts = spratrs(1:2:end); fnsh_pnts = fnsh_pnts(2:end-0);
    t_scatter_stat = nan(length(fnsh_pnts),3+1+1+3); clearvars spratrs
    for idx = 1:length(strt_pnts)
        Photon_Path = t_model.MC.examplePaths(:,strt_pnts(idx)+1:fnsh_pnts(idx)-1);
        if Photon_Path(3,end) <= 0.0 % diffuse reflectance
            % figure(17), plot3(Photon_Path(1,:),Photon_Path(2,:),-Photon_Path(3,:),'r'), hold on
            t_scatter_stat(idx,:) = [...
                Photon_Path(1,end), ...
                Photon_Path(2,end), ...
                Photon_Path(3,end), ...
                sum(sqrt(sum((Photon_Path(1:3,2:end-0)-Photon_Path(1:3,1:end-1)).^2,1))), ...
                Photon_Path(4,end), ...
                (Photon_Path(1,1)), ...
                (Photon_Path(2,1)), ...
                (Photon_Path(3,1)), ...
                ];
        elseif Photon_Path(3,end) >= t_model.G.Lz % Transmittance
            % figure(17), plot3(Photon_Path(1,:),Photon_Path(2,:),-Photon_Path(3,:),'g'), hold on
        else % Absorbance
            % figure(17), plot3(Photon_Path(1,:),Photon_Path(2,:),-Photon_Path(3,:),'k'), hold on
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

% histcount
p_ot = scatter_stat(:,1:3);                                                % detedtor positions
s = scatter_stat(:,4);                                                     % photon pathlengths
w = scatter_stat(:,5);                                                     % detedtor weights
p_in = scatter_stat(:,6:8);                                                % source positions
d = sqrt((p_in(:,1)-p_ot(:,1)).^2+(p_in(:,2)-p_ot(:,2)).^2);               % d -> source to detector distance
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
        ((X-double(parameters{ 5})         ).^2)./(double(parameters{11}).^2) + ...
        ((Y-double(parameters{ 6})         ).^2)./(double(parameters{12}).^2) + ...
        ((Z-double(parameters{ 7})-zSurface).^2)./(double(parameters{13}).^2))  ...
        <=1) = double(parameters{3}); % HbO2 tissue
M(sqrt( ...
        ((X-double(parameters{ 8})         ).^2)./(double(parameters{11}).^2) + ...
        ((Y-double(parameters{ 9})         ).^2)./(double(parameters{12}).^2) + ...
        ((Z-double(parameters{10})-zSurface).^2)./(double(parameters{13}).^2))  ...
        <=1) = double(parameters{4}); % HbO2 tissue
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
mediaProperties(j).mua   = double(var{1});  % Absorption coefficient [cm^-1]
mediaProperties(j).mus   = double(var{2});  % Scattering coefficient [cm^-1]
mediaProperties(j).g     = double(var{3});  % Henyey-Greenstein scattering anisotropy
mediaProperties(j).n     = double(var{4});  % The refractive index

j = 3;
mediaProperties(j).name  = 'HbO2 tissue';
mediaProperties(j).mua   = double(var{5});  % Absorption coefficient [cm^-1]
mediaProperties(j).mus   = double(var{6});  % Scattering coefficient [cm^-1]
mediaProperties(j).g     = double(var{7});  % Henyey-Greenstein scattering anisotropy
mediaProperties(j).n     = double(var{8});  % The refractive index

j = 4;
mediaProperties(j).name  = 'dxHb tissue';
mediaProperties(j).mua   = double(var{9});  % Absorption coefficient [cm^-1]
mediaProperties(j).mus   = double(var{10}); % Scattering coefficient [cm^-1]
mediaProperties(j).g     = double(var{11}); % Henyey-Greenstein scattering anisotropy
mediaProperties(j).n     = double(var{12}); % The refractive index
end
function [opt_brst,opt_Hb] = give_breast_mua_mus(lambda_in,StO2)
% loading breast tissue data
brst_db = readtable('DB/breast.normal.csv');
brst_strct(1).lambda = 690; brst_strct(1).mua = 0.030; brst_strct(1).mus = 12.0*20; brst_strct(1).g = 0.95; brst_strct(1).n = 1.4;
brst_strct(2).lambda = 825; brst_strct(2).mua = 0.040; brst_strct(2).mus = 11.0*20; brst_strct(2).g = 0.95; brst_strct(2).n = 1.4;
% mix. of HbO2 & Hb
g   = 0.95;
if     lambda_in==660
    mua = StO2*0.15   + (1-StO2)*1.64  ;
    mup = StO2*1.3844 + (1-StO2)*1.1566;
    mus = mup/(1-g);
    n = 1.41;
elseif lambda_in==940
    mua = StO2*0.65   + (1-StO2)*0.43  ;
    mup = StO2*1.3187 + (1-StO2)*1.1124;
    mus = mup/(1-g);
    n = 1.38;
else
    keyboard
end
opt_Hb.mua = mua*10; % mm^(-1) to cm^(-1)
opt_Hb.mus = mus*10; % mm^(-1) to cm^(-1)
opt_Hb.g = g;
opt_Hb.n = n;
clearvars mua mus mup g n

% breast tissue
[x,iv,~] = unique(brst_db.Var1); v = brst_db.Var2(iv)*10; opt_brst.mua = interp1(x,v,lambda_in,'pchip','extrap'); clearvars x v iv brst_db
bs_brst = -log(brst_strct(2).mus/brst_strct(1).mus)/log(brst_strct(2).lambda/brst_strct(1).lambda); as_brst = brst_strct(1).mus/brst_strct(1).lambda^(-bs_brst); opt_brst.mus = as_brst*(lambda_in)^(-bs_brst); clearvars as_brst bs_brst
opt_brst.g = brst_strct(2).g;
opt_brst.n = brst_strct(1).n;
clearvars brst_strct
end
