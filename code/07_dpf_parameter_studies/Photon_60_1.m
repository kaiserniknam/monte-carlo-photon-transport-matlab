function [] = Photon_60_1 ()
% Repository group: 07_dpf_parameter_studies
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: the same as 33, but for
% K. M. Yoo, Feng Liu, and R. R. Alfano, "Imaging through a scattering wall using absorption," Opt. Lett. 16, 1068-1070 (1991)
% data analysis

clc
close all

c = 3e8;
wvlngth = 620; % wavelength in nm
z_air = 0.0; % the thickness of air layer
% Optical & size properties
Lx = 29.1; Ly = 29.1; Lz = 01.0; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
n = 1.3323; g = 0.9;
set_of_mua = [10/1.83,10/3.4,10/13.5,10/550]; % @ 620 nm
set_of_mus = 10./0.3;
flight_time = linspace(0,900,901);
set_of_hist = nan(length(set_of_mua),length(flight_time)-1);


for i_a = 1:length(set_of_mua)
    for i_s = 1:length(set_of_mus)
        mua = set_of_mua(i_a);
        mus = set_of_mus(i_s);
        % the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_60_mua_',sprintf('%.2f',set_of_mua(i_a)),'_mus_',sprintf('%.2f',set_of_mus(i_s)),'_1e7_g_0.70.mat'];
        the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_60_mua_',sprintf('%.2f',set_of_mua(i_a)),'_mus_',sprintf('%.2f',set_of_mus(i_s)),'_g_0.70.mat'];
        t_db = load(the_filename); clearvars the_filename
        disp(['mua = ',num2str(mua,'%.2f'),', mus = ',num2str(mus,'%.0f')])
        dlta = 0.1;
        set_of_hist(i_a,:) = histcounts(t_db.s(Lz-eps<=t_db.z_ot)./100./c.*n.*1e12,flight_time);
        set_of_hist(i_a,:) = histcounts(t_db.s(Lz-eps<=t_db.z_ot&t_db.x_ot<=dlta&-dlta<=t_db.x_ot&t_db.y_ot<=dlta&-dlta<=t_db.y_ot)./100./c.*n.*1e12,flight_time);
        subplot(4,2,10-i_a*2),
        rectangle('Position',[45 0 35 1],'FaceColor','y','EdgeColor','none'), hold on
        plot(1/2*(flight_time(1:end-1)+flight_time(2:end-0)),set_of_hist(i_a,:)./max(set_of_hist(i_a,:)),'LineWidth',2,'DisplayName',['la = ',num2str(10/mua),' mm']); hold on
        % plot([50 50],[0 1],'r-','HandleVisibility','off','LineWidth',2)
        % plot([75 75],[0 1],'r-','HandleVisibility','off','LineWidth',2), hold off
        xlabel('Time (ps)'), set(gca,'xtick',[0 300 600 900])
        ylabel([{'Normalized.'},{'Intensity'}]), set(gca,'ytick',[0 1])
        % title(['la = ',num2str(10/mua),' mm']),
        legend('show','Location','northeast')
        axis([0 900 0 1])
        set(gca,'fontsize',18),
        clearvars mua mus the_filename t_db
    end
end
end

function [x_in, y_in, z_in, x_ot, y_ot, z_ot, s, w, no_of_photons, M_raw] = do_simulation (opt_bckg,opt_sgnl,cmd_size,sig_pos,sig_size,zSurface,wvlngth,nPhotonsReq,nExamplePaths,beam_X,beam_Y,beam_phi,beam_tht)
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
model.MC.lightSource.sourceType = 4;   % radial-factorizable beam
model.MC.lightSource.focalPlaneIntensityDistribution.radialDistr = 0;   % 0 = Top-hat
model.MC.lightSource.focalPlaneIntensityDistribution.radialWidth = 0.1; % [cm] radius
model.MC.lightSource.angularIntensityDistribution.radialDistr = 0;      % 0 = Top-hat angular
model.MC.lightSource.angularIntensityDistribution.radialWidth = 1.5e-3; % [rad] viewing half-angle :contentReference[oaicite:2]{index=2}

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
    disp(['progress ... ',num2str(size(scatter_stat,1)/nExamplePaths*100),'%'])
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
