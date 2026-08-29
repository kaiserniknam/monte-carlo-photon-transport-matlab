function [] = Photon_08 ()
% Repository group: 01_monte_carlo_foundations
% Version role: primary simulation or study entry point
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% the study of Beer-Lambert law using MCmatlab
% implementing Monte Carlo Simulation - Find source to detector photon
% pathways when mu_s changes
clc

% healthy breast
var = {0,0,0,1,0.50,0.50,0.50};
for mu_s = [50,100,150,200,250,332.7]
    [d, N, N_wght, x_mean, y_mean, y_var, no_of_surf_photons, no_of_dpth_photons] = do_simulation (var,mu_s);

    figure(1), plot(d,N,'LineWidth',2,'DisplayName',['\mu_s = ',num2str(mu_s)])
    xlabel('d (cm)'), ylabel('probability of detected photons at s'), title('photo # vs. d')
    set(gca,'fontsize',16), axis tight, grid on, legend show, hold on

    figure(2), plot(d,log10(N_wght),'LineWidth',2,'DisplayName',['\mu_s = ',num2str(mu_s)])
    xlabel('d (cm)'), ylabel('weight of detected photons at d'), title('log_{10}(I/I_0) vs. d')
    set(gca,'fontsize',16), axis tight, grid on, legend show, hold on

    figure(3), plot(d,log10(N_wght)./d','LineWidth',2,'DisplayName',['\mu_s = ',num2str(mu_s)])
    xlabel('d (cm)'), ylabel('weight of detected photons at d / d'), title('log_{10}(I/I_0)/d vs. d')
    set(gca,'fontsize',16), axis tight, grid on, legend show, hold on

    figure(4), errorbar(d,y_mean,y_var,'LineWidth',2,'DisplayName',['\mu_s = ',num2str(mu_s)])
    xlabel('d (cm)'), ylabel('s (cm)'), title('s vs. d')
    set(gca,'fontsize',16), axis tight, grid on, legend show, hold on

    figure(5), plot(d,y_mean./x_mean,'LineWidth',2,'DisplayName',['\mu_s = ',num2str(mu_s)])
    xlabel('d (cm)'), ylabel('s/d (cm)'), title('~DPF')
    set(gca,'fontsize',16), axis tight, grid on, legend show, hold on
end
end

function [d, N, N_wght, x_mean, y_mean, y_var, no_of_surf_photons, no_of_dpth_photons] = do_simulation (var,mu_s)
% Detectors: Placed at different distances and directions on the surface of the block
model = MCmatlab.model;
model.G.nx                = 121;  % Number of bins in the x direction
model.G.ny                = 121;  % Number of bins in the y direction
model.G.nz                = 60;   % Number of bins in the z direction
model.G.Lx                = 12.1; % [cm] x size of simulation cuboid
model.G.Ly                = 12.1; % [cm] y size of simulation cuboid
model.G.Lz                = 6.0;  % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.mediaPropParams     = {mu_s};                % A cell array that you can use to contain all sorts of inputs you would like to use inside the mediaPropertiesFunc
model.G.geomFunc            = @geometryDefinition;  % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.
model.G.geomFuncParams      = var;                  % A cell array that you can use to contain all sorts of inputs you would like to use inside the geomFunc
% model = plot(model,'G');

% Monte Carlo simulation
model.MC.simulationTimeRequested  = .1;    % [min] Time duration of the simulation
model.MC.matchedInterfaces        = true;  % Assumes all refractive indices are the same
model.MC.boundaryType             = 1;     % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.nPhotonsRequested = 1e5;          % The time to run the MC simualtion for, in minutes. The number of photos launched will vary from run to run.
model.MC.nExamplePaths = 1e5;              % The code store the paths of the first N photons for subsequent visualization during the plotting.
model.MC.wavelength               = 800;   % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.lightSource.sourceType   = 0;    % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.silentMode = false;

% For a pencil beam, the "focus" is just a point that the beam goes
% through, here set to be the center of the cuboid:
model.MC.lightSource.xFocus       = 0.0;  % [cm] x position of focus
model.MC.lightSource.yFocus       = 0.0;  % [cm] y position of focus
model.MC.lightSource.zFocus       = 0.0;  % [cm] z position of focus
model.MC.lightSource.theta        = 0;    % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = 0;    % [rad] Azimuthal angle of beam center axis

% These lines will run the Monte Carlo simulation with the provided parameters and subsequently plot the results:
% figure(10)
scatter_stat = [];
no_of_surf_photons = 0;
no_of_dpth_photons = 0;
while size(scatter_stat,1) < 50000
    t_model = runMonteCarlo(model);
%     model = plot(model,'MC');
    % do calc
    spratrs = find(isnan(t_model.MC.examplePaths(1,:)));
    strt_pnts = spratrs(2:2:end); strt_pnts = strt_pnts(1:end-1);
    fnsh_pnts = spratrs(1:2:end); fnsh_pnts = fnsh_pnts(2:end-0);
    t_scatter_stat = nan(length(fnsh_pnts),5);
    for idx = 1:length(strt_pnts)
        Photon_Path = t_model.MC.examplePaths(:,strt_pnts(idx)+1:fnsh_pnts(idx)-1);
        if isnan(sum(sum(Photon_Path)))
            keyboard
        else
            if Photon_Path(3,end) <= 0.0
%                 plot3(Photon_Path(1,:),Photon_Path(2,:),-Photon_Path(3,:),'r'), hold on
                t_scatter_stat(idx,:) = [...
                    Photon_Path(1,end), ...
                    Photon_Path(2,end), ...
                    Photon_Path(3,end), ...
                    sum(sqrt(sum(((Photon_Path(1:3,2:end-0)-Photon_Path(1:3,1:end-1))).^2,1))), ...
                    Photon_Path(4,end)
                    ];
            else
%                 plot3(Photon_Path(1,:),Photon_Path(2,:),-Photon_Path(3,:),'k'), hold on
            end
        end
        clearvars Photon_Path
    end
    clearvars strt_pnts fnsh_pnts idx spratrs

    surf_photons = ~isnan(t_scatter_stat(:,1));
    no_of_surf_photons = no_of_surf_photons + sum( surf_photons);
    no_of_dpth_photons = no_of_dpth_photons + sum(~surf_photons);
    t_scatter_stat = t_scatter_stat(surf_photons,:);

    scatter_stat = [scatter_stat ; t_scatter_stat];
    clearvars surf_photons t_scatter_stat t_model
end
% figure(10), hold off, axis equal, view([-90 0]), grid on
% axis([min(model.G.x) max(model.G.x) min(model.G.y) max(model.G.y) -6 0])
% xlabel('x (cm)'), ylabel('y (cm)'), zlabel('z (cm)'), set(gca,'fontsize',16)

% histcount
x = sqrt(scatter_stat(:,1).^2+scatter_stat(:,2).^2);     % d
y = scatter_stat(:,4);                                   % s
N_xbin = 250; Xedges = linspace(min(x),max(x),N_xbin+1); % # of x_bins
[N,~,Ind] = histcounts(x,Xedges); N = N./sum(N);
N_wght = accumarray(Ind,scatter_stat(:,5),[],@sum); N_wght = N_wght/model.MC.nPhotonsRequested;
d = 1/2*(Xedges(1:end-1)+Xedges(2:end-0));
x_mean = accumarray(Ind,x,[],@mean);
y_mean = accumarray(Ind,y,[],@mean);
y_var  = accumarray(Ind,y,[],@var);
y_var(y_mean==0) = nan; y_mean(y_mean==0) = nan;
end
function M = geometryDefinition (X,Y,Z,parameters)
% Geometry function(s) (see readme for details)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% mediaPropertiesFunc) at each voxel location.
zSurface = 0.0;
M = ones(size(X)); % Air
if     double(parameters{1}) == 0
    M(Z>=zSurface) = 2; % Breast tissue
elseif double(parameters{1}) == 1
    M(Z>=zSurface) = 2; % Breast tissue
    M(sqrt( ...
        (X-double(parameters{2})).^2/double(parameters{5})^2+ ...
        (Y-double(parameters{3})).^2/double(parameters{6})^2+ ...
        (Z-double(parameters{4})).^2/double(parameters{7})^2) ...
        <=1) = 3;
else
    return
end
end
function mediaProperties = mediaPropertiesFunc (var)
% Media Properties function (see readme for details)
% The media properties function defines all the optical and thermal
% properties of the media involved by filling out and returning a
% "mediaProperties" array of "mediumProperties" objects with various
% properties. The j indices are the numbers that are referred to in the
% geometry function (in this case, 1 for "air" and 2 for "standard tissue")
% See the readme file or the examples for a list of properties you may
% specify. Most properties may be specified as a numeric constant or as
% function handles.
%
% The function must take one input; the cell array containing any special
% parameters you might specify above in the model file, for example
% parameters that you might loop over in a for loop. In most simulations
% this "parameters" cell array is empty. Dependence on wavelength is shown
% in examples 4 and 23. Dependence on excitation fluence rate FR,
% temperature T or fractional heat damage FD can be specified as in
% examples 12-15.

% Always leave the following line in place to initialize the
  % mediaProperties array:
  mediaProperties = MCmatlab.mediumProperties;
% Modeling optical fluence and diffuse reflectance distribution in normal
% and cancerous breast tissues exposed to planar and Gaussian NIR beam
% shapes using Monte Carlo simulation (2023)

  % Put in your own media property definitions below at 800 nm
  j=1;
  mediaProperties(j).name  = 'air';
  mediaProperties(j).mua   = 1e-4; % Absorption coefficient [cm^-1]
  mediaProperties(j).mus   = 1e-5; % Scattering coefficient [cm^-1]
  mediaProperties(j).g     = 1; % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = 1; % The refractive index

  j=2;
  mediaProperties(j).name  = 'breast tissue';
  mediaProperties(j).mua   = 0.55; % Absorption coefficient [cm^-1]
  mediaProperties(j).mus   = double(var{1}); % Scattering coefficient [cm^-1]
  mediaProperties(j).g     = 0.965; % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = 1.4; % The refractive index

  j=3;
  mediaProperties(j).name  = 'breast tumor';
  mediaProperties(j).mua   = 0.34; % Absorption coefficient [cm^-1]
  mediaProperties(j).mus   = 233.0; % Scattering coefficient [cm^-1]
  mediaProperties(j).g     = 0.962; % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = 1.4; % The refractive index
end
