function [] = Photon_09 ()
% Repository group: 01_monte_carlo_foundations
% Version role: primary simulation or study entry point
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% the study of Beer-Lambert law using MCmatlab
% implementing Monte Carlo Simulation - try yo replicate "Estimation of
% optical pathlength through tissue from direct time of flight measurement" figures
% Figure 1
clc

mu_s = [20,60,100]/log(10);
mu_a = [0.456, 0.334, 0.263, 0.217, 0.0867, 0.0434];
rest = nan(length(mu_s),length(mu_a));
for i_s = 1:length(mu_s)
    for i_a = 1:length(mu_a)
        [model] = do_simulation (mu_s(i_s),mu_a(i_a));
        rest(i_s,i_a) = sum(sum(model.MC.normalizedFluenceRate(:,:,end))./sum(sum(model.MC.normalizedFluenceRate(:,:,1))));
        clear model
    end
end
close all
figure(1)
plot(mu_a,-log10(rest(1,:)),'LineWidth',2,'Marker','+','DisplayName',['\mu_s = ',num2str(mu_s(1))]), hold on
plot(mu_a,-log10(rest(2,:)),'LineWidth',2,'Marker','x','DisplayName',['\mu_s = ',num2str(mu_s(2))]), hold on
plot(mu_a,-log10(rest(3,:)),'LineWidth',2,'Marker','s','DisplayName',['\mu_s = ',num2str(mu_s(3))]), hold on
hold off, grid on, xlim([0 0.5]), title('Figure 1'), xlabel('\mu_a cm^{-1}'), ylabel('Attenuation')
set(gca,'FontSize',16)
end

function [o_model] = do_simulation (mu_s,mu_a)
% Detectors: Placed at different distances and directions on the surface of the block
model = MCmatlab.model;
model.G.nx                = 1201;  % Number of bins in the x direction
model.G.ny                = 1201;  % Number of bins in the y direction
model.G.nz                = 100;   % Number of bins in the z direction
model.G.Lx                = 12.01; % [cm] x size of simulation cuboid
model.G.Ly                = 12.01; % [cm] y size of simulation cuboid
model.G.Lz                = 1.0;  % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.mediaPropParams     = {mu_s,mu_a};          % A cell array that you can use to contain all sorts of inputs you would like to use inside the mediaPropertiesFunc
model.G.geomFunc            = @geometryDefinition;  % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.
% model = plot(model,'G');

% Monte Carlo simulation
model.MC.simulationTimeRequested  = .1;    % [min] Time duration of the simulation
model.MC.matchedInterfaces        = true;  % Assumes all refractive indices are the same
model.MC.boundaryType             = 1;     % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.nPhotonsRequested = 5e4;          % The time to run the MC simualtion for, in minutes. The number of photos launched will vary from run to run.
model.MC.nExamplePaths = 1e3;              % The code store the paths of the first N photons for subsequent visualization during the plotting.
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
o_model = runMonteCarlo(model);
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
M(Z>=zSurface) = 2; % standard tissue
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
  mediaProperties(j).mua   = double(var{2}); % Absorption coefficient [cm^-1]
  mediaProperties(j).mus   = double(var{1}); % Scattering coefficient [cm^-1]
  mediaProperties(j).g     = 0.965; % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = 1.4; % The refractive index
end
