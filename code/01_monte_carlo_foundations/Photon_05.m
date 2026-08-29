function [] = Photon_05 ()
% Repository group: 01_monte_carlo_foundations
% Version role: primary simulation or study entry point
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% the study of Beer-Lambert law using MCmatlab
% implementing Monte Carlo Simulation - The relation between log10(I/I0) vs. mu_s & mu_a

g = 0.9; n = 1.4;
mua = 0.05:0.1:0.95;
mus = 0.05:0.1:0.95;
detpos = [51, 51, 95, 2;
            51, 49, 95, 2;
            49, 51, 95, 2;
            51, 53, 95, 2;
            53, 51, 95, 2;

            51, 43, 95, 2;
            43, 51, 95, 2;
            51, 59, 95, 2;
            59, 51, 95, 2;

            51, 37, 95, 2;
            37, 51, 95, 2;
            51, 65, 95, 2;
            65, 51, 95, 2;

            51, 31, 95, 2;
            31, 51, 95, 2;
            51, 71, 95, 2;
            71, 51, 95, 2;
            ];
d = zeros(size(detpos,1),1);
for i_d = 2:17, d(i_d) = sqrt(sum((detpos(i_d,1:3)-detpos(1,1:3)).^2)); end
clearvars i_d

r_f = nan(length(mua),length(mus),size(detpos,1));
for i_a = 1:length(mua)
    for i_s = 1:length(mus)
        r_f(i_a,i_s,:) = get_measuremetns (mua(i_a), mus(i_s), g, n, detpos);
    end
end
save('Photon_Packet_Transport_in_a_Scattering_and_Absorbing_Medium_05.mat',"r_f","detpos","n","g","mus","mua")


close all

% Heatmaps of mu_a & mu_s
figure(1), subplot(2,3,2)
imagesc(mua,mus,log10(squeeze(r_f(:,:,01)'))), xlabel('\mu_a'), ylabel('\mu_s'), colormap jet, colorbar
title('log_{10}(I/I_0) at light source'), set(gca,"FontSize",16), axis equal, axis tight, axis xy
figure(1), subplot(2,3,1)
imagesc(mua,mus,log10(squeeze(mean(r_f(:,:,02:05),3)'))), xlabel('\mu_a'), ylabel('\mu_s'), colormap jet, colorbar
title('norm. log_{10}(I/I_0) at d = 2 mm from lit'), set(gca,"FontSize",16), axis equal, axis tight, axis xy
figure(1), subplot(2,3,3)
imagesc(mua,mus,log10(squeeze(mean(r_f(:,:,06:09),3)'))), xlabel('\mu_a'), ylabel('\mu_s'), colormap jet, colorbar
title('norm. log_{10}(I/I_0) at d = 8 mm from lit'), set(gca,"FontSize",16), axis equal, axis tight, axis xy
figure(1), subplot(2,2,3)
imagesc(mua,mus,log10(squeeze(mean(r_f(:,:,10:13),3)'))), xlabel('\mu_a'), ylabel('\mu_s'), colormap jet, colorbar
title('norm. log_{10}(I/I_0) at d = 14 mm from lit'), set(gca,"FontSize",16), axis equal, axis tight, axis xy
figure(1), subplot(2,2,4)
imagesc(mua,mus,log10(squeeze(mean(r_f(:,:,14:17),3)'))), xlabel('\mu_a'), ylabel('\mu_s'), colormap jet, colorbar
title('norm. log_{10}(I/I_0) at d = 20 mm from lit'), set(gca,"FontSize",16), axis equal, axis tight, axis xy

% Heatmaps of normalized mu_a & mu_s
figure(2), subplot(2,3,2)
imagesc(mua,mus,log10(squeeze(r_f(:,:,01)'))), xlabel('\mu_a'), ylabel('\mu_s'), colormap jet, colorbar
title('log_{10}(I/I_0) at light source'), set(gca,"FontSize",16), axis equal, axis tight, axis xy
figure(2), subplot(2,3,1)
imagesc(mua,mus,log10(squeeze(mean(r_f(:,:,02:05),3)'))/2), xlabel('\mu_a'), ylabel('\mu_s'), colormap jet, colorbar
title('log_{10}(I/I_0) at d = 2 mm from lit'), set(gca,"FontSize",16), axis equal, axis tight, axis xy
figure(2), subplot(2,3,3)
imagesc(mua,mus,log10(squeeze(mean(r_f(:,:,06:09),3)'))/8), xlabel('\mu_a'), ylabel('\mu_s'), colormap jet, colorbar
title('log_{10}(I/I_0) at d = 8 mm from lit'), set(gca,"FontSize",16), axis equal, axis tight, axis xy
figure(2), subplot(2,2,3)
imagesc(mua,mus,log10(squeeze(mean(r_f(:,:,10:13),3)'))/14), xlabel('\mu_a'), ylabel('\mu_s'), colormap jet, colorbar
title('log_{10}(I/I_0) at d = 14 mm from lit'), set(gca,"FontSize",16), axis equal, axis tight, axis xy
figure(2), subplot(2,2,4)
imagesc(mua,mus,log10(squeeze(mean(r_f(:,:,14:17),3)'))/20), xlabel('\mu_a'), ylabel('\mu_s'), colormap jet, colorbar
title('log_{10}(I/I_0) at d = 20 mm from lit'), set(gca,"FontSize",16), axis equal, axis tight, axis xy

% mu_s vs distance
figure(3)
for i_mua = 4
    for i_d = 1:size(detpos,1)
        subplot(1,2,1), plot(mus,r_f(i_mua,:,i_d),'Marker','o','MarkerFaceColor',[1,0,d(i_d)./max(d)],'MarkerEdgeColor','k'), hold on
    end
    hold off
end
end

function [flu_mat] = get_measuremetns (mua, mus, g, n, detpos)
% Detectors: Placed at different distances and directions on the surface of the block

MCmatlab.closeMCmatlabFigures();
model = MCmatlab.model;
model.G.nx                = 101; % Number of bins in the x direction
model.G.ny                = 101; % Number of bins in the y direction
model.G.nz                = 100; % Number of bins in the z direction
model.G.Lx                = 10.1; % [cm] x size of simulation cuboid
model.G.Ly                = 10.1; % [cm] y size of simulation cuboid
model.G.Lz                = 10.0; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.mediaPropParams = {mua, mus, g, n}; % Media properties defined as a function at the end of this file
model.G.geomFunc            = @geometryDefinition;  % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.
model = plot(model,'G');

% Monte Carlo simulation
model.MC.simulationTimeRequested  = .1; % [min] Time duration of the simulation
model.MC.matchedInterfaces        = true; % Assumes all refractive indices are the same
model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.nExamplePaths = 100;
% model.MC.wavelength               = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.lightSource.sourceType   = 0; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.silentMode = true;

% For a pencil beam, the "focus" is just a point that the beam goes
% through, here set to be the center of the cuboid:
model.MC.lightSource.xFocus       = 0; % [cm] x position of focus
model.MC.lightSource.yFocus       = 0; % [cm] y position of focus
model.MC.lightSource.zFocus       = 0.5; % [cm] z position of focus

model.MC.lightSource.theta        = 0; % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = 0; % [rad] Azimuthal angle of beam center axis


% These lines will run the Monte Carlo simulation with the provided
% parameters and subsequently plot the results:
model = runMonteCarlo(model);
% model = plot(model,'MC');
flu_mat = nan(size(detpos,1),1);
for i_abs = 1:size(detpos,1)
    flu_mat(i_abs) = model.MC.normalizedFluenceRate(detpos(i_abs,1),detpos(i_abs,2),101-detpos(i_abs,3))/100;
end

% model = plot(model,'MC');
end
function M = geometryDefinition (X,Y,Z,parameters)
% Geometry function(s) (see readme for details)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% mediaPropertiesFunc) at each voxel location.
  M = ones(size(X)); % Air
  M(:,:, 6:100) = 2;  % Material
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

  % Put in your own media property definitions below:
  j=1;
  mediaProperties(j).name  = 'air';
  mediaProperties(j).mua   = 1e-8; % Absorption coefficient [cm^-1]
  mediaProperties(j).mus   = 1e-8; % Scattering coefficient [cm^-1]
  mediaProperties(j).g     = 1; % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = 1; % The refractive index

  j=2;
  mediaProperties(j).name  = 'standard tissue';
  mediaProperties(j).mua   = double(var{1}); % Absorption coefficient [cm^-1]
  mediaProperties(j).mus   = double(var{2}); % Scattering coefficient [cm^-1]
  mediaProperties(j).g     = double(var{3}); % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = double(var{4}); % The refractive index
end
