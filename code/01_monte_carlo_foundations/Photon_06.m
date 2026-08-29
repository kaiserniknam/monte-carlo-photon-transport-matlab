function [] = Photon_06 ()
% Repository group: 01_monte_carlo_foundations
% Version role: primary simulation or study entry point
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% the study of Beer-Lambert law using MCmatlab
% implementing Monte Carlo Simulation - Breast cancer (compressed & normal) vs Healthy Breast

% healthy breast
var = {0,0,0,1,0.50,0.50,0.50}; [srf_flu_ob,crs_abs_ob,PathsMap_ob,crs_geo_ob] = do_simulation (var);
% breast with ordinary tumor
var = {1,0,0,1,0.40,0.40,0.40}; [srf_flu_ot,crs_abs_ot,PathsMap_ot,crs_geo_ot] = do_simulation (var);
% breast with aordinary tumor
var = {1,0,0,1,0.40,0.40,0.25}; [srf_flu_ct,crs_abs_ct,PathsMap_ct,crs_geo_ct] = do_simulation (var);

close all
% geometries
figure(1), subplot(1,3,1), imagesc(linspace(-3.01/2,+3.01/2,301),0.01:0.01:1.5,crs_geo_ob.'), axis ij, colormap parula, caxis([1 3])
xlabel('y (cm)'), ylabel('z (cm)'), title('Healthy Breast'), set(gca,'fontsize',14), axis equal, axis tight
figure(1), subplot(1,3,2), imagesc(linspace(-3.01/2,+3.01/2,301),0.01:0.01:1.5,crs_geo_ot.'), axis ij, colormap parula, caxis([1 3])
xlabel('y (cm)'), ylabel('z (cm)'), title('Ordinary Tumor'), set(gca,'fontsize',14), axis equal, axis tight
figure(1), subplot(1,3,3), imagesc(linspace(-3.01/2,+3.01/2,301),0.01:0.01:1.5,crs_geo_ct.'), axis ij, colormap parula, caxis([1 3])
xlabel('y (cm)'), ylabel('z (cm)'), title('Compressed Tumor'), set(gca,'fontsize',14), axis equal, axis tight

% surface intensities
figure(2), subplot(1,3,1), imagesc(linspace(-3.01/2,+3.01/2,301),linspace(-3.01/2,+3.01/2,301),log10(srf_flu_ob.')), axis ij, colormap jet, caxis ([-10 2]), h = colorbar; h.Label.String = 'log_{10}(I/I_0)'; h.Label.FontSize = 16;
xlabel('x (cm)'), ylabel('y (cm)'), title('Healthy Breast'), set(gca,'fontsize',14), axis equal, axis tight
figure(2), subplot(1,3,2), imagesc(linspace(-3.01/2,+3.01/2,301),linspace(-3.01/2,+3.01/2,301),log10(srf_flu_ot.')), axis ij, colormap jet, caxis ([-10 2]), h = colorbar; h.Label.String = 'log_{10}(I/I_0)'; h.Label.FontSize = 16;
xlabel('x (cm)'), ylabel('y (cm)'), title('Ordinary Tumor'), set(gca,'fontsize',14), axis equal, axis tight
figure(2), subplot(1,3,3), imagesc(linspace(-3.01/2,+3.01/2,301),linspace(-3.01/2,+3.01/2,301),log10(srf_flu_ct.')), axis ij, colormap jet, caxis ([-10 2]), h = colorbar; h.Label.String = 'log_{10}(I/I_0)'; h.Label.FontSize = 16;
xlabel('x (cm)'), ylabel('y (cm)'), title('Compressed Tumor'), set(gca,'fontsize',14), axis equal, axis tight

% Delta surface intensities
figure(3), subplot(1,3,1), imagesc(linspace(-3.01/2,+3.01/2,301),linspace(-3.01/2,+3.01/2,301),abs(log10(srf_flu_ob.')-log10(srf_flu_ot.'))), axis ij, colormap jet, caxis([0 7]), h = colorbar; h.Label.String = '\Deltalog_{10}(I/I_0)'; h.Label.FontSize = 16;
xlabel('x (cm)'), ylabel('y (cm)'), title('Healthy - Tumor'), set(gca,'fontsize',14), axis equal, axis tight
figure(3), subplot(1,3,2), imagesc(linspace(-3.01/2,+3.01/2,301),linspace(-3.01/2,+3.01/2,301),abs(log10(srf_flu_ob.')-log10(srf_flu_ct.'))), axis ij, colormap jet, caxis([0 7]), h = colorbar; h.Label.String = '\Deltalog_{10}(I/I_0)'; h.Label.FontSize = 16;
xlabel('x (cm)'), ylabel('y (cm)'), title('Healthy - Compress'), set(gca,'fontsize',14), axis equal, axis tight
figure(3), subplot(1,3,3), imagesc(linspace(-3.01/2,+3.01/2,301),linspace(-3.01/2,+3.01/2,301),abs(log10(srf_flu_ot.')-log10(srf_flu_ct.'))), axis ij, colormap jet, caxis([0 7]), h = colorbar; h.Label.String = '\Deltalog_{10}(I/I_0)'; h.Label.FontSize = 16;
xlabel('x (cm)'), ylabel('y (cm)'), title('Tumor - Compress'), set(gca,'fontsize',14), axis equal, axis tight

% surface intensities
figure(4), subplot(2,3,1), imagesc(linspace(-3.01/2,+3.01/2,301),0.01:0.01:1.5,log10(crs_abs_ob.')), axis ij, colormap jet, caxis([-9 4]), h = colorbar; h.Label.String = 'log_{10}(I/I_0)'; h.Label.FontSize = 16;
hold on, contour(linspace(-3.01/2,+3.01/2,301),0.01:0.01:1.5,(crs_geo_ob.'),3,'color','k','linestyle','-.'), hold off, axis ij
xlabel('y (cm)'), ylabel('z (cm)'), title('Healthy Breast'), set(gca,'fontsize',14), axis equal, axis tight
figure(4), subplot(2,3,2), imagesc(linspace(-3.01/2,+3.01/2,301),0.01:0.01:1.5,log10(crs_abs_ot.')), axis ij, colormap jet, caxis([-9 4]), h = colorbar; h.Label.String = 'log_{10}(I/I_0)'; h.Label.FontSize = 16;
hold on, contour(linspace(-3.01/2,+3.01/2,301),0.01:0.01:1.5,(crs_geo_ot.'),3,'color','k','linestyle','-.'), hold off, axis ij
xlabel('y (cm)'), ylabel('z (cm)'), title('Ordinary Tumor'), set(gca,'fontsize',14), axis equal, axis tight
figure(4), subplot(2,3,3), imagesc(linspace(-3.01/2,+3.01/2,301),0.01:0.01:1.5,log10(crs_abs_ct.')), axis ij, colormap jet, caxis([-9 4]), h = colorbar; h.Label.String = 'log_{10}(I/I_0)'; h.Label.FontSize = 16;
hold on, contour(linspace(-3.01/2,+3.01/2,301),0.01:0.01:1.5,(crs_geo_ct.'),3,'color','k','linestyle','-.'), hold off, axis ij
xlabel('y (cm)'), ylabel('z (cm)'), title('Compressed Tumor'), set(gca,'fontsize',14), axis equal, axis tight

figure(4), subplot(2,3,4), imagesc(linspace(-3.01/2,+3.01/2,301),0.01:0.01:1.5,abs(log10(crs_abs_ob.')-log10(crs_abs_ot.'))), axis ij, colormap jet, caxis([0 6]), h = colorbar; h.Label.String = '\Deltalog_{10}(I/I_0)'; h.Label.FontSize = 16;
hold on, contour(linspace(-3.01/2,+3.01/2,301),0.01:0.01:1.5,(crs_geo_ob.'),3,'color','w','linestyle','-.'), hold off, axis ij
xlabel('y (cm)'), ylabel('z (cm)'), title('Healthy - Tumor'), set(gca,'fontsize',14), axis equal, axis tight
figure(4), subplot(2,3,5), imagesc(linspace(-3.01/2,+3.01/2,301),0.01:0.01:1.5,abs(log10(crs_abs_ob.')-log10(crs_abs_ct.'))), axis ij, colormap jet, caxis([0 6]), h = colorbar; h.Label.String = '\Deltalog_{10}(I/I_0)'; h.Label.FontSize = 16;
hold on, contour(linspace(-3.01/2,+3.01/2,301),0.01:0.01:1.5,(crs_geo_ot.'),3,'color','w','linestyle','-.'), hold off, axis ij
xlabel('y (cm)'), ylabel('z (cm)'), title('Healthy - Compress'), set(gca,'fontsize',14), axis equal, axis tight
figure(4), subplot(2,3,6), imagesc(linspace(-3.01/2,+3.01/2,301),0.01:0.01:1.5,abs(log10(crs_abs_ot.')-log10(crs_abs_ct.'))), axis ij, colormap jet, caxis([0 6]), h = colorbar; h.Label.String = '\Deltalog_{10}(I/I_0)'; h.Label.FontSize = 16;
hold on, contour(linspace(-3.01/2,+3.01/2,301),0.01:0.01:1.5,(crs_geo_ct.'),3,'color','w','linestyle','-.'), hold off, axis ij
xlabel('y (cm)'), ylabel('z (cm)'), title('Tumor - Compress'), set(gca,'fontsize',14), axis equal, axis tight

figure(5),
plot(linspace(-3.01/2,+3.01/2,301),log10(crs_abs_ob(:,51)),'DisplayName','Healthy Breast','LineWidth',2), hold on
plot(linspace(-3.01/2,+3.01/2,301),log10(crs_abs_ot(:,51)),'DisplayName','Ordinary Tumor','LineWidth',2), hold on
plot(linspace(-3.01/2,+3.01/2,301),log10(crs_abs_ct(:,51)),'DisplayName','Compressed Tumor','LineWidth',2), hold on
xlabel('y (cm)'), ylabel('log_{10}(I/I_0)'), title('recorded intensities'), set(gca,'fontsize',16), axis tight, legend('show')

% scattering paths
figure(6), subplot(2,3,1), contourf(linspace(-3.01/2,+3.01/2,301),0.01:0.01:1.5,log10(squeeze(PathsMap_ob(151,:,:)).')), axis ij, colormap jet, h = colorbar; h.Label.String = 'Pathway (Log_{10})'; h.Label.FontSize = 16; caxis([-2 4])
hold on, contour(linspace(-3.01/2,+3.01/2,301),0.01:0.01:1.5,(crs_geo_ob.'),3,'color','k','linestyle','-.'), hold off, axis ij
xlabel('y (cm)'), ylabel('z (cm)'), title('Healthy Breast'), set(gca,'fontsize',14), axis equal, axis tight
figure(6), subplot(2,3,2), contourf(linspace(-3.01/2,+3.01/2,301),0.01:0.01:1.5,log10(squeeze(PathsMap_ot(151,:,:)).')), axis ij, colormap jet, h = colorbar; h.Label.String = 'Pathway (Log_{10})'; h.Label.FontSize = 16; caxis([-2 4])
hold on, contour(linspace(-3.01/2,+3.01/2,301),0.01:0.01:1.5,(crs_geo_ot.'),3,'color','k','linestyle','-.'), hold off, axis ij
xlabel('y (cm)'), ylabel('z (cm)'), title('Ordinary Tumor'), set(gca,'fontsize',14), axis equal, axis tight
figure(6), subplot(2,3,3), contourf(linspace(-3.01/2,+3.01/2,301),0.01:0.01:1.5,log10(squeeze(PathsMap_ct(151,:,:)).')), axis ij, colormap jet, h = colorbar; h.Label.String = 'Pathway (Log_{10})'; h.Label.FontSize = 16; caxis([-2 4])
hold on, contour(linspace(-3.01/2,+3.01/2,301),0.01:0.01:1.5,(crs_geo_ct.'),3,'color','k','linestyle','-.'), hold off, axis ij
xlabel('y (cm)'), ylabel('z (cm)'), title('Compressed Tumor'), set(gca,'fontsize',14), axis equal, axis tight

figure(6), subplot(2,3,4), contourf(linspace(-3.01/2,+3.01/2,301),0.01:0.01:1.5,abs(log10(squeeze(PathsMap_ob(151,:,:)).')-log10(squeeze(PathsMap_ot(151,:,:)).'))), axis ij, colormap jet, h = colorbar; h.Label.String = 'Pathway (Log_{10})'; h.Label.FontSize = 16;
hold on, contour(linspace(-3.01/2,+3.01/2,301),0.01:0.01:1.5,(crs_geo_ob.'),3,'color','k','linestyle','-.'), hold off, axis ij
xlabel('y (cm)'), ylabel('z (cm)'), title('Healthy - Tumor'), set(gca,'fontsize',14), axis equal, axis tight
figure(6), subplot(2,3,5), contourf(linspace(-3.01/2,+3.01/2,301),0.01:0.01:1.5,abs(log10(squeeze(PathsMap_ob(151,:,:)).')-log10(squeeze(PathsMap_ct(151,:,:)).'))), axis ij, colormap jet, h = colorbar; h.Label.String = 'Pathway (Log_{10})'; h.Label.FontSize = 16;
hold on, contour(linspace(-3.01/2,+3.01/2,301),0.01:0.01:1.5,(crs_geo_ot.'),3,'color','k','linestyle','-.'), hold off, axis ij
xlabel('y (cm)'), ylabel('z (cm)'), title('Healthy - Compress'), set(gca,'fontsize',14), axis equal, axis tight
figure(6), subplot(2,3,6), contourf(linspace(-3.01/2,+3.01/2,301),0.01:0.01:1.5,abs(log10(squeeze(PathsMap_ot(151,:,:)).')-log10(squeeze(PathsMap_ct(151,:,:)).'))), axis ij, colormap jet, h = colorbar; h.Label.String = 'Pathway (Log_{10})'; h.Label.FontSize = 16;
hold on, contour(linspace(-3.01/2,+3.01/2,301),0.01:0.01:1.5,(crs_geo_ct.'),3,'color','k','linestyle','-.'), hold off, axis ij
xlabel('y (cm)'), ylabel('z (cm)'), title('Tumor - Compress'), set(gca,'fontsize',14), axis equal, axis tight
end

function [srf_flu,crs_abs,PathsMap,crs_geo] = do_simulation (var)
% Detectors: Placed at different distances and directions on the surface of the block

MCmatlab.closeMCmatlabFigures();
model = MCmatlab.model;
model.G.nx                = 301;  % Number of bins in the x direction
model.G.ny                = 301;  % Number of bins in the y direction
model.G.nz                = 150;  % Number of bins in the z direction
model.G.Lx                = 3.01; % [cm] x size of simulation cuboid
model.G.Ly                = 3.01; % [cm] y size of simulation cuboid
model.G.Lz                = 1.50; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc            = @geometryDefinition;  % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.
model.G.geomFuncParams      = var; % A cell array that you can use to contain all sorts of inputs you would like to use inside the geomFunc
model = plot(model,'G');

% Monte Carlo simulation
model.MC.simulationTimeRequested  = .1;   % [min] Time duration of the simulation
model.MC.matchedInterfaces        = true; % Assumes all refractive indices are the same
model.MC.boundaryType             = 1;    % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.nPhotonsRequested = 5e4;         % The time to run the MC simualtion for, in minutes. The number of photos launched will vary from run to run.
model.MC.nExamplePaths = 5e4;             % The code store the paths of the first N photons for subsequent visualization during the plotting.
model.MC.wavelength               = 800;  % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.lightSource.sourceType   = 0;    % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.silentMode = false;

% For a pencil beam, the "focus" is just a point that the beam goes
% through, here set to be the center of the cuboid:
model.MC.lightSource.xFocus       = 0.0;  % [cm] x position of focus
model.MC.lightSource.yFocus       = 0.5;  % [cm] y position of focus
model.MC.lightSource.zFocus       = 0.5;  % [cm] z position of focus
model.MC.lightSource.theta        = 0;    % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = 0;    % [rad] Azimuthal angle of beam center axis


% These lines will run the Monte Carlo simulation with the provided
% parameters and subsequently plot the results:
model = runMonteCarlo(model);
% model = plot(model,'MC');
crs_geo = squeeze(model.G.M_raw(find(model.G.x==0),:,:));
srf_flu = model.MC.normalizedFluenceRate(:,:,6)/100;
crs_abs = squeeze(model.MC.normalizedAbsorption(find(model.G.x==0),:,:));
examplePaths = model.MC.examplePaths;
idx = 1; examplePaths(idx,:) = round(examplePaths(idx,:)/model.G.dx) + find(model.G.x==0);              examplePaths(idx,examplePaths(idx,:)==0) = nan; examplePaths(idx,examplePaths(idx,:)>size(model.G.M_raw,idx)) = nan;
idx = 2; examplePaths(idx,:) = round(examplePaths(idx,:)/model.G.dy) + find(model.G.y==0);              examplePaths(idx,examplePaths(idx,:)==0) = nan; examplePaths(idx,examplePaths(idx,:)>size(model.G.M_raw,idx)) = nan;
idx = 3; examplePaths(idx,:) = round(examplePaths(idx,:)/model.G.dz) + find(model.G.z==min(model.G.z)); examplePaths(idx,examplePaths(idx,:)==0) = nan; examplePaths(idx,examplePaths(idx,:)>size(model.G.M_raw,idx)) = nan;
examplePaths(:,isnan(sum(examplePaths(1:3,:),1))) = [];
PathsMap = zeros(size(model.G.M_raw));
for i_examplePaths = 1:size(examplePaths,2)
    PathsMap(examplePaths(1,i_examplePaths),examplePaths(2,i_examplePaths),examplePaths(3,i_examplePaths)) = ...
    PathsMap(examplePaths(1,i_examplePaths),examplePaths(2,i_examplePaths),examplePaths(3,i_examplePaths)) + ...
        examplePaths(4,i_examplePaths);
end
end
function M = geometryDefinition (X,Y,Z,parameters)
% Geometry function(s) (see readme for details)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% mediaPropertiesFunc) at each voxel location.
zSurface = 0.5;
M = ones(size(X)); % Air
if     double(parameters{1}) == 0
    M(Z>zSurface) = 2; % Breast tissue
elseif double(parameters{1}) == 1
    M(Z>zSurface) = 2; % Breast tissue
    M(sqrt( ...
        (X-double(parameters{2})).^2/double(parameters{5})^2+ ...
        (Y-double(parameters{3})).^2/double(parameters{6})^2+ ...
        (Z-double(parameters{4})).^2/double(parameters{7})^2)<=1) = 3;
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
  mediaProperties(j).mus   = 332.7; % Scattering coefficient [cm^-1]
  mediaProperties(j).g     = 0.965; % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = 1.4; % The refractive index

  j=3;
  mediaProperties(j).name  = 'breast tumor';
  mediaProperties(j).mua   = 0.34; % Absorption coefficient [cm^-1]
  mediaProperties(j).mus   = 233.0; % Scattering coefficient [cm^-1]
  mediaProperties(j).g     = 0.962; % Henyey-Greenstein scattering anisotropy
  mediaProperties(j).n     = 1.4; % The refractive index
end
