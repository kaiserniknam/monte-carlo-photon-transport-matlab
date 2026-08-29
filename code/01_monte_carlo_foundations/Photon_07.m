function [] = Photon_07 ()
% Repository group: 01_monte_carlo_foundations
% Version role: primary simulation or study entry point
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% the study of Beer-Lambert law using MCmatlab
% implementing Monte Carlo Simulation - Find source to detector photon pathways

% healthy breast
var = {0,0,0,1,0.50,0.50,0.50}; [model,examplePaths] = do_simulation (var);

close all
% geometries
spratrs = find(isnan(examplePaths(1,:)));
strt_pnts = spratrs(2:2:end); strt_pnts = strt_pnts(1:end-1);
fnsh_pnts = spratrs(1:2:end); fnsh_pnts = fnsh_pnts(2:end-0);
scatter_stst = nan(length(fnsh_pnts),4);

for idx = 1:length(strt_pnts)
    Photon_Path = examplePaths(:,strt_pnts(idx)+1:fnsh_pnts(idx)-1);
    if isnan(sum(sum(Photon_Path)))
        keyboard
    else
        if Photon_Path(3,end) <= 0.0
            plot3(Photon_Path(1,:),Photon_Path(2,:),-Photon_Path(3,:),'r'), hold on
            scatter_stst(idx,:) = [...
                Photon_Path(1,end), ...
                Photon_Path(2,end), ...
                Photon_Path(3,end), ...
                sum(sqrt(sum(((Photon_Path(1:3,2:end-0)-Photon_Path(1:3,1:end-1))).^2,1)))];
        else
            plot3(Photon_Path(1,:),Photon_Path(2,:),-Photon_Path(3,:),'k'), hold on
        end
    end
    clearvars Photon_Path
end
figure(1), hold off, axis equal, view([-90 0]), grid on
axis([min(model.G.x) max(model.G.x) min(model.G.y) max(model.G.y) -6 0])
xlabel('x (cm)'), ylabel('y (cm)'), zlabel('z (cm)'), set(gca,'fontsize',16)
clearvars fnsh_pnts strt_pnts idx spratrs

surf_photons = ~isnan(scatter_stst(:,1));
no_of_surf_photons = sum( surf_photons);
no_of_dpth_photons = sum(~surf_photons);
scatter_stst = scatter_stst(surf_photons,:);
clearvars surf_photons

% scatterplot
figure(2)
plot(sqrt(scatter_stst(:,1).^2+scatter_stst(:,2).^2),scatter_stst(:,4),'ko','MarkerFaceColor','m','MarkerEdgeColor','k')
hold on, plot([0 max(sqrt(scatter_stst(:,1).^2+scatter_stst(:,2).^2))],[0 max(sqrt(scatter_stst(:,1).^2+scatter_stst(:,2).^2))],'k-.'), hold off
xlabel('d (cm)'), ylabel('s (cm)'), title('d vs. s')
set(gca,'fontsize',16), axis tight, grid on


% histcounts2
x = sqrt(scatter_stst(:,1).^2+scatter_stst(:,2).^2);     % d
y = scatter_stst(:,4);                                   % s
N_xbin = 100; Xedges = linspace(min(x),max(x),N_xbin+1); % # of x_bins
N_ybin = 101; Yedges = linspace(min(y),max(y),N_ybin+1); % # of x_bins
[N,~,~] = histcounts2(x,y,Xedges,Yedges);
figure(3), contourf(1/2*(Xedges(1:end-1)+Xedges(2:end-0)),1/2*(Yedges(1:end-1)+Yedges(2:end-0)),log10(N')), h = colorbar; h.Label.String = 'log_{10}(freq)';
hold on, plot([0 max(sqrt(scatter_stst(:,1).^2+scatter_stst(:,2).^2))],[0 max(sqrt(scatter_stst(:,1).^2+scatter_stst(:,2).^2))],'k-.'), hold off
xlabel('d (cm)'), ylabel('s (cm)'), title('d vs. s')
set(gca,'fontsize',16), axis tight, grid on
clearvars x y N_xbin N_ybin N h Xedges Yedges

% histcount
x = sqrt(scatter_stst(:,1).^2+scatter_stst(:,2).^2);     % d
y = scatter_stst(:,4);                                   % s
N_xbin = 150; Xedges = linspace(min(x),max(x),N_xbin+1); % # of x_bins
[N,~,Ind] = histcounts(x,Xedges); N = N./sum(N);
figure(4), plot(1/2*(Xedges(1:end-1)+Xedges(2:end-0)),N,'b','LineWidth',2)
xlabel('d (cm)'), ylabel('probability of detected photons at s'), title('d vs. s')
set(gca,'fontsize',16), axis tight, grid on
x_mean = accumarray(Ind,x,[],@mean);
y_mean = accumarray(Ind,y,[],@mean);
y_var  = accumarray(Ind,y,[],@var);
y_var(y_mean==0) = nan; y_mean(y_mean==0) = nan;
figure(5), errorbar(1/2*(Xedges(1:end-1)+Xedges(2:end-0)),y_mean,y_var,'b','LineWidth',2)
xlabel('d (cm)'), ylabel('s (cm)'), title('d vs. s')
set(gca,'fontsize',16), axis tight, grid on
figure(6), plot(1/2*(Xedges(1:end-1)+Xedges(2:end-0)),y_mean./x_mean,'b','LineWidth',2)
xlabel('d (cm)'), ylabel('s/d (cm)'), title('d vs. s')
set(gca,'fontsize',16), axis tight, grid on
end

function [model,examplePaths] = do_simulation (var)
% Detectors: Placed at different distances and directions on the surface of the block
MCmatlab.closeMCmatlabFigures();
model = MCmatlab.model;
model.G.nx                = 121;  % Number of bins in the x direction
model.G.ny                = 121;  % Number of bins in the y direction
model.G.nz                = 60;   % Number of bins in the z direction
model.G.Lx                = 12.1; % [cm] x size of simulation cuboid
model.G.Ly                = 12.1; % [cm] y size of simulation cuboid
model.G.Lz                = 6.0;  % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc            = @geometryDefinition;  % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.
model.G.geomFuncParams      = var;                  % A cell array that you can use to contain all sorts of inputs you would like to use inside the geomFunc
model = plot(model,'G');

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


% These lines will run the Monte Carlo simulation with the provided
% parameters and subsequently plot the results:
model = runMonteCarlo(model);
% model = plot(model,'MC');
examplePaths = model.MC.examplePaths;
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
function [binned_values,binned_counts,x_edges,y_edges] = do_2Dbining (scatter_stst,num_bins_x,num_bins_y)
% Extract x, y, and value columns
x = scatter_stst(:, 1);
y = scatter_stst(:, 2);
values = scatter_stst(:, 3);

% Define the bin edges (you can define specific ranges if necessary)
x_edges = linspace(min(x), max(x), num_bins_x+1);
y_edges = linspace(min(y), max(y), num_bins_y+1);

% Bin the data using histcounts2 to find which bin each point belongs to
[binned_counts, ~, ~, binX, binY] = histcounts2(x, y, x_edges, y_edges);

% Initialize a matrix to hold the binned values
binned_values = zeros(num_bins_x, num_bins_y);  % Size is num_bins_y x num_bins_x

% Loop through all data points and assign values to corresponding bins
for idx = 1:length(values)
    % Get bin index for the current x, y location
    x_bin_idx = binX(idx);
    y_bin_idx = binY(idx);
    binned_values(y_bin_idx, x_bin_idx) = binned_values(y_bin_idx, x_bin_idx) + values(idx);
end
end
