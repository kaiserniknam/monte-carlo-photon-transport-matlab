function [] = Photon_16 ()
% Repository group: 01_monte_carlo_foundations
% Version role: primary simulation or study entry point
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% the study of Beer-Lambert law using MCmatlab
% implementing Monte Carlo Simulation - Find source to detector photon (DPF) pathways when mu_s/mu_a changes
% the same as
% "Photon_Packet_Transport_in_a_Scattering_and_Absorbing_Medium_13" but to validate this code
% a layer of air (0.5 cm) is considered to see what happenes: sanity check

clc
close all

% healthy breast
var = {0,0,0,1,0.50,0.50,0.50};
g = 0.9; n = 1.4; do_count = 0;
mua = [0.1];
mus = [100];
% mua = [0.1];
% mus = [50,100,150,200,250,300];
mua = [0.1,0.2,0.3,0.4,0.5];
mus = [100];

for i_a = 1:length(mua)
    for i_s = 1:length(mus)
        do_count = do_count + 1;
        [d, s, w, no_of_surf_photons, no_of_dpth_photons] = do_simulation (mua(i_a),mus(i_s),g,n,var);

        N_bin = 250;
        s_edges = linspace(min(s),max(s),N_bin+1); % # s bins
        s_c = 1/2*(s_edges(1:end-1)+s_edges(2:end-0))';
        [f_s,~,ind_s] = histcounts(s,s_edges); f_s = f_s./sum(f_s); clearvars s_edges
        d_edges = linspace(min(d),max(d),N_bin+1); % # d bins
        d_c = 1/2*(d_edges(1:end-1)+d_edges(2:end-0))';
        [f_d,~,ind_d] = histcounts(d,d_edges); f_d = f_d./sum(f_d); clearvars d_edges

        w_mean = accumarray(ind_d,w,[],@sum); w_mean = w_mean/(no_of_surf_photons+no_of_dpth_photons);
        s_mean = accumarray(ind_d,s,[],@mean);
        s_var  = accumarray(ind_d,s,[],@var);
        s_var(s_mean==0) = nan; s_mean(s_mean==0) = nan;
        clearvars ind_s ind_d

        figure(1),
        if do_count==1
            plot([0 max(d)],[0 max(d)],'k-.','DisplayName','1:1 line'), hold on
        end
        plot(d,s,'o','MarkerFaceColor',get_color(do_count),'MarkerEdgeColor','k','DisplayName',['\mu_s = ',num2str(mus(i_s)),', \mu_a = ',num2str(mua(i_a))])
        xlabel('d (cm)'), ylabel('s (cm)'), title('d vs. s')
        set(gca,'fontsize',24), axis tight, grid on

        figure(2), plot(s,-log10(w),'o','MarkerFaceColor',get_color(do_count),'MarkerEdgeColor','k','DisplayName',['\mu_s = ',num2str(mus(i_s)),', \mu_a = ',num2str(mua(i_a))])
        xlabel('s (cm)'), ylabel('OD (unitless)'), title('OD vs. s (a sanity check)')
        set(gca,'fontsize',24), axis tight, grid on, legend('show','Location','northwest'), hold on

        figure(3), plot(d_c,-log10(w_mean),'LineWidth',2,'DisplayName',['\mu_s = ',num2str(mus(i_s)),', \mu_a = ',num2str(mua(i_a))])
        xfit = d_c; yfit = -log10(w_mean); xfit = xfit(~isinf(yfit)); yfit = yfit(~isinf(yfit));
        Bfit = [ones(size(xfit)) xfit]\yfit; % Estimate Parameters
        yhat = [ones(size(xfit)) xfit]*Bfit; % Calculate Regression Line
        if length(mua)==1 && length(mus)==1
            hold on, plot(xfit,yhat,'k-.','LineWidth',2,'DisplayName','linear fit'), clearvars xfit yfit Bfit yhat
        end
        xlabel('d (cm)'), ylabel('OD (unitless)'), title('OD vs. d')
        set(gca,'fontsize',24), axis tight, grid on, legend('show','Location','northwest'), hold on

        figure(4), loglog(d_c,-log10(w_mean)./d_c,'LineWidth',2,'DisplayName',['\mu_s = ',num2str(mus(i_s)),', \mu_a = ',num2str(mua(i_a))])
        xlabel('d (cm)'), ylabel('OD/d (cm^{-1})'), title('OD/d vs. d')
        set(gca,'fontsize',24), axis tight, grid on, legend show, hold on

        figure(5), errorbar(d_c,s_mean,s_var,'LineWidth',2,'DisplayName',['\mu_s = ',num2str(mus(i_s)),', \mu_a = ',num2str(mua(i_a))])
        xlabel('d (cm)'), ylabel('s (cm)'), title('s vs. d')
        set(gca,'fontsize',24), axis tight, grid on, legend('show','Location','northwest'), hold on

        figure(6), plot(d_c,s_mean,'LineWidth',2,'DisplayName',['\mu_s = ',num2str(mus(i_s)),', \mu_a = ',num2str(mua(i_a))])
        xfit = d_c; yfit = s_mean; xfit = xfit(~isnan(yfit)); yfit = yfit(~isnan(yfit));
        Bfit = [ones(size(xfit)) xfit]\yfit; % Estimate Parameters
        yhat = [ones(size(xfit)) xfit]*Bfit; % Calculate Regression Line
        if length(mua)==1 && length(mus)==1
            hold on, plot(xfit,yhat,'k-.','LineWidth',2,'DisplayName','linear fit'), clearvars xfit yfit Bfit yhat
        end
        xlabel('d (cm)'), ylabel('s (cm)'), title('s vs. d')
        set(gca,'fontsize',24), axis tight, grid on, legend('show','Location','northwest'), hold on

        figure(7), plot(d_c,s_mean./d_c,'LineWidth',2,'DisplayName',['\mu_s = ',num2str(mus(i_s)),', \mu_a = ',num2str(mua(i_a))])
        xlabel('d (cm)'), ylabel('s/d (cm)'), title('~ DPF')
        set(gca,'fontsize',24), axis tight, grid on, legend show, hold on

        clearvars d s w d_c s_c s_mean w_mean s_var f_d f_s no_of_dpth_photons no_of_surf_photons
    end
end
end

function [d, s, w, no_of_surf_photons, no_of_dpth_photons] = do_simulation (mua,mus,g,n,var)
z_surf = 0.5;
model = MCmatlab.model;
model.G.nx                = 121;  % Number of bins in the x direction
model.G.ny                = 121;  % Number of bins in the y direction
model.G.nz                = 60;   % Number of bins in the z direction
model.G.Lx                = 12.1; % [cm] x size of simulation cuboid
model.G.Ly                = 12.1; % [cm] y size of simulation cuboid
model.G.Lz                = 6.0;  % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.mediaPropParams     = {mua,mus,g,n};                % A cell array that you can use to contain all sorts of inputs you would like to use inside the mediaPropertiesFunc
model.G.geomFunc            = @geometryDefinition;  % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.
model.G.geomFuncParams      = var;                  % A cell array that you can use to contain all sorts of inputs you would like to use inside the geomFunc
% model = plot(model,'G');
clearvars mua mus g n var

% Monte Carlo simulation
model.MC.simulationTimeRequested  = .1;    % [min] Time duration of the simulation
model.MC.matchedInterfaces        = true;  % Assumes all refractive indices are the same
model.MC.boundaryType             = 1;     % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.nPhotonsRequested = 1e5;          % The time to run the MC simualtion for, in minutes. The number of photos launched will vary from run to run.
model.MC.nExamplePaths = 1e5;              % The code store the paths of the first N photons for subsequent visualization during the plotting.
model.MC.wavelength               = 800;   % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.lightSource.sourceType   = 0;    % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.silentMode = false;

% For a pencil beam, the "focus" is just a point that the beam goes through, here set to be the center of the cuboid:
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
    % t_model = plot(t_model,'MC');
    % do calc
    spratrs = find(isnan(t_model.MC.examplePaths(1,:)));
    strt_pnts = spratrs(2:2:end); strt_pnts = strt_pnts(1:end-1);
    fnsh_pnts = spratrs(1:2:end); fnsh_pnts = fnsh_pnts(2:end-0);
    t_scatter_stat = nan(length(fnsh_pnts),5); clearvars spratrs
    for idx = 1:length(strt_pnts)
        Photon_Path = t_model.MC.examplePaths(:,strt_pnts(idx)+1:fnsh_pnts(idx)-1);
        if isnan(sum(sum(Photon_Path)))
            keyboard
        else
            if Photon_Path(3,end) <= z_surf
                [~,i_sort] = find(Photon_Path(3,:)>=z_surf,1,'first');
                Photon_Path = Photon_Path(:,i_sort+0:end); clearvars i_sort
                [~,i_sort] = find(Photon_Path(3,:)>=z_surf,1,'last');
                Photon_Path = Photon_Path(:,1:i_sort+1); clearvars i_sort
                % plot3(Photon_Path(1,:),Photon_Path(2,:),-Photon_Path(3,:),'r'), hold on
                t_scatter_stat(idx,:) = [...
                    Photon_Path(1,end), ...
                    Photon_Path(2,end), ...
                    Photon_Path(3,end), ...
                    sum(sqrt(sum(((Photon_Path(1:3,2:end-0)-Photon_Path(1:3,1:end-1))).^2,1))), ...
                    Photon_Path(4,end)
                    ];
            else
                % plot3(Photon_Path(1,:),Photon_Path(2,:),-Photon_Path(3,:),'k'), hold on
            end
        end
        clearvars Photon_Path
    end
    clearvars strt_pnts fnsh_pnts idx t_model

    surf_photons = ~isnan(t_scatter_stat(:,1));
    no_of_surf_photons = no_of_surf_photons + sum( surf_photons);
    no_of_dpth_photons = no_of_dpth_photons + sum(~surf_photons);
    t_scatter_stat = t_scatter_stat(surf_photons,:);

    scatter_stat = [scatter_stat ; t_scatter_stat];
    clearvars surf_photons t_scatter_stat
end
% figure(10), hold off, axis equal, view([-90 0]), grid on
% axis([min(model.G.x) max(model.G.x) min(model.G.y) max(model.G.y) -6 0])
% xlabel('x (cm)'), ylabel('y (cm)'), zlabel('z (cm)'), set(gca,'fontsize',24)

% histcount
d = sqrt(scatter_stat(:,1).^2+scatter_stat(:,2).^2);     % d -> source to detector distance
s = scatter_stat(:,4);                                   % s -> true distance
w = scatter_stat(:,5);                                   % w -> intensity
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
    M(Z>=zSurface) = 2; % standard tissue
elseif double(parameters{1}) == 1
    M(Z>=zSurface) = 2; % standard tissue
    M(sqrt( ...
        ((X-double(parameters{2})).^2)./(double(parameters{5}).^2) + ...
        ((Y-double(parameters{3})).^2)./(double(parameters{6}).^2) + ...
        ((Z-double(parameters{4})).^2)./(double(parameters{7}).^2))  ...
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
mediaProperties(j).mua   = double(var{1}); % Absorption coefficient [cm^-1]
mediaProperties(j).mus   = double(var{2}); % Scattering coefficient [cm^-1]
mediaProperties(j).g     = double(var{3}); % Henyey-Greenstein scattering anisotropy
mediaProperties(j).n     = double(var{4}); % The refractive index
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
