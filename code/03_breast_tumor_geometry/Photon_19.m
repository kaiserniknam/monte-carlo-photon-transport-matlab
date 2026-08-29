  function [] = Photon_19 ()
% Repository group: 03_breast_tumor_geometry
% Version role: primary simulation or study entry point
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: s vs. d for breast tissue
% ref. #1: Assessment of the size, position, and optical properties of breast tumors in vivo by noninvasive optical methods
% ref. #2: Quantitative Absorption and Scattering Spectra in Thick Tissues Using Broadband Diffuse Optical Spectroscopy
% ref. #3: Tumor location of the central and nipple portion is associated with impaired survival for women with breast cancer
% ref. #4: Broadband Optical Mammography: Chromophore Concentration and Hemoglobin Saturation Contrast in Breast Cancer
% a template

clc
close all

z_air = 0.0; % the thickness of air layer
wvlnt_sample = 800; % sample wavelength (nm)
nPhotonsReq = 1e5; % number of requested photon
nExamplePaths = 50000; % number of example photon paths
% Optical properties of breast and tumor [1-4]
[~,~,~,~,brst_strct,tumr_strct] = give_breast_mua_mus(wvlnt_sample,false);
Lx = 19.1; Ly = 19.1; Lz = 10.0; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
cx = 0.00; cy = 0.00; cz = 1.50; % Tumor depth from the skin (in cm)(0.5 to 5.0 cm)
rx = 1.00; ry = 1.00; rz = 1.00; % Tumor radius within the breast (in cm)(1.0 to 2.5 cm).
do_count = 0;

for idx = 1:length(brst_strct)
    % set_of_mus = linspace(brst_strct(idx).mus*0.8,brst_strct(idx).mus*1.2,9);
    % for i_mus = 1:length(set_of_mus)
    do_count = do_count + 1;
    mua_brst = brst_strct(idx).mua; g_brst = brst_strct(idx).g; n_brst = brst_strct(idx).n; mus_brst = brst_strct(idx).mus/(1-g_brst);
    mua_tumr = tumr_strct(idx).mua; g_tumr = tumr_strct(idx).g; n_tumr = tumr_strct(idx).n; mus_tumr = tumr_strct(idx).mus/(1-g_tumr); wvlnt = tumr_strct(idx).lambda;
    [d, s, w, a, no_of_surf_photons, no_of_dpth_photons] = do_simulation ( ...
        [mua_brst,mus_brst,g_brst,n_brst], ...
        [mua_tumr,mus_tumr,g_tumr,n_tumr], ...
        [Lx,Ly,Lz],[cx,cy,cz],[rx,ry,rz],z_air,wvlnt,nPhotonsReq,nExamplePaths);

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
    plot(d,s,'o','MarkerFaceColor',get_color(do_count),'MarkerEdgeColor','k','DisplayName',['\mu_{s,b} = ',num2str(mus_brst),', \mu_{a,b} = ',num2str(mua_brst),', \mu_{s,t} = ',num2str(mus_tumr),', \mu_{a,t} = ',num2str(mua_tumr)])
    xlabel('d (cm)'), ylabel('s (cm)'), title('d vs. s')
    set(gca,'fontsize',24)% colormap(my_colormap); , axis tight, grid on
    x_min = xlim;
    plot([min(x_min) max(x_min)],(-log10(0.9)/mua_brst*100).*[1 1],'-.','Color',get_color(do_count),'HandleVisibility','off','LineWidth',2)
    clearvars x_min

    u = round(-log10(w)./s,4); u_unique = unique(u); freq = nan(size(u_unique));
    for i_f = 1:length(u_unique), freq(i_f) = sum(u==u_unique(i_f)); end
    [~,i_x] = max(freq); s_f = s; w_f = w; s_f(u~=u_unique(i_x)) = []; w_f(u~=u_unique(i_x)) = [];
    figure(2), plot(s_f,-log10(w_f),'o','MarkerFaceColor',get_color(do_count),'MarkerEdgeColor','k','DisplayName',['\mu_{s,b} = ',num2str(mus_brst),', \mu_{a,b} = ',num2str(mua_brst),', \mu_{s,t} = ',num2str(mus_tumr),', \mu_{a,t} = ',num2str(mua_tumr)]), hold on
    xlabel('s (cm)'), ylabel('OD (unitless)'), title('OD vs. s (a sanity check)')
    set(gca,'fontsize',24), axis tight, grid on, legend('show','Location','northwest'), hold on
    clearvars u u_unique freq i_f i_x s_f w_f

    figure(3), plot(d_c,-log10(w_mean),'LineWidth',2,'DisplayName',['\mu_{s,b} = ',num2str(mus_brst),', \mu_{a,b} = ',num2str(mua_brst),', \mu_{s,t} = ',num2str(mus_tumr),', \mu_{a,t} = ',num2str(mua_tumr)],'Color',get_color(do_count)), hold on
    if isscalar(mua_brst) && isscalar(mus_brst)
        xfit = d_c; yfit = -log10(w_mean); xfit = xfit(~isinf(yfit)); yfit = yfit(~isinf(yfit));
        Bfit = [ones(size(xfit)) xfit]\yfit; % Estimate Parameters
        yhat = [ones(size(xfit)) xfit]*Bfit; % Calculate Regression Line
        hold on, plot(xfit,yhat,'-.','LineWidth',2,'DisplayName','linear fit','Color',get_color(do_count)), clearvars xfit yfit Bfit yhat
    end
    xlabel('d (cm)'), ylabel('OD (unitless)'), title('OD vs. d')
    set(gca,'fontsize',24), axis tight, grid on, legend('show','Location','northwest'), hold on

    figure(4), loglog(d_c,-log10(w_mean)./d_c,'LineWidth',2,'DisplayName',['\mu_{s,b} = ',num2str(mus_brst),', \mu_{a,b} = ',num2str(mua_brst),', \mu_{s,t} = ',num2str(mus_tumr),', \mu_{a,t} = ',num2str(mua_tumr)],'Color',get_color(do_count)), hold on
    xlabel('d (cm)'), ylabel('OD/d (cm^{-1})'), title('OD/d vs. d')
    set(gca,'fontsize',24), axis tight, grid on, legend show, hold on

    figure(5), errorbar(d_c,s_mean,s_var,'LineWidth',2,'DisplayName',['\mu_{s,b} = ',num2str(mus_brst),', \mu_{a,b} = ',num2str(mua_brst),', \mu_{s,t} = ',num2str(mus_tumr),', \mu_{a,t} = ',num2str(mua_tumr)],'Color',get_color(do_count)), hold on
    xlabel('d (cm)'), ylabel('s (cm)'), title('s vs. d')
    set(gca,'fontsize',24), axis tight, grid on, legend('show','Location','northwest'), hold on

    figure(6), plot(d_c,s_mean,'LineWidth',2,'DisplayName',['\mu_{s,b} = ',num2str(mus_brst),', \mu_{a,b} = ',num2str(mua_brst),', \mu_{s,t} = ',num2str(mus_tumr),', \mu_{a,t} = ',num2str(mua_tumr)],'Color',get_color(do_count)), hold on
    if isscalar(mua_brst) && isscalar(mus_brst)% colormap(my_colormap);
        xfit = d_c; yfit = s_mean; xfit = xfit(~isnan(yfit)); yfit = yfit(~isnan(yfit));
        Bfit = [ones(size(xfit)) xfit]\yfit; % Estimate Parameters
        yhat = [ones(size(xfit)) xfit]*Bfit; % Calculate Regression Line
        hold on, plot(xfit,yhat,'-.','LineWidth',2,'DisplayName','linear fit','Color',get_color(do_count)), clearvars xfit yfit Bfit yhat
    end
    xlabel('d (cm)'), ylabel('s (cm)'), title('s vs. d')
    set(gca,'fontsize',24), axis tight, grid on, legend('show','Location','northwest'), hold on

    figure(7), plot(d_c,s_mean./d_c,'LineWidth',2,'DisplayName',['\mu_{s,b} = ',num2str(mus_brst),', \mu_{a,b} = ',num2str(mua_brst),', \mu_{s,t} = ',num2str(mus_tumr),', \mu_{a,t} = ',num2str(mua_tumr)],'Color',get_color(do_count)), hold on
    xlabel('d (cm)'), ylabel('s/d (cm)'), title('~ DPF')
    set(gca,'fontsize',24), axis tight, grid on, legend show, hold on

    figure(8), plot(d_c,s_mean./d_c,'LineWidth',2,'DisplayName',['\mu_{s,b} = ',num2str(mus_brst),', \mu_{a,b} = ',num2str(mua_brst),', \mu_{s,t} = ',num2str(mus_tumr),', \mu_{a,t} = ',num2str(mua_tumr)],'Color',get_color(do_count)), hold on
    musp_brst = mus_brst*(1-g_brst);
    figure(8), plot(d_c,(sqrt(3*musp_brst)/2/sqrt(mua_brst)).*(d_c.*sqrt(3*musp_brst*mua_brst))./(d_c.*sqrt(3*musp_brst*mua_brst)+1),'LineWidth',2,'LineStyle','-.','Color',get_color(do_count),'HandleVisibility','off'), hold on
    xlabel('d (cm)'), ylabel('s/d (cm)'), title('~ DPF')
    set(gca,'fontsize',24), axis tight, grid on, legend show, hold on

    clearvars d s a w d_c s_c s_mean w_mean i_mean s_var f_d f_s no_of_dpth_photons no_of_surf_photons
    clearvars mua_brst g_brst n_brst mus_brst
    clearvars mua_tumr g_tumr n_tumr mus_tumr
end
end

function [d, s, w, a, no_of_surf_photons, no_of_dpth_photons] = do_simulation (opt_bckg,opt_sgnl,cmd_size,sig_pos,sig_size,zSurface,wvlngth,nPhotonsReq,nExamplePaths)
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
model.MC.lightSource.sourceType   = 0;    % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.silentMode = false;
model.MC.lightSource.xFocus       = 0.0;       % [cm] x position of focus
model.MC.lightSource.yFocus       = 0.0;       % [cm] y position of focus
model.MC.lightSource.zFocus       = zSurface;  % [cm] z position of focus
model.MC.lightSource.theta        = 0.0;       % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = 0.0;       % [rad] Azimuthal angle of beam center axis
clearvars zSurface

%% These lines will run the Monte Carlo simulation with the provided parameters and subsequently plot the results:
% figure(17)
scatter_stat = [];
no_of_surf_photons = 0;
no_of_dpth_photons = 0;
while size(scatter_stat,1) < nExamplePaths
    t_model = runMonteCarlo(model);
    % t_model = plot(t_model,'MC');
    % do calc
    spratrs = find(isnan(t_model.MC.examplePaths(1,:)));
    strt_pnts = spratrs(2:2:end); strt_pnts = strt_pnts(1:end-1);
    fnsh_pnts = spratrs(1:2:end); fnsh_pnts = fnsh_pnts(2:end-0);
    t_scatter_stat = nan(length(fnsh_pnts),3+1+1+3); clearvars spratrs
    for idx = 1:length(strt_pnts)
        Photon_Path = t_model.MC.examplePaths(:,strt_pnts(idx)+1:fnsh_pnts(idx)-1);
        if isnan(sum(sum(Photon_Path)))
            keyboard
        else
            if Photon_Path(3,end) <= 0.0
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
d = sqrt(scatter_stat(:,1).^2+scatter_stat(:,2).^2);     % d -> source to detector distance
s = scatter_stat(:,4);                                   % s -> true distance
w = scatter_stat(:,5);                                   % w -> intensity
a = scatter_stat(:,8)./sqrt(scatter_stat(:,6).^2 + scatter_stat(:,7).^2 + scatter_stat(:,8).^2); % a -> angle
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
function [mua_brst,mua_tumr,mus_brst,mus_tumr,brst_strct,tumr_strct] = give_breast_mua_mus(lambda_in,do_fig)
% data
tumr_db = readtable('DB/breast.tumor.csv');
brst_db = readtable('DB/breast.normal.csv');
brst_strct(1).lambda = 690; brst_strct(1).mua = 0.030; brst_strct(1).mus = 12.0; brst_strct(1).g = 0.95; brst_strct(1).n = 1.4;
brst_strct(2).lambda = 825; brst_strct(2).mua = 0.040; brst_strct(2).mus = 11.0; brst_strct(2).g = 0.95; brst_strct(2).n = 1.4;
tumr_strct(1).lambda = 690; tumr_strct(1).mua = 0.084; tumr_strct(1).mus = 15.0; tumr_strct(1).g = 0.95; tumr_strct(1).n = 1.4;
tumr_strct(2).lambda = 825; tumr_strct(2).mua = 0.085; tumr_strct(2).mus = 12.7; tumr_strct(2).g = 0.95; tumr_strct(2).n = 1.4;
% fit
bs_brst = -log(brst_strct(2).mus/brst_strct(1).mus)/log(brst_strct(2).lambda/brst_strct(1).lambda); as_brst = brst_strct(1).mus/brst_strct(1).lambda^(-bs_brst); mus_brst = as_brst*(lambda_in)^(-bs_brst); % clearvars as_brst bs_brst
bs_tumr = -log(tumr_strct(2).mus/tumr_strct(1).mus)/log(tumr_strct(2).lambda/tumr_strct(1).lambda); as_tumr = tumr_strct(1).mus/tumr_strct(1).lambda^(-bs_tumr); mus_tumr = as_tumr*(lambda_in)^(-bs_tumr); % clearvars as_tumr bs_tumr
[x,iv,~] = unique(brst_db.Var1); v = brst_db.Var2(iv)*10; mua_brst = interp1(x,v,lambda_in,'pchip','extrap'); clearvars x v iv
[x,iv,~] = unique(tumr_db.Var1); v = tumr_db.Var2(iv)*10; mua_tumr = interp1(x,v,lambda_in,'pchip','extrap'); clearvars x v iv
if do_fig
    subplot(1,2,1)
    plot(brst_db.Var1,brst_db.Var2*10,'DisplayName','breast: model','LineStyle','-','Color','b','LineWidth',1.5), hold on
    plot([brst_strct(:).lambda],[brst_strct(:).mua],'Marker','s','MarkerFaceColor','b','MarkerEdgeColor','k','DisplayName','breast: data','MarkerSize',16,'LineStyle','none'), hold on
    plot(lambda_in,mua_brst,'Marker','h','MarkerFaceColor','b','MarkerEdgeColor','k','DisplayName','breast: output','MarkerSize',16,'LineStyle','none'), hold on
    plot(tumr_db.Var1,tumr_db.Var2*10,'DisplayName','tumor: model', 'LineStyle','-','Color','r','LineWidth',1.5), hold on
    plot([tumr_strct(:).lambda],[tumr_strct(:).mua],'Marker','s','MarkerFaceColor','r','MarkerEdgeColor','k','DisplayName','tumor: data','MarkerSize',16,'LineStyle','none'), hold on
    plot(lambda_in,mua_tumr,'Marker','h','MarkerFaceColor','r','MarkerEdgeColor','k','DisplayName','tumor: output','MarkerSize',16,'LineStyle','none'), hold on
    xlabel('wavelength (\lambda, nm)'), ylabel('\mu_a (cm^{-1})'), title('\mu_a vs. wavelength')
    axis tight, axis square, set(gca,'fontsize',18)
    legend('show','Location','northwest')
    hold off

    subplot(1,2,2)
    plot(brst_db.Var1,as_brst.*(brst_db.Var1.^(-bs_brst)),'DisplayName','breast: model','LineStyle','-','Color','b','LineWidth',1.5), hold on
    plot([brst_strct(:).lambda],[brst_strct(:).mus],'Marker','s','MarkerFaceColor','b','MarkerEdgeColor','k','DisplayName','breast: data','MarkerSize',16,'LineStyle','none'), hold on
    plot(lambda_in,mus_brst,'Marker','h','MarkerFaceColor','b','MarkerEdgeColor','k','DisplayName','breast: output','MarkerSize',16,'LineStyle','none'), hold on
    plot(tumr_db.Var1,as_tumr.*(tumr_db.Var1.^(-bs_tumr)),'DisplayName','tumor: model', 'LineStyle','-','Color','r','LineWidth',1.5), hold on
    plot([tumr_strct(:).lambda],[tumr_strct(:).mus],'Marker','s','MarkerFaceColor','r','MarkerEdgeColor','k','DisplayName','tumor: data','MarkerSize',16,'LineStyle','none'), hold on
    plot(lambda_in,mus_tumr,'Marker','h','MarkerFaceColor','r','MarkerEdgeColor','k','DisplayName','tumor: output','MarkerSize',16,'LineStyle','none'), hold on
    xlabel('wavelength (\lambda, nm)'), ylabel('\mu''_s (cm^{-1})'), title('\mu''_s vs. wavelength')
    axis tight, axis square, set(gca,'fontsize',18)
    legend('show','Location','northeast')
    hold off
end
end
