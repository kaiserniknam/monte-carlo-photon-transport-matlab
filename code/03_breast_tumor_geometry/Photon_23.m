function [] = Photon_23 ()
% Repository group: 03_breast_tumor_geometry
% Version role: primary simulation or study entry point
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: s vs. d for breast tissue (trace photon path inside tumor)
% ref. #1: Assessment of the size, position, and optical properties of breast tumors in vivo by noninvasive optical methods
% ref. #2: Quantitative Absorption and Scattering Spectra in Thick Tissues Using Broadband Diffuse Optical Spectroscopy
% ref. #3: Tumor location of the central and nipple portion is associated with impaired survival for women with breast cancer
% ref. #4: Broadband Optical Mammography: Chromophore Concentration and Hemoglobin Saturation Contrast in Breast Cancer

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

for idx = length(brst_strct)
    do_count = do_count + 1;
    mua_brst = brst_strct(idx).mua; g_brst = brst_strct(idx).g; n_brst = brst_strct(idx).n; mus_brst = brst_strct(idx).mus/(1-g_brst);
    mua_tumr = tumr_strct(idx).mua; g_tumr = tumr_strct(idx).g; n_tumr = tumr_strct(idx).n; mus_tumr = tumr_strct(idx).mus/(1-g_tumr); wvlnt = tumr_strct(idx).lambda;

    [d_tumor, s_tumor, w_tumor, c_tumor, no_of_surf_photons_tumor, no_of_dpth_photons_tumor] = do_simulation ( ...
        [mua_brst,mus_brst,g_brst,n_brst], ...
        [mua_tumr,mus_tumr,g_tumr,n_tumr], ...
        [Lx,Ly,Lz],[cx,cy,cz],[rx,ry,rz],z_air,wvlnt,nPhotonsReq,nExamplePaths);

    [d_brest, s_brest, w_brest, c_brest, no_of_surf_photons_brest, no_of_dpth_photons_brest] = do_simulation ( ...
        [mua_brst,mus_brst,g_brst,n_brst], ...
        [mua_brst,mus_brst,g_brst,n_brst], ...
        [Lx,Ly,Lz],[cx,cy,cz],[rx,ry,rz],z_air,wvlnt,nPhotonsReq,nExamplePaths);

    d_brest_pass = d_brest(c_brest==1); d_brest_nops = d_brest(c_brest==0);
    s_brest_pass = s_brest(c_brest==1); s_brest_nops = s_brest(c_brest==0);
    w_brest_pass = w_brest(c_brest==1); w_brest_nops = w_brest(c_brest==0);
    d_tumor_pass = d_tumor(c_tumor==1); d_tumor_nops = d_tumor(c_tumor==0);
    s_tumor_pass = s_tumor(c_tumor==1); s_tumor_nops = s_tumor(c_tumor==0);
    w_tumor_pass = w_tumor(c_tumor==1); w_tumor_nops = w_tumor(c_tumor==0);

    figure(1), subplot(1,2,1), the_seed = 0:5:150;
    [N,edges] = histcounts(s_brest_pass,the_seed,'Normalization','probability');
    plot(1/2*(edges(1:end-1)+edges(2:end-0)),(N),'LineWidth',2,'DisplayName','no tumor'), hold on, clearvars N edges
    [N,edges] = histcounts(s_tumor_pass,the_seed,'Normalization','probability');
    plot(1/2*(edges(1:end-1)+edges(2:end-0)),(N),'LineWidth',2,'DisplayName','with tumor'), hold on, clearvars N edges
    xlabel('s (cm)'), ylabel('freq (unitless)'), hold off
    title([{'hist. of photon pathlength when passing tumor loc.'},{['\mu_{no-tumor} = ',num2str(mean(s_brest_pass),'%.2f'),' vs . \mu_{tumor} = ',num2str(mean(s_tumor_pass),'%.2f'),' (p-val = ',num2str(ranksum(s_brest_pass,s_tumor_pass)),')']}])
    set(gca,'fontsize',16), axis tight, grid on, legend('show','Location','northeast')
    figure(1), subplot(1,2,2), the_seed = 0:5:150;
    [N,edges] = histcounts(s_brest_nops,the_seed,'Normalization','probability');
    plot(1/2*(edges(1:end-1)+edges(2:end-0)),(N),'LineWidth',2,'DisplayName','no tumor'), hold on, clearvars N edges
    [N,edges] = histcounts(s_tumor_nops,the_seed,'Normalization','probability');
    plot(1/2*(edges(1:end-1)+edges(2:end-0)),(N),'LineWidth',2,'DisplayName','with tumor'), hold on, clearvars N edges
    xlabel('s (cm)'), ylabel('freq (unitless)'), hold off
    title([{'hist. of photon pathlength when passing no tumor loc.'},{['\mu_{no-tumor} = ',num2str(mean(s_brest_nops),'%.2f'),' vs . \mu_{tumor} = ',num2str(mean(s_tumor_nops),'%.2f'),' (p-val = ',num2str(ranksum(s_brest_nops,s_tumor_nops)),')']}])
    set(gca,'fontsize',16), axis tight, grid on, legend('show','Location','northeast')
    clearvars s_tumor_nops s_brest_nops s_tumor_pass s_brest_pass

    figure(2), the_seed = 0:0.05:1;
    [N,edges] = histcounts(w_brest_pass,the_seed,'Normalization','probability');
    plot(1/2*(edges(1:end-1)+edges(2:end-0)),(N),'LineWidth',2,'DisplayName','no tumor'), hold on, clearvars N edges
    [N,edges] = histcounts(w_tumor_pass,the_seed,'Normalization','probability');
    plot(1/2*(edges(1:end-1)+edges(2:end-0)),(N),'LineWidth',2,'DisplayName','with tumor'), hold on, clearvars N edges
    xlabel('w (%)'), ylabel('freq (unitless)'), hold off
    title([{'histogram of attenuation when passing tumor loc.'},{['\mu_{no-tumor} = ',num2str(mean(w_brest_pass),'%.2f'),' vs . \mu_{tumor} = ',num2str(mean(w_tumor_pass),'%.2f'),' (p-val = ',num2str(ranksum(w_brest_pass,w_tumor_pass)),')']}])
    set(gca,'fontsize',24), axis tight, grid on, legend('show','Location','northwest')
    % figure(2), subplot(1,2,2), the_seed = 0:0.05:1;
    % [N,edges] = histcounts(w_brest_nops,the_seed,'Normalization','probability');
    % plot(1/2*(edges(1:end-1)+edges(2:end-0)),(N),'LineWidth',2,'DisplayName','no tumor'), hold on, clearvars N edges
    % [N,edges] = histcounts(w_tumor_nops,the_seed,'Normalization','probability');
    % plot(1/2*(edges(1:end-1)+edges(2:end-0)),(N),'LineWidth',2,'DisplayName','with tumor'), hold on, clearvars N edges
    % xlabel('w (%)'), ylabel('freq (unitless)'), hold off
    % title([{'hist. of attenuation when passing no tumor loc.'},{['\mu_{no-tumor} = ',num2str(mean(w_brest_nops),'%.2f'),' vs . \mu_{tumor} = ',num2str(mean(w_tumor_nops),'%.2f'),' (p-val = ',num2str(ranksum(w_brest_nops,w_tumor_nops)),')']}])
    % set(gca,'fontsize',16), axis tight, grid on, legend('show','Location','northwest')

    figure(3), subplot(1,2,1), N_bin = 150;
    d_t = d_brest_pass; w_t = w_brest_pass; n_t = no_of_surf_photons_brest+no_of_dpth_photons_brest;
    d_edges = linspace(min(d_t),max(d_t),N_bin+1); % # d bins
    d_c = 1/2*(d_edges(1:end-1)+d_edges(2:end-0))';
    [f_d,~,ind_d] = histcounts(d_t,d_edges); f_d = f_d./sum(f_d); clearvars d_edges
    w_mean = accumarray(ind_d,w_t,[],@sum); w_mean = w_mean/n_t;
    figure(3), plot(d_c,-log10(w_mean),'LineWidth',2,'DisplayName','no tumor'), hold on
    clearvars d_t w_t n_t
    d_t = d_tumor_pass; w_t = w_tumor_pass; n_t = no_of_surf_photons_brest+no_of_dpth_photons_brest;
    d_edges = linspace(min(d_t),max(d_t),N_bin+1); % # d bins
    d_c = 1/2*(d_edges(1:end-1)+d_edges(2:end-0))';
    [f_d,~,ind_d] = histcounts(d_t,d_edges); f_d = f_d./sum(f_d); clearvars d_edges
    w_mean = accumarray(ind_d,w_t,[],@sum); w_mean = w_mean/n_t;
    figure(3), plot(d_c,-log10(w_mean),'LineWidth',2,'DisplayName','with tumor'), hold on
    clearvars d_t w_t n_t
    xlabel('d (cm)'), ylabel('OD (unitless)'), title('OD vs. d')
    title('OD vs. d when passing tumor loc.')
    set(gca,'fontsize',24), axis tight, grid on, legend('show','Location','southeast'), hold on

    figure(3), subplot(1,2,2), N_bin = 150;
    d_t = d_brest_nops; w_t = w_brest_nops; n_t = no_of_surf_photons_brest+no_of_dpth_photons_brest;
    d_edges = linspace(min(d_t),max(d_t),N_bin+1); % # d bins
    d_c = 1/2*(d_edges(1:end-1)+d_edges(2:end-0))';
    [f_d,~,ind_d] = histcounts(d_t,d_edges); f_d = f_d./sum(f_d); clearvars d_edges
    w_mean = accumarray(ind_d,w_t,[],@sum); w_mean = w_mean/n_t;
    figure(3), plot(d_c,-log10(w_mean),'LineWidth',2,'DisplayName','no tumor'), hold on
    clearvars d_t w_t n_t
    d_t = d_tumor_nops; w_t = w_tumor_nops; n_t = no_of_surf_photons_brest+no_of_dpth_photons_brest;
    d_edges = linspace(min(d_t),max(d_t),N_bin+1); % # d bins
    d_c = 1/2*(d_edges(1:end-1)+d_edges(2:end-0))';
    [f_d,~,ind_d] = histcounts(d_t,d_edges); f_d = f_d./sum(f_d); clearvars d_edges
    w_mean = accumarray(ind_d,w_t,[],@sum); w_mean = w_mean/n_t;
    figure(3), plot(d_c,-log10(w_mean),'LineWidth',2,'DisplayName','with tumor'), hold on
    clearvars d_t w_t n_t
    xlabel('d (cm)'), ylabel('OD (unitless)'), title('OD vs. d')
    title('OD vs. d when passing no tumor loc.')
    set(gca,'fontsize',24), axis tight, grid on, legend('show','Location','southeast'), hold on

    clearvars d s a w d_c s_c s_mean w_mean i_mean s_var f_d f_s no_of_dpth_photons no_of_surf_photons
    clearvars mua_brst g_brst n_brst mus_brst
    clearvars mua_tumr g_tumr n_tumr mus_tumr
end
end

function [d, s, w, c, no_of_surf_photons, no_of_dpth_photons] = do_simulation (opt_bckg,opt_sgnl,cmd_size,sig_pos,sig_size,zSurface,wvlngth,nPhotonsReq,nExamplePaths)
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
clearvars opt_bckg opt_sgnl cmd_size dl

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
    t_scatter_stat = nan(length(fnsh_pnts),3+1+1+1); clearvars spratrs
    for idx = 1:length(strt_pnts)
        Photon_Path = t_model.MC.examplePaths(:,strt_pnts(idx)+1:fnsh_pnts(idx)-1);
        if isnan(sum(sum(Photon_Path)))
            keyboard
        else
        inside_sig = sum(sqrt( ...
            ((Photon_Path(1,:)-sig_pos(1)         ).^2)./((sig_size(1)).^2) + ...
            ((Photon_Path(2,:)-sig_pos(2)         ).^2)./((sig_size(2)).^2) + ...
            ((Photon_Path(3,:)-sig_pos(3)-zSurface).^2)./((sig_size(3)).^2)) <= 1);
            if Photon_Path(3,end) <= 0.0
                t_scatter_stat(idx,:) = [...
                    Photon_Path(1,end), ...
                    Photon_Path(2,end), ...
                    Photon_Path(3,end), ...
                    sum(sqrt(sum((Photon_Path(1:3,2:end-0)-Photon_Path(1:3,1:end-1)).^2,1))), ...
                    Photon_Path(4,end), ...
                    inside_sig>0, ...
                    ];
                if inside_sig>0
                    % plot3(Photon_Path(1,:),Photon_Path(2,:),-Photon_Path(3,:),'g'), hold on
                else
                    % plot3(Photon_Path(1,:),Photon_Path(2,:),-Photon_Path(3,:),'r'), hold on
                end
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
c = scatter_stat(:,6);                                   % a -> 1 if inside the signal
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
