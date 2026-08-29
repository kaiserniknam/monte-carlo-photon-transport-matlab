function [] = Photon_64 ()
% Repository group: 07_dpf_parameter_studies
% Version role: primary simulation or study entry point
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: w & s vs. d over mu_s & mu_a: more photons, more combinations
% study the effect of d on s, s_max, and s_min (why s vs. d is limited?) for my TiO2 experiment
% based on 37

clc
close all

z_air = 0.0; % the thickness of air layer
nPhotonsReq = 1e5; % number of requested photon
nExamplePaths = 500000; % number of example photon paths
% Optical & size properties
Lx = 29.1; Ly = 29.1; Lz = 05.0; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
beam_phi = 0;     % polar angle of beam in cm
beam_tht = 0;     % azimuthal angle of beam in cm
beam_X = 0.0;     % X deviation of beam in cm
beam_Y = 0.0;     % Y deviation of beam in cm

n = 1.33; g = 0.93;
set_of_mua = 0.0275;
set_of_mus = 55.0;
dlta_d = 0.17;

for i_a = 1:length(set_of_mua)
    for i_s = 1:length(set_of_mus)
        mua = set_of_mua(i_a);
        mus = set_of_mus(i_s);
        the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_64_mua_',sprintf('%.2f',set_of_mua(i_a)),'_mus_',sprintf('%.2f',set_of_mus(i_s)),'.mat' ];
        if ~exist(the_filename,'file')
            [x, y, z, d, s, w, c, a, no_of_photons, M_raw,SetExamplePaths] = do_simulation ( ...
                [mua,mus,g   ,n  ], ...
                [nan,nan,nan ,nan], ...
                [Lx,Ly,Lz],[nan,nan,nan],[nan,nan,nan],z_air,001,nPhotonsReq,nExamplePaths,beam_X,beam_Y,beam_phi,beam_tht);

            save  (the_filename, ...
                'x', 'y', 'z', 'd', 's', 'w', 'c', 'a', 'no_of_photons', 'M_raw', 'SetExamplePaths', '-v7.3')

            clearvars x y z d s w c a no_of_photons M_raw SetExamplePaths
            clearvars mua mus
        else
            t_db = load(the_filename);

            % removing noisy points
            u = round(-log(t_db.w)./t_db.s,4); u_unique = unique(u); freq = nan(size(u_unique));
            for i_f = 1:length(u_unique), freq(i_f) = sum(u==u_unique(i_f)); end
            [~,i_fx] = max(freq);
            t_db.x = t_db.x(u==u_unique(i_fx));
            t_db.y = t_db.y(u==u_unique(i_fx));
            t_db.z = t_db.z(u==u_unique(i_fx));
            t_db.d = t_db.d(u==u_unique(i_fx));
            t_db.s = t_db.s(u==u_unique(i_fx));
            t_db.w = t_db.w(u==u_unique(i_fx));
            t_db.c = t_db.c(u==u_unique(i_fx));
            t_db.a = t_db.a(u==u_unique(i_fx));
            t_db.SetExamplePaths = t_db.SetExamplePaths(u==u_unique(i_fx));

            t_db.no_of_photons = t_db.no_of_photons-sum(u~=u_unique(i_fx));
            clearvars u u_unique freq i_f i_fx the_filename

            % 1-D sorting
            d_diff_edges = 0:dlta_d:sqrt((Lx/2).^2+(Ly/2).^2); clearvars Lx Ly Lz
            [~,~,ind_diff] = histcounts(sqrt(t_db.x(t_db.c==0).^2+t_db.y(t_db.c==0).^2+t_db.z(t_db.c==0).^2),d_diff_edges); % d_diffuse bins
            clearvars d_trns_edges
            clearvars beam_phi beam_tht beam_X beam_Y

            % I vs. s (scatterplot)
            TheColor = 'b'; TheOutFun = @(x)(-log(x)); TheLineStyle = 'none'; TheMarker = 'o'; x_label = 's (cm)'; y_label = 'OD (a.u.)';
            i_subplot = 1; TheCode = 0; TheTitle = 'OD vs. s';
            x_temp = t_db.s(t_db.c==TheCode);
            y_temp = t_db.w(t_db.c==TheCode); y_temp = TheOutFun(y_temp);
            subplot(1,2,i_subplot)
            plot(x_temp,y_temp,'Color',TheColor,'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','w','MarkerFaceColor',TheColor,'MarkerSize',8), hold on
            xlim([0 170]), xlabel(x_label), ylabel(y_label), title(TheTitle), axis([0 170 0 5])
            set(gca,'fontsize',16), axis square, grid on, hold on
            clearvars i_subplot index_in x_temp y_temp TheCode TheTitle mdl
            clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label TheLegend TheColor

            % I vs. d (scatterplot)
            TheColor = 'b'; TheLineStyle = '-'; TheMarker = 'none'; x_label = 'd (cm)'; y_label = 'OD (a.u.)';
            fun_x = @mean; fun_y = @sum; TheOutFun = @(x)(-log(x)); TheLegend = 'simulation';
            i_fig = 1; i_subplot = 2; TheCode = 0; TheTitle = 'OD vs. d'; index_in = ind_diff;
            x_temp = t_db.d(t_db.c==TheCode);
            y_temp = t_db.w(t_db.c==TheCode);
            [~,~] = make_a_subplot(i_fig,i_subplot,index_in,x_temp,y_temp,fun_x,fun_y,TheColor,TheLegend,TheOutFun,TheLineStyle,TheMarker,true,dlta_d,t_db.no_of_photons);
            xlim([0 8]), xlabel(x_label), ylabel(y_label), title(TheTitle), axis([0 8 0 15])
            set(gca,'fontsize',16), axis square, grid on, hold on
            clearvars i_subplot index_in x_temp y_temp TheCode TheTitle mdl
            clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label TheLegend TheColor
            % ---- read
            fname = "/home/kaiser/Dropbox/KSR/University of Houston/Experiment/2.NIRs Project/008 My DFP/code/db/OD_grid_2025-12-29-11-40-40-750nm-OD-1gL.txt";     % 1.00 g/L
            T = readtable(fname, "FileType","text", "Delimiter","\t"); clearvars fname
            db_1_00 = T.Variables; clearvars T fname % assumes fixed column order
            [x_bind,y_bind] = get_data(db_1_00);
            plot(x_bind,y_bind,'Color','r','DisplayName','data','Marker','o','LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',10,'LineWidth',1.5), hold on
            clearvars db_1_00 x_bind y_bind fname



            % banana-shaped path from a light source to a photodetector
            idx_involved = 4:3:25; N_d = 51;
            set_of_paths = nan(length(idx_involved),N_d-1,3);
            for i_d = idx_involved % 1:length(d_diff_edges)-1
                set_of_dist = linspace(0,d_diff_edges(i_d),N_d);
                idx_in = ...
                    d_diff_edges(i_d-1)<=[t_db.SetExamplePaths(:).PathDistnc]&...
                                         [t_db.SetExamplePaths(:).PathDistnc]<d_diff_edges(i_d+1)&t_db.c.'==0;
                SetExamplePaths = t_db.SetExamplePaths(idx_in); clearvars idx_in
                the_collected_coords = nan(length(SetExamplePaths),length(set_of_dist)-1,3);
                for i_path = 1:length(SetExamplePaths)
                    x_coords = SetExamplePaths(i_path).PathCoords(1,:);
                    y_coords = SetExamplePaths(i_path).PathCoords(2,:);
                    z_coords = SetExamplePaths(i_path).PathCoords(3,:);
                    [x_coords,y_coords,z_coords] = align_X_direction(x_coords,y_coords,z_coords);
                    [~,~,idx_in] = histcounts(x_coords,set_of_dist); % d_diffuse bins
                    x_coords(idx_in==0) = [];
                    y_coords(idx_in==0) = [];
                    z_coords(idx_in==0) = [];
                    idx_in(idx_in==0) = [];
                    tmp_vector = [accumarray(idx_in.',x_coords.',[],@mean,nan) accumarray(idx_in.',y_coords.',[],@mean,nan) accumarray(idx_in.',z_coords.',[],@mean,nan)];
                    the_collected_coords(i_path,1:size(tmp_vector,1),:) = tmp_vector;
                    clearvars x_coords y_coords z_coords idx_in tmp_vector
                end
                set_of_paths(i_d,:,:) = nanmean(the_collected_coords(:,:,:),1);
                clearvars set_of_dist SetExamplePaths the_collected_coords i_path
            end
            clearvars i_d N_d set_of_dist
            set_of_depth = nan(size(set_of_paths,1),2);
            for i_d = idx_involved
                figure(2), subplot(1,2,1)
                plot3([squeeze(set_of_paths(i_d,:,1))],[squeeze(set_of_paths(i_d,:,2))],[-squeeze(set_of_paths(i_d,:,3))],'DisplayName',['d = ',num2str(d_diff_edges(i_d))],'LineWidth',2), hold on
                set_of_depth(i_d,1) = round(d_diff_edges(i_d),1);
                set_of_depth(i_d,2) = max(squeeze(set_of_paths(i_d,:,3)));
            end
            xlabel('x (cm)'), ylabel('y (cm)'), zlabel('z (cm)'), set(gca,'fontsize',24)
            title('banana-shape pathways')
            axis equal, grid on, view([00 00])

            figure(2), subplot(1,2,2)
            set_of_depth(isnan(set_of_depth(:,1))&isnan(set_of_depth(:,2)),:) = [];
            plot(set_of_depth(:,1),-set_of_depth(:,2),'DisplayName','none','LineWidth',2), hold on
            xlabel('d (cm)'), ylabel('z (cm)'), set(gca,'fontsize',24)
            title('effective depth vs. SDS')
            grid on, axis square
        end
        clearvars t_db the_filename mua mus ind_diff ind_trns
    end
end
end

function [x, y, z, d, s, w, c, a, no_of_photons, M_raw, SetExamplePaths] = do_simulation (opt_bckg,opt_sgnl,cmd_size,sig_pos,sig_size,zSurface,wvlngth,nPhotonsReq,nExamplePaths,beam_X,beam_Y,beam_phi,beam_tht)
%% Geometry definition
SetExamplePaths = struct();
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
counter = 0;
while size(scatter_stat,1) < nExamplePaths
    t_model = runMonteCarlo(model);
    M_raw = t_model.G.M_raw;
    % t_model = plot(t_model,'MC');
    % do calc
    spratrs = find(isnan(t_model.MC.examplePaths(1,:)));
    strt_pnts = spratrs(2:2:end); strt_pnts = strt_pnts(1:end-1);
    fnsh_pnts = spratrs(1:2:end); fnsh_pnts = fnsh_pnts(2:end-0);
    t_scatter_stat = nan(length(fnsh_pnts),3+1+1+3+1); clearvars spratrs
    for idx = 1:length(strt_pnts)
        Photon_Path = t_model.MC.examplePaths(:,strt_pnts(idx)+1:fnsh_pnts(idx)-1);
        if Photon_Path(3,end) <= 0.0
            counter = counter + 1;
            % plot3(Photon_Path(1,:),Photon_Path(2,:),-Photon_Path(3,:),'r'), hold on
            SetExamplePaths(counter).PathCoords = Photon_Path;
            SetExamplePaths(counter).PathLength = sum(sqrt(sum((Photon_Path(1:3,2:end-0)-Photon_Path(1:3,1:end-1)).^2,1)));
            SetExamplePaths(counter).PathDistnc = sqrt(Photon_Path(1,end).^2+Photon_Path(2,end).^2+Photon_Path(3,end).^2);
            SetExamplePaths(counter).PathWeight = Photon_Path(4,end);
            t_scatter_stat(idx,:) = [...
                Photon_Path(1,end), ...
                Photon_Path(2,end), ...
                Photon_Path(3,end), ...
                sum(sqrt(sum((Photon_Path(1:3,2:end-0)-Photon_Path(1:3,1:end-1)).^2,1))), ...
                Photon_Path(4,end), ...
                (Photon_Path(1,end)-Photon_Path(1,end-1)), ...
                (Photon_Path(2,end)-Photon_Path(2,end-1)), ...
                (Photon_Path(3,end)-Photon_Path(3,end-1)), ...
                (0), ...
                ];
        else
            % plot3(Photon_Path(1,:),Photon_Path(2,:),-Photon_Path(3,:),'k'), hold on
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
x = scatter_stat(:,1);                                                 % x -> source to detector distance
y = scatter_stat(:,2);                                                 % y -> source to detector distance
z = scatter_stat(:,3);                                                 % y -> source to detector distance
d = sqrt((scatter_stat(:,1)-beam_X).^2+(scatter_stat(:,2)-beam_Y).^2); % d -> source to detector distance
s = scatter_stat(:,4);                                                 % s -> true distance
w = scatter_stat(:,5);                                                 % w -> intensity
a = scatter_stat(:,8)./sqrt(scatter_stat(:,6).^2 + scatter_stat(:,7).^2 + scatter_stat(:,8).^2); % a -> angle
c = scatter_stat(:,9);                                                 % c -> code: diffuse reflectance, 0; transmittive, 1
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

function [x_bind,y_bind] = make_a_subplot(i_fig,i_subplot,index_in,x_temp,y_temp,fun_x,fun_y,TheColor,TheLegend,TheOutFun,TheLineStyle,TheMarker,blnNorm,dlta_d,no_of_photons)
x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
y_bind = accumarray(index_in,y_temp,[],fun_y,nan);
if blnNorm
    y_bind = (y_bind./(pi.*dlta_d.*dlta_d + 2.*pi.*x_bind.*dlta_d)) ./ ...
        (no_of_photons./(pi.*dlta_d.*dlta_d));
    figure(i_fig)
    subplot(1,2,i_subplot)
    plot(x_bind,TheOutFun(y_bind),'Color',TheColor,'DisplayName',TheLegend,'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','w','MarkerFaceColor',TheColor,'MarkerSize',8,'LineWidth',1.5), hold on
end
end
function [out,Ns] = get_average(SetExamplePaths,d_edges)
N_bin = length(d_edges)-1;
out = nan(3,N_bin);
for i_d = 1:N_bin
    the_collected_coords = [];
    for idx = 1:length(SetExamplePaths)
        x_coords = SetExamplePaths(idx).PathCoords(1,:);
        y_coords = SetExamplePaths(idx).PathCoords(2,:);
        z_coords = SetExamplePaths(idx).PathCoords(3,:);
        [x_coords,y_coords,z_coords] = align_X_direction(x_coords,y_coords,z_coords);
        idx_in = (d(i_d)<=x_coords)&(x_coords<=d(i_d+1));
        the_collected_coords = [the_collected_coords , [mean(x_coords(idx_in)) ; mean(y_coords(idx_in)) ; mean(z_coords(idx_in))]];
        clearvars idx_in x_coords y_coords z_coords
    end
    if ~isempty(mean(the_collected_coords,2))
        out(:,i_d) = nanmean(the_collected_coords,2);
    end
    % if i_d==N_bin, keyboard; end
    clearvars idx the_collected_coords
end
clearvars d_strt d_fnsh d i_d
out(:,isnan(out(1,:))) = [];
Ns = length(SetExamplePaths);
% better visualization
% if ~isempty(out)
%     f_l = nan(2,3,size(out,2));
%     for idx = 1:length(SetExamplePaths)
%         x_coords = SetExamplePaths(idx).PathCoords(1,:);
%         y_coords = SetExamplePaths(idx).PathCoords(2,:);
%         z_coords = SetExamplePaths(idx).PathCoords(3,:);
%         [x_coords,y_coords,z_coords] = align_X_direction(x_coords,y_coords,z_coords);
%         f_l(1,:,idx) = [x_coords(1),  y_coords(1),  z_coords(1)  ];
%         f_l(2,:,idx) = [x_coords(end),y_coords(end),z_coords(end)];
%         clearvars x_coords y_coords z_coords
%     end
%     new_out = nan(3,size(out,2)+2);
%     new_out(:,2:end-1) = out; clearvars out
%     new_out(:,1) =   squeeze(mean(f_l(1,  :,:),3));
%     new_out(:,end) = squeeze(mean(f_l(end,:,:),3));
%     out = new_out;
% end
end
function [x_coords,y_coords,z_coords] = align_X_direction(x_coords,y_coords,z_coords)
tht = atan2(y_coords(end)-y_coords(1),x_coords(end)-x_coords(1));
% perform the rotation
R = [cos(tht), sin(tht);-sin(tht), cos(tht)];
rot_pnts = R*[x_coords; y_coords];
x_coords = rot_pnts(1,:);
y_coords = rot_pnts(2,:);
end
function [x_bind,y_bind] = get_data(db_1_00)
addpath('/home/kaiser/Dropbox/KSR/University of Houston/Experiment/2.NIRs Project/008 My DFP/code/Characteristics/')

% ---- parameters
R_L    = 200e3;
R_LED  = 150;
lambda = 750;

% IMPORTANT: set this to your actual spatial pitch
pitch_cm = 1.0;   % <-- CHANGE if each grid step is not 1 cm

% Columns assumed (based on your code):
% 1=time_ms, 3=voltage_V, 4=i_row, 5=i_col, 6=recording_flag (1=ON)

% ---- extract intensity map
N = 13;
% Preallocate max possible vector length: N*N points
The_Vector = nan(N*N, 1+1);
k = 0;
for i_row = 1:N
    for i_col = 1:N
        sig = db_1_00(db_1_00(:,4)==i_row & db_1_00(:,5)==i_col & db_1_00(:,6)==1, 3);
        sig(sig<eps) = nan;
        if isempty(sig) || all(isnan(sig))
            Io_1 = nan;
        else
            % Vmean = mean(sig, "omitnan");
            Vmean = mean(sig);
            [~, Io_1] = do_FDS100(Vmean, R_L, lambda);
        end
        dist_steps = sqrt((i_row-7).^2 + (i_col-7).^2);
        dist_cm    = dist_steps * pitch_cm;
        clearvars dist_steps sig Vmean

        k = k + 1;
        The_Vector(k,:) = [dist_cm, Io_1];
        clearvars dist_cm Io_25 Io_1 Io_4
    end
end
The_Vector = The_Vector(1:k,:); % trim
clearvars N k i_row i_col pitch_cm

% ---- source intensity I0
if lambda == 750
    [~, I0, ~] = do_LED750L(R_LED);
else
    [~, I0, ~] = do_LED850LN(R_LED);
end
The_Vector(:,2) = The_Vector(:,2) ./ I0;
clearvars I0 R_L R_LED lambda

OD = The_Vector(:,2:end);
OD = -log(OD);

x_bind = The_Vector(:,1);
y_bind = OD;
end
