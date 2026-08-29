function [] = Photon_45 ()
% Repository group: 05_digital_breast_phantoms
% Version role: primary simulation or study entry point
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: simulating phantoms prepared by Diego
% Similar to Photon_44: but in better intervals

clc
close all

wvlnt = 800; % sample wavelength (nm)
nPhotonsReq = 1e5;  % number of requested photon
nExamplePaths = 1000000; % number of example photon paths
% Optical properties of breast and tumor [1-4]
mua_glnd = 0.55; mus_glnd = 332.7; g_glnd = 0.965; n_glnd = 1.4;
mua_adps = 0.82; mus_adps = 313.8; g_adps = 0.976; n_adps = 1.4;
set_of_percnt = [25,50,75]; % percent of adipose
set_of_versns = [0,1,2]; % version of each density
set_of_HCT_Cn = 0:2.5:40; % HCT (in %)
set_of_beam_X = [-1.0,0.0,+1.0]; % beam_X location (in cm)
set_of_Tumr_Z = [1.0,1.25,1.50,1.75]; % Tumor_Z location (in cm)
Nz = 1020; Ny = 257; Nx = 323; dl = 0.2/10; r_lowres = 0.5; dl = dl./r_lowres;
the_Tumr_X = 0; the_Tumr_Y = 1; the_Tumr_r = 0.75; % Tumor location and size all in cm
                the_beam_Y = 1; the_beam_Z = 0; % Source location in cm

for i_percent = 1:length(set_of_percnt)
    for i_versns = 1:length(set_of_versns)
        for i_HCT_Cn = 1:length(set_of_HCT_Cn)
            for i_beam_X = 1:length(set_of_beam_X)
                for i_Tumr_Z = 1:length(set_of_Tumr_Z)
                    the_percnt = set_of_percnt(i_percent);
                    the_verson = set_of_versns(i_versns);
                    the_HCT_Cn = set_of_HCT_Cn(i_HCT_Cn);
                    the_beam_X = set_of_beam_X(i_beam_X);
                    the_Tumr_Z = set_of_Tumr_Z(i_Tumr_Z);

                    the_path = ['/home/kaiser/Phantoms/Diego/Phantoms/CupB_Pd',num2str(the_percnt),'_',num2str(the_verson),'_y_-z_x.obj'];
                    fid = fopen(the_path, 'r');
                    data = fread(fid,'float32');
                    fclose(fid);
                    data(        data<=0.0) = 0.0; % Air
                    data(0.0<data&data<0.9) = 0.1; % Fat (Adipose Tissue)
                    data(0.9<=data        ) = 0.2; % Fibroglandular Tissue
                    data = round(data*10);
                    % Reshape the data into a 3D array
                    TheImage = reshape(data, [Nz, Ny, Nx]);
                    TheImage = round(imresize3(TheImage,r_lowres));
                    TheImage = permute(TheImage,[1,3,2]); % TheImage = TheImage(:,:,end:-1:1);
                    clearvars data the_path
                    % put a tumor if necessary
                    if the_HCT_Cn>0
                        [X,Y,Z] = ndgrid( ...
                            ((1:size(TheImage,1))-size(TheImage,1)/2).*dl, ...
                            ((1:size(TheImage,2))-size(TheImage,2)/2).*dl, ...
                            ((1:size(TheImage,3))                   ).*dl);
                        TheImage( ...
                            (X-the_Tumr_X).^2+ ...
                            (Y-the_Tumr_Y).^2+ ...
                            (Z-the_Tumr_Z).^2<=the_Tumr_r) = 3;
                        clearvars X Y Z
                    end
                    % calculate optical properties of tumor
                    [mua_tumr,mus_tumr,g_tumr,~] = calc_muas_based_HCT(the_HCT_Cn,wvlnt);
                    n_tumr = refractive_index_water(wvlnt).*((Specific_Refractive_Increment_beta(wvlnt).*the_HCT_Cn./3)+1);

                    [p_in, i_in, m_in, w_in, p_ot, i_ot, m_ot, w_ot, s, dl, no_of_photons, M_raw] = do_simulation ( ...
                        [mua_glnd,mus_glnd,g_glnd,n_glnd], ...
                        [mua_adps,mus_adps,g_adps,n_adps], ...
                        [mua_tumr,mus_tumr,g_tumr,n_tumr], ...
                        TheImage,dl,wvlnt,nPhotonsReq,nExamplePaths,the_beam_X,the_beam_Y,the_beam_Z);

                    save  (['Photon_45_Pd',num2str(the_percnt),'_',num2str(the_verson),'_',num2str(the_HCT_Cn),'_',num2str(the_beam_X),'_',num2str(the_Tumr_Z),'_y_-z_x.mat' ], ...
                        'p_in', 'i_in', 'm_in', 'w_in', ...
                        'p_ot', 'i_ot', 'm_ot', 'w_ot', ...
                        'mua_glnd', 'mus_glnd', 'g_glnd', 'n_glnd', ...
                        'mua_adps', 'mus_adps', 'g_adps', 'n_adps', ...
                        'mua_tumr', 'mus_tumr', 'g_tumr', 'n_tumr', ...
                        'the_percnt', 'the_verson', 'the_HCT_Cn', 'the_beam_X', 'the_beam_Y', 'the_beam_Z', 'the_Tumr_X', 'the_Tumr_Y', 'the_Tumr_Z', ...
                        's','dl','no_of_photons','M_raw')

                    clearvars p_in i_in m_in w_in
                    clearvars p_ot i_ot m_ot w_ot
                    clearvars s no_of_photons M_raw
                    clearvars fid the_percnt the_verson the_HCT_Cn the_beam_X the_Tumr_Z TheImage
                    clearvars mua_tumr mus_tumr g_tumr n_tumr
                end
            end
        end
    end
end
end

function [p_in, i_in, m_in, w_in, p_ot, i_ot, m_ot, w_ot, s, dl, no_of_photons, M_raw] = do_simulation (opt_glnd, opt_adps, opt_tumr, M_in, dl, wvlngth, nPhotonsReq, nExamplePaths,the_beam_X,the_beam_Y,the_beam_Z)
%% Geometry definition
model = MCmatlab.model;
model.G.nx                = size(M_in,1);  % Number of bins in the x direction
model.G.ny                = size(M_in,2);  % Number of bins in the y direction
model.G.nz                = size(M_in,3);  % Number of bins in the z direction
model.G.Lx                = size(M_in,1)*dl; % [cm] x size of simulation cuboid
model.G.Ly                = size(M_in,2)*dl; % [cm] y size of simulation cuboid
model.G.Lz                = size(M_in,3)*dl; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc;  % Media properties defined as a function at the end of this file
model.G.mediaPropParams     = {opt_adps(1),opt_adps(2),opt_adps(3),opt_adps(4),opt_glnd(1),opt_glnd(2),opt_glnd(3),opt_glnd(4),opt_tumr(1),opt_tumr(2),opt_tumr(3),opt_tumr(4)}; % A cell array that you can use to contain all sorts of inputs you would like to use inside the mediaPropertiesFunc
model.G.geomFunc            = @geometryDefinition;   % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.
model.G.geomFuncParams      = {M_in}; % A cell array that you can use to contain all sorts of inputs you would like to use inside the geomFunc
model = plot(model,'G');
clearvars opt_glnd opt_adps opt_tumr M_in

%% Monte Carlo simulation
model.MC.simulationTimeRequested  = 0.25;    % [min] Time duration of the simulation
model.MC.matchedInterfaces        = false;   % Assumes all refractive indices are the same
model.MC.boundaryType             = 1;       % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.nPhotonsRequested = nPhotonsReq;    % The time to run the MC simualtion for, in minutes. The number of photos launched will vary from run to run.
model.MC.nExamplePaths = nExamplePaths;      % The code stores the paths of the first N photons for subsequent visualization during the plotting.
model.MC.wavelength = wvlngth;               % [nm] Excitation wavelength, used for determination of optical properties for excitation light
clearvars wvlngth nPhotonsReq

%% For a pencil beam, the "focus" is just a point that the beam goes through, here set to be the center of the cuboid:
model.MC.lightSource.sourceType   = 0;    % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.silentMode = false;
beam_X = the_beam_X; beam_Y = the_beam_Y; beam_Z = the_beam_Z; beam_tht = 0; beam_phi = 0;
[ix_beam,iy_beam,iz_beam] = findIndex(beam_X,beam_Y,beam_Z,model.G.x,model.G.y,model.G.z);
while model.G.M_raw(ix_beam,iy_beam,iz_beam)==1
    beam_Z = iz_beam*dl-dl;
    iz_beam = iz_beam+1;
end
disp([ix_beam,iy_beam,iz_beam-1])
model.MC.lightSource.xFocus       = beam_X;   % [cm] x position of focus
model.MC.lightSource.yFocus       = beam_Y;   % [cm] y position of focus
model.MC.lightSource.zFocus       = beam_Z;   % [cm] z position of focus
model.MC.lightSource.theta        = beam_tht; % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = beam_phi; % [rad] Azimuthal angle of beam center axis
clearvars the_beam_X the_beam_Y the_beam_Z beam_X beam_Y beam_Z beam_phi beam_tht ix_beam iy_beam iz_beam

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
    t_scatter_stat = nan(length(fnsh_pnts),3+3+2+3+3+2+1); clearvars spratrs
    for idx = 1:length(strt_pnts)
        Photon_Path = t_model.MC.examplePaths(:,strt_pnts(idx)+1:fnsh_pnts(idx)-1);

        [ix_e,iy_e,iz_e] = findIndex(Photon_Path(1,end),  Photon_Path(2,end),  Photon_Path(3,end),  t_model.G.x,t_model.G.y,t_model.G.z); M_e = double(t_model.G.M_raw(ix_e,iy_e,iz_e)); clearvars ix_e iy_e iz_e
        [ix_1,iy_1,iz_1] = findIndex(Photon_Path(1,end-1),Photon_Path(2,end-1),Photon_Path(3,end-1),t_model.G.x,t_model.G.y,t_model.G.z); M_1 = double(t_model.G.M_raw(ix_1,iy_1,iz_1)); clearvars ix_1 iy_1 iz_1
        while (M_e==1&&M_1==1&&size(Photon_Path,2)>2)
            clearvars M_e M_1
            Photon_Path = Photon_Path(:,1:end-1);
            [ix_e,iy_e,iz_e] = findIndex(Photon_Path(1,end),  Photon_Path(2,end),  Photon_Path(3,end),  t_model.G.x,t_model.G.y,t_model.G.z); M_e = double(t_model.G.M_raw(ix_e,iy_e,iz_e)); clearvars ix_e iy_e iz_e
            [ix_1,iy_1,iz_1] = findIndex(Photon_Path(1,end-1),Photon_Path(2,end-1),Photon_Path(3,end-1),t_model.G.x,t_model.G.y,t_model.G.z); M_1 = double(t_model.G.M_raw(ix_1,iy_1,iz_1)); clearvars ix_1 iy_1 iz_1
        end
        clearvars M_e M_1

        [ix_f,iy_f,iz_f] = findIndex(Photon_Path(1,1),  Photon_Path(2,1),  Photon_Path(3,1),  t_model.G.x,t_model.G.y,t_model.G.z); M_f = double(t_model.G.M_raw(ix_f,iy_f,iz_f));
        [ix_e,iy_e,iz_e] = findIndex(Photon_Path(1,end),Photon_Path(2,end),Photon_Path(3,end),t_model.G.x,t_model.G.y,t_model.G.z); M_e = double(t_model.G.M_raw(ix_e,iy_e,iz_e));
        if size(Photon_Path,2)>2
            t_scatter_stat(idx,:) = [...
                Photon_Path(1,1), ...
                Photon_Path(2,1), ...
                Photon_Path(3,1), ...
                ix_f, ...
                iy_f, ...
                iz_f, ...
                (M_f), ...
                Photon_Path(4,1), ...
                Photon_Path(1,end), ...
                Photon_Path(2,end), ...
                Photon_Path(3,end), ...
                ix_e, ...
                iy_e, ...
                iz_e, ...
                (M_e), ...
                Photon_Path(4,end), ...
                sum(sqrt(sum((Photon_Path(1:3,2:end-0)-Photon_Path(1:3,1:end-1)).^2,1))), ...
                ];
        end
        clearvars Photon_Path M_f M_e ix_f iy_f iz_f ix_e iy_e iz_e
    end
    clearvars strt_pnts fnsh_pnts idx t_model

    idx_photons = ~isnan(t_scatter_stat(:,1));
    t_scatter_stat = t_scatter_stat(idx_photons,:); clearvars idx_photons

    no_of_photons = no_of_photons + size(t_scatter_stat,1);
    scatter_stat = [scatter_stat ; t_scatter_stat];
    disp('****************************************************')
    disp(['progress ... ',num2str(size(scatter_stat,1)/nExamplePaths*100),'%'])
    disp('****************************************************')
    clearvars surf_photons t_scatter_stat
end

clearvars model
% histcount
p_in = scatter_stat(:,1:3);                                                % source positions
i_in = scatter_stat(:,4:6);                                                % source indices
m_in = scatter_stat(:,7);                                                  % source media
w_in = scatter_stat(:,8);                                                  % source weights
p_ot = scatter_stat(:, 9:11);                                              % detedtor positions
i_ot = scatter_stat(:,12:14);                                              % detedtor indices
m_ot = scatter_stat(:,15);                                                 % detedtor media
w_ot = scatter_stat(:,16);                                                 % detedtor weights
s    = scatter_stat(:,17);                                                 % true photon pathlength
end
function M = geometryDefinition (X,Y,Z,parameters)
% Geometry function(s) (see readme for details)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% mediaPropertiesFunc) at each voxel location.
M = parameters{1} + 1;
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
mediaProperties(j).name  = 'adipose tissue';
mediaProperties(j).mua   = double(var{1}); % Absorption coefficient [cm^-1]
mediaProperties(j).mus   = double(var{2}); % Scattering coefficient [cm^-1]
mediaProperties(j).g     = double(var{3}); % Henyey-Greenstein scattering anisotropy
mediaProperties(j).n     = double(var{4}); % The refractive index

j = 3;
mediaProperties(j).name  = 'fibroglandular tissue';
mediaProperties(j).mua   = double(var{5}); % Absorption coefficient [cm^-1]
mediaProperties(j).mus   = double(var{6}); % Scattering coefficient [cm^-1]
mediaProperties(j).g     = double(var{7}); % Henyey-Greenstein scattering anisotropy
mediaProperties(j).n     = double(var{8}); % The refractive index

j = 4;
mediaProperties(j).name  = 'tumor/blood tissue';
mediaProperties(j).mua   = double(var{9}); % Absorption coefficient [cm^-1]
mediaProperties(j).mus   = double(var{10}); % Scattering coefficient [cm^-1]
mediaProperties(j).g     = double(var{11}); % Henyey-Greenstein scattering anisotropy
mediaProperties(j).n     = double(var{12}); % The refractive index
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
function [ix,iy,iz] = findIndex(x,y,z,X,Y,Z)
[~,ix] = min(abs(X-x));
[~,iy] = min(abs(Y-y));
[~,iz] = min(abs(Z-z));
end
function [mua,mus,g,musp] = calc_muas_based_HCT(HCT,lambda)
% Martina Meinke, Gerhard Müller, Jürgen Helfmann, and Moritz Friebel, "Empirical model functions to calculate hematocrit-dependent optical properties of human blood," Appl. Opt. 46, 1742-1753 (2007)
[mua_st,mus_st,g_st,musp_st] = Data_of_Standard_Optical_Parameters (lambda);
mua = nan(size(lambda));
mus = nan(size(lambda));
g   = nan(size(lambda));
musp= nan(size(lambda));
for i_lambda = 1:length(lambda)
    % mu_a
    if     ((250<=lambda(i_lambda)&&lambda(i_lambda)<=400)||(430<=lambda(i_lambda)&&lambda(i_lambda)<=600))&&(0.84<=HCT&&HCT<=42.1)
        mua(i_lambda) = 0.1233.*mua_st(i_lambda).*HCT;
    elseif  (600<=lambda(i_lambda)&&lambda(i_lambda)<=1100)&&(0.84<=HCT&&HCT<=42.1)
        mua(i_lambda) = 0.1206.*mua_st(i_lambda).*HCT;
    elseif  (400< lambda(i_lambda)&&lambda(i_lambda)< 430 )&&(0.84<=HCT&&HCT<=42.1)
        mua(i_lambda) = mean([0.1233 0.1206]).*mua_st(i_lambda).*HCT;
    else
        mua(i_lambda) = nan;
    end
    % mu_s
    if      (250<=lambda(i_lambda)&&lambda(i_lambda)<=1100)&&(0.84<=HCT&&HCT<=17.1)
        mus(i_lambda) = (-0.0015.*HCT.^2+0.1268.*HCT).*mus_st(i_lambda);
    elseif  (600<=lambda(i_lambda)&&lambda(i_lambda)<=1100)&&(17.1<=HCT&&HCT<=42.1)
        mus(i_lambda) = (-0.0015.*HCT.^2+0.1268.*HCT).*mus_st(i_lambda);
    elseif  (250<=lambda(i_lambda)&&lambda(i_lambda)< 600 )&&(17.1<=HCT&&HCT<=42.1)
        mus(i_lambda) = (-0.0015.*HCT.^2+0.1268.*HCT).*mus_st(i_lambda);
    else
        mus(i_lambda) = nan;
    end
    % mu_sp
    if      (250<=lambda(i_lambda)&&lambda(i_lambda)<=1100)&&(0.84<=HCT&&HCT<=17.1)
        musp(i_lambda) = 0.1167.*HCT.*musp_st(i_lambda);
    elseif  (600<=lambda(i_lambda)&&lambda(i_lambda)<=1100)&&(17.1<=HCT&&HCT<=42.1)
        musp(i_lambda) = 0.1167.*HCT.*musp_st(i_lambda);
    elseif  (250<=lambda(i_lambda)&&lambda(i_lambda)< 600 )&&(17.1<=HCT&&HCT<=42.1)
        musp(i_lambda) = 0.1167.*HCT.*musp_st(i_lambda);
    else
        musp(i_lambda) = nan;
    end
    % g
    g = (1-musp./mus);
    % if      (600<=lambda(i_lambda)&&lambda(i_lambda)<=1100)&&(0.84<=HCT&&HCT<=42.1)
    %     g(i_lambda) = (((-2.684e-6).*HCT.^2)+((-2.373e-4).*HCT)+1.003).*g_st(i_lambda);
    % elseif  (250<=lambda(i_lambda)&&lambda(i_lambda)< 600 )&&(0.84<=HCT&&HCT<=42.1)
    %     g(i_lambda) = (((-2.684e-6).*HCT.^2)+((-2.373e-4).*HCT)+1.003).*g_st(i_lambda);
    % else
    %     g(i_lambda) = nan;
    % end
end
end
function [mua_st,mus_st,g_st,musp_st] = Data_of_Standard_Optical_Parameters (lmbda)
% Martina Meinke, Gerhard Müller, Jürgen Helfmann, and Moritz Friebel, "Empirical model functions to calculate hematocrit-dependent optical properties of human blood," Appl. Opt. 46, 1742-1753 (2007)
db = [
250	7.52	27.9	0.877	3.43
255	7.93	27.2	0.88	3.27
260	8.6	26.6	0.883	3.1
265	9.41	25.7	0.883	3.01
270	9.82	25.4	0.882	2.99
275	9.93	25.2	0.884	2.93
280	9.74	25.3	0.889	2.8
285	9.19	25.5	0.899	2.59
290	8.24	25.9	0.912	2.28
295	7.04	26.7	0.9268	1.96
300	6	27.4	0.9389	1.68
305	5.52	27.7	0.9449	1.53
310	5.56	27.5	0.9442	1.54
315	5.95	26.9	0.9398	1.62
320	6.5	26	0.9335	1.73
325	7.13	25.4	0.928	1.83
330	7.69	25.2	0.9239	1.92
335	8.1	24.9	0.9208	1.97
340	8.38	24.7	0.9192	1.99
345	8.48	24.6	0.9191	1.99
350	8.37	24.6	0.9215	1.93
355	8.07	24.7	0.9257	1.84
360	7.69	24.7	0.9308	1.71
365	7.4	24.7	0.9345	1.61
370	7.37	24.5	0.9356	1.58
375	7.72	24	0.9323	1.62
380	8.58	23.2	0.9242	1.76
385	9.98	22.2	0.9105	1.99
390	11.86	21.1	0.892	2.29
395	14.21	20.1	0.869	2.64
400	16.86	19.1	0.845	2.96
405	19.66	18.3	0.824	3.23
410	21.8	17.9	0.81	3.39
415	22.66	18.6	0.812	3.48
420	21.54	19.8	0.831	3.35
425	18.62	21.8	0.861	3.03
430	14.93	23.8	0.893	2.55
435	11.52	25.3	0.92	2.02
440	8.95	26.8	0.9405	1.6
445	7.11	28	0.9543	1.28
450	5.79	29	0.9629	1.08
455	4.83	30	0.9685	0.946
460	4.11	30.5	0.9724	0.841
465	3.56	31	0.9756	0.757
470	3.13	31.3	0.9777	0.699
475	2.79	31.5	0.9793	0.652
480	2.52	31.8	0.9806	0.619
490	2.2	32	0.982	0.575
500	2.04	32.2	0.9831	0.545
510	1.98	32.2	0.9835	0.531
520	2.51	31.5	0.9838	0.51
530	3.97	29.5	0.9794	0.609
540	5.05	28.2	0.9755	0.691
550	4.4	29	0.9779	0.642
560	3.6	30.1	0.9804	0.59
570	4.56	28.8	0.9777	0.641
580	4.56	28.9	0.9771	0.662
590	1.9	32.1	0.9827	0.556
600	0.478	33.9	0.9854	0.496
610	0.17	34.2	0.9858	0.487
620	0.0812	34.3	0.9861	0.477
630	0.0496	34.2	0.9863	0.469
640	0.0348	34.3	0.9865	0.462
660	0.0251	34.1	0.9868	0.45
670	0.0239	33.9	0.9871	0.439
690	0.0243	33.6	0.9872	0.43
700	0.0246	33.4	0.9872	0.427
720	0.0284	32.9	0.9871	0.426
740	0.036	32.4	0.987	0.422
760	0.0461	31.8	0.9868	0.42
780	0.0558	31.2	0.9867	0.415
800	0.0641	30.8	0.9868	0.407
820	0.0762	30.5	0.9867	0.407
840	0.0853	30.5	0.9867	0.405
860	0.0953	29.7	0.9864	0.403
880	0.104	29.8	0.9861	0.413
900	0.106	28.4	0.9857	0.405
920	0.111	28.1	0.9854	0.411
940	0.117	27.7	0.9847	0.423
960	0.125	26.8	0.984	0.429
980	0.133	26.1	0.9836	0.428
1000	0.128	25.8	0.9837	0.421
1020	0.118	25.7	0.9839	0.412
1040	0.104	25.3	0.9837	0.412
1060	0.0879	25	0.9838	0.406
1080	0.0795	24.6	0.9841	0.392
1100	0.074	24.5	0.9842	0.387
];
mua_st  = interp1(db(:,1),db(:,2),lmbda,'linear','extrap'); mua_st  =  mua_st*10;
mus_st  = interp1(db(:,1),db(:,3),lmbda,'linear','extrap'); mus_st  =  mus_st*10;
g_st    = interp1(db(:,1),db(:,4),lmbda,'linear','extrap'); g_st    =       g_st;
musp_st = interp1(db(:,1),db(:,5),lmbda,'linear','extrap'); musp_st = musp_st*10;
end
function [n_water] = refractive_index_water(lambda_nm)
% Computes the refractive index of water as a function of wavelength in nm
% George M. Hale and Marvin R. Querry, "Optical Constants of Water in the 200-nm to 200-μm Wavelength Region," Appl. Opt. 12, 555-563 (1973)
% from https://refractiveindex.info/?book=H2O&page=Hale&shelf=main&utm_source=chatgpt.com
db = [...
    0.200	1.396
    0.225	1.373
    0.250	1.362
    0.275	1.354
    0.300	1.349
    0.325	1.346
    0.350	1.343
    0.375	1.341
    0.400	1.339
    0.425	1.338
    0.450	1.337
    0.475	1.336
    0.500	1.335
    0.525	1.334
    0.550	1.333
    0.575	1.333
    0.600	1.332
    0.625	1.332
    0.650	1.331
    0.675	1.331
    0.700	1.331
    0.725	1.330
    0.750	1.330
    0.775	1.330
    0.800	1.329
    0.825	1.329
    0.850	1.329
    0.875	1.328
    0.900	1.328
    0.925	1.328
    0.950	1.327
    0.975	1.327
    1.0	1.327
    1.2	1.324
    1.4	1.321
    1.6	1.317
    1.8	1.312
    2.0	1.306
    2.2	1.296
    2.4	1.279
    2.6	1.242
    2.65	1.219
    2.70	1.188
    2.75	1.157
    2.80	1.142
    2.85	1.149
    2.90	1.201
    2.95	1.292
    3.00	1.371
    3.05	1.426
    3.10	1.467
    3.15	1.483
    3.20	1.478
    3.25	1.467
    3.30	1.450
    3.35	1.432
    3.40	1.420
    3.45	1.410
    3.50	1.400
    3.6	1.385
    3.7	1.374
    3.8	1.364
    3.9	1.357
    4.0	1.351
    4.1	1.346
    4.2	1.342
    4.3	1.338
    4.4	1.334
    4.5	1.332
    4.6	1.330
    4.7	1.330
    4.8	1.330
    4.9	1.328
    5.0	1.325
    5.1	1.322
    5.2	1.317
    5.3	1.312
    5.4	1.305
    5.5	1.298
    5.6	1.289
    5.7	1.277
    5.8	1.262
    5.9	1.248
    6.0	1.265
    6.1	1.319
    6.2	1.363
    6.3	1.357
    6.4	1.347
    6.5	1.339
    6.6	1.334
    6.7	1.329
    6.8	1.324
    6.9	1.321
    7.0	1.317
    7.1	1.314
    7.2	1.312
    7.3	1.309
    7.4	1.307
    7.5	1.304
    7.6	1.302
    7.7	1.299
    7.8	1.297
    7.9	1.294
    8.0	1.291
    8.2	1.286
    8.4	1.281
    8.6	1.275
    8.8	1.269
    9.0	1.262
    9.2	1.255
    9.4	1.247
    9.6	1.239
    9.8	1.229
    10.0	1.218
    10.5	1.185
    11.0	1.153
    11.5	1.126
    12.0	1.111
    12.5	1.123
    13.0	1.146
    13.5	1.177
    14.0	1.210
    14.5	1.241
    15.0	1.270
    15.5	1.297
    16.0	1.325
    16.5	1.351
    17.0	1.376
    17.5	1.401
    18.0	1.423
    18.5	1.443
    19.0	1.461
    19.5	1.476
    20.0	1.480
    21.0	1.487
    22	1.500
    23	1.511
    24	1.521
    25	1.531
    26	1.539
    27	1.545
    28	1.549
    29	1.551
    30	1.551
    32	1.546
    34	1.536
    36	1.527
    38	1.522
    40	1.519
    42	1.522
    44	1.530
    46	1.541
    48	1.555
    50	1.587
    60	1.703
    70	1.821
    80	1.886
    90	1.924
    100	1.957
    110	1.966
    120	2.004
    130	2.036
    140	2.056
    150	2.069
    160	2.081
    170	2.094
    180	2.107
    190	2.119
    200	2.130
    ];
n_water = interp1(db(:,1).*1000,db(:,2),lambda_nm,'linear','extrap');
end
function [beta_st] = Specific_Refractive_Increment_beta(lmbda)
% Moritz Friebel and Martina Meinke, "Model function to calculate the refractive index of native hemoglobin in the wavelength range of 250-1100 nm dependent on concentration," Appl. Opt. 45, 2838-2842 (2006)
db = [...
250	0.00221
255	0.002155
260	0.002105
265	0.002069
270	0.002048
275	0.002042
280	0.002044
285	0.002047
290	0.002047
295	0.002037
300	0.00202
305	0.001999
310	0.001998
320	0.002007
330	0.002021
340	0.00201
350	0.001989
355	0.001985
360	0.001983
365	0.001912
370	0.00186
375	0.001816
380	0.001774
385	0.001732
390	0.001694
395	0.001668
400	0.001664
405	0.001701
410	0.001799
415	0.001985
420	0.002117
425	0.002195
430	0.002273
435	0.002227
440	0.00221
445	0.002184
450	0.002156
455	0.002131
460	0.002109
465	0.002092
470	0.002078
475	0.002067
480	0.002056
485	0.002045
490	0.002033
495	0.002019
500	0.002005
510	0.002009
520	0.001983
530	0.001966
540	0.001981
550	0.001998
560	0.001992
570	0.001988
580	0.002004
590	0.002015
600	0.001988
610	0.001967
620	0.001964
630	0.00196
640	0.001954
660	0.001958
680	0.00197
700	0.001992
720	0.001979
740	0.001955
760	0.001958
780	0.00196
800	0.001939
820	0.00192
840	0.001935
860	0.001951
880	0.001982
900	0.001998
920	0.002011
940	0.002015
960	0.002021
980	0.002017
1000	0.002052
1020	0.002049
1040	0.002044
1060	0.00204
1080	0.002044
1100	0.002056
];
beta_st = interp1(db(:,1),db(:,2),lmbda,'linear','extrap');
end
