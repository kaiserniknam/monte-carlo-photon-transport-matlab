function [] = Photon_45_4 ()
% Repository group: 05_digital_breast_phantoms
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: simulating phantoms prepared by Diego
% Similar to Photon_44: but in better intervals
% analyzing simulated data: general saving but include p-values

clc
close all

% Optical properties of breast and tumor [1-4]
% mua_glnd = 0.55; mus_glnd = 332.7; g_glnd = 0.965; n_glnd = 1.4;
% mua_adps = 0.82; mus_adps = 313.8; g_adps = 0.976; n_adps = 1.4;
set_of_percnt = [25,50,75]; % percent of adipose
set_of_versns = [0,1,2]; % version of each density
set_of_beam_X = [-1.0,0.0,+1.0]; % beam_X location (in cm)
set_of_Tumr_Z = [1.0,1.25,1.50,1.75]; % Tumor_Z location (in cm)
set_of_HCT_Cn = 0:2.5:40; % HCT (in %)
N_bins = 201; Lx = 20.4; Ly = 6.48; Lz = 5.16;
d_edges = linspace(0,sqrt((Lx/2).^2+(Ly/2).^2+(Lz/2).^2),N_bins+1);
set_of_I_s = nan(length(set_of_percnt),length(set_of_versns),length(set_of_beam_X),length(set_of_Tumr_Z),length(set_of_HCT_Cn),5,N_bins);
counter = 0;

for i_percent = 1:length(set_of_percnt)
    for i_versns = 1:length(set_of_versns)
        for i_beam_X = 1:length(set_of_beam_X)
            for i_Tumr_Z = 1:length(set_of_Tumr_Z)
                for i_HCT_Cn = 1:length(set_of_HCT_Cn)
                    counter = counter + 1; disp(num2str(counter))
                    the_percnt = set_of_percnt(i_percent);
                    the_verson = set_of_versns(i_versns);
                    the_beam_X = set_of_beam_X(i_beam_X);
                    the_Tumr_Z = set_of_Tumr_Z(i_Tumr_Z);
                    the_HCT_Cn = set_of_HCT_Cn(i_HCT_Cn);
                    data = load(['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_45_Pd',num2str(the_percnt),'_',num2str(the_verson),'_',num2str(the_HCT_Cn),'_',num2str(the_beam_X),'_',num2str(the_Tumr_Z),'_y_-z_x.mat']);
                    refc = load(['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_45_Pd',num2str(the_percnt),'_',num2str(the_verson),'_',num2str(0),         '_',num2str(the_beam_X),'_',num2str(the_Tumr_Z),'_y_-z_x.mat']);
                    fun_x = @mean; fun_y = @sum; TheOutFun = @(x)(-log(x));
                    x_temp = sqrt(...
                        (data.p_in(:,1)-data.p_ot(:,1)).^2+ ...
                        (data.p_in(:,2)-data.p_ot(:,2)).^2+ ...
                        (data.p_in(:,3)-data.p_ot(:,3)).^2);
                    y_temp = data.w_ot/data.no_of_photons;
                    [~,~,index_in] = histcounts(x_temp,d_edges);
                    % x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
                    % y_bind = accumarray(index_in,y_temp,[],fun_y,nan);

                    x_refc = sqrt(...
                        (refc.p_in(:,1)-refc.p_ot(:,1)).^2+ ...
                        (refc.p_in(:,2)-refc.p_ot(:,2)).^2+ ...
                        (refc.p_in(:,3)-refc.p_ot(:,3)).^2);
                    y_refc = refc.w_ot/refc.no_of_photons;
                    [~,~,index_rf] = histcounts(x_refc,d_edges);
                    U = [index_rf;index_in]; [unique_U,~,~] = unique(U); clearvars U
                    p_bind = nan(length(unique_U),5);
                    for i_U = 1:size(unique_U,1)
                        idx_from_in = index_in==unique_U(i_U); p_bind(i_U,1) = fun_x(x_temp(idx_from_in)); p_bind(i_U,2) = fun_y(y_temp(idx_from_in));
                        idx_from_rf = index_rf==unique_U(i_U); p_bind(i_U,3) = fun_x(x_refc(idx_from_rf)); p_bind(i_U,4) = fun_y(y_refc(idx_from_rf));
                        if sum(idx_from_in)>0&&sum(idx_from_rf)
                            p_bind(i_U,5) = ranksum(y_temp(idx_from_in),y_refc(idx_from_rf));
                        end
                        clearvars idx_from_in idx_from_rf
                    end
                    set_of_I_s(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,1:5,1:size(p_bind,1)) = p_bind.';

                    clearvars x_refc y_refc index_rf unique_U i_U refc p_bind
                    clearvars the_percnt the_verson the_beam_X the_Tumr_Z the_HCT_Cn data
                    clearvars index_in x_temp y_temp x_bind y_bind fun_x fun_y TheOutFun
                end
            end
        end
    end
end
clearvars i_beam_X i_HCT_Cn i_percent i_Tumr_Z i_versns
save(['Photon_45_4_',num2str(N_bins),'.mat'],'set_of_percnt','set_of_versns','set_of_beam_X','set_of_Tumr_Z','set_of_HCT_Cn','N_bins','d_edges','set_of_I_s')
end
