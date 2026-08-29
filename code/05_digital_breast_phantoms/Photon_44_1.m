function [] = Photon_44_1 ()
% Repository group: 05_digital_breast_phantoms
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: simulating phantoms prepared by Diego
% Similar to Photon_42: comparisons between Tumor and No Tumor cases, with the tumor composed of blood. The source is placed at four different locations. Tumor radius is 1.5/2 cm. Simulations are run for different hematocrit (HCT) levels.
% analyzing simulated data: only validation of data

clc
close all

% Optical properties of breast and tumor [1-4]
% mua_glnd = 0.55; mus_glnd = 332.7; g_glnd = 0.965; n_glnd = 1.4;
% mua_adps = 0.82; mus_adps = 313.8; g_adps = 0.976; n_adps = 1.4;
set_of_percnt = [25,50,75]; % percent of adipose
set_of_versns = [0,1,2]; % version of each density
set_of_beam_X = [-6,-3,0,+3]; % beam_X location (in cm)
set_of_Tumr_Z = [0.75,1.25,1.75]; % Tumor_Z location (in cm)
set_of_HCT_Cn = [0,10,20,30,40]; % HCT (in %)
set_of_endmat = nan(length(set_of_percnt),length(set_of_versns),length(set_of_beam_X),length(set_of_Tumr_Z),length(set_of_HCT_Cn),3);
set_of_endmap = nan(length(set_of_percnt),length(set_of_versns),length(set_of_beam_X),length(set_of_Tumr_Z),length(set_of_HCT_Cn),510,162);
counter = 0;

for i_percent = 1:length(set_of_percnt)
    for i_versns = 1:length(set_of_versns)
        for i_beam_X = 1:length(set_of_beam_X)
            for i_Tumr_Z = 1:length(set_of_Tumr_Z)
                for i_HCT_Cn = 1:length(set_of_HCT_Cn)

                    counter = counter + 1;
                    the_percnt = set_of_percnt(i_percent);
                    the_verson = set_of_versns(i_versns);
                    the_beam_X = set_of_beam_X(i_beam_X);
                    the_Tumr_Z = set_of_Tumr_Z(i_Tumr_Z);
                    the_HCT_Cn = set_of_HCT_Cn(i_HCT_Cn);
                    data = load(['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_44_Pd',num2str(the_percnt),'_',num2str(the_verson),'_',num2str(the_HCT_Cn),'_',num2str(the_beam_X),'_',num2str(the_Tumr_Z),'_y_-z_x.mat']);

                    % Display information on the percentage of adipose tissue (1), and the percentage of photons that ended in each tissue type: air, adipose, and fibrogranular.
                    set_of_endmat(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,1) = sum(data.m_ot==1)/length(data.m_ot);
                    set_of_endmat(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,2) = sum(data.m_ot==2)/length(data.m_ot);
                    set_of_endmat(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_HCT_Cn,3) = sum(data.m_ot==3)/length(data.m_ot);
                    disp([num2str(counter),' - Pd = ',num2str(the_percnt),' ~ ',num2str(100*sum(reshape(data.M_raw==2,[],1))/(sum(reshape(data.M_raw==2,[],1))+sum(reshape(data.M_raw==3,[],1))),'%.2f'),', version = ',num2str(the_verson),', beam X = ',num2str(the_beam_X),', tumor Z = ',num2str(the_Tumr_Z),', HCT = ',num2str(the_HCT_Cn)])

                    clearvars data the_percnt the_verson the_HCT_Cn the_beam_X the_Tumr_Z
                end
            end
        end
    end
end
clearvars i_tissue i_percent i_versns i_beam_X i_Tumr_Z i_HCT_Cn
end
