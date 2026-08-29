function [] = Photon_53_1 ()
% Repository group: 07_dpf_parameter_studies
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% This version computes weight (w) and pathlength (s) versus source-detector distance (d),
% across a range of μₐ and μₛ values. It is adapted from version 33 for Arnab,
% with simulations performed over multiple tissue thicknesses.

clc
close all

% optical & geometry properties
% z_air = 0.0; % the thickness of air layer
% Optical & size properties
Lx = 29.1; Ly = 29.1; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
% beam_phi = 0;     % polar angle of beam in cm
% beam_tht = 0;     % azimuthal angle of beam in cm
% beam_X = 0.0;     % X deviation of beam in cm
% beam_Y = 0.0;     % Y deviation of beam in cm

n = 1.4; g = 0.90;
set_of_mua = 0.4;
set_of_mus = 5.0:5.0:25.;
set_of_Lzs = 0.5:0.5:5.0;
d_detector = 0.8;

set_of_Ts = nan(length(set_of_mua),length(set_of_mus),length(set_of_Lzs));
for i_a = 1:length(set_of_mua)
    for i_s = 1:length(set_of_mus)
        for i_z = 1:length(set_of_Lzs)
            mua = set_of_mua(i_a);
            mus = set_of_mus(i_s);
            Lz  = set_of_Lzs(i_z);
            TheColor = [(i_s-1)/(length(set_of_mus)-1) 0 1-(i_s-1)/(length(set_of_mus)-1)];

            % read dbase
            the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_53_mua_',sprintf('%.2f',set_of_mua(i_a)),'_mus_',sprintf('%.2f',set_of_mus(i_s)),'_Lz_',num2str(set_of_Lzs(i_z)),'.mat' ];
            t_db = load(the_filename); clearvars the_filename
            disp(['mua = ',num2str(mua,'%.2f'),', mus = ',num2str(mus,'%.0f'),', Lz = ',num2str(Lz,'%.2f')])

            % removing noisy points
            u = round(-log(t_db.w)./t_db.s,4); u_unique = unique(u); freq = nan(size(u_unique));
            for i_f = 1:length(u_unique), freq(i_f) = sum(u==u_unique(i_f)); end
            [~,i_fx] = max(freq);
            t_db.p_in = t_db.p_in(u==u_unique(i_fx),:);
            t_db.p_ot = t_db.p_ot(u==u_unique(i_fx),:);
            t_db.s = t_db.s(u==u_unique(i_fx));
            t_db.w = t_db.w(u==u_unique(i_fx));
            t_db.d = t_db.d(u==u_unique(i_fx));
            t_db.c = t_db.c(u==u_unique(i_fx));
            disp(num2str([length(u) sum(u==u_unique(i_fx)) sum(u==u_unique(i_fx))/length(u)*100]))
            t_db.no_of_photons = t_db.no_of_photons-sum(u~=u_unique(i_fx));
            clearvars u u_unique freq i_f i_fx

            % (transmittance)
            n_d = sum(t_db.c==1&t_db.d<=d_detector/2);
            n_t = t_db.no_of_photons;
            set_of_Ts(i_a,i_s,i_z) = n_d/n_t;
            clearvars n_d n_t

            clearvars t_db the_filename mua mus ind_diff ind_trns
        end
        figure(1)
        plot(set_of_Lzs,squeeze(set_of_Ts(i_a,i_s,:)),'LineWidth',1,'Color',TheColor,'DisplayName',['\mu_a = ',num2str(set_of_mua(i_a)),' cm^{-1}, \mu_s=',num2str(set_of_mus(i_s)),' cm^{-1}: simulation'],'Marker','square','MarkerSize',12,'MarkerFaceColor',TheColor,'MarkerEdgeColor','k','LineStyle','-'), hold on
        % plot(set_of_Lzs,exp(-set_of_Lzs./sqrt(3.*set_of_mua(i_a).*(set_of_mua(i_a)+set_of_mus(i_s).*(1-g)))),'LineWidth',1,'Color',TheColor,'DisplayName',['\mu_a = ',num2str(set_of_mua(i_a)),' cm^{-1}, \mu_s=',num2str(set_of_mus(i_s)),' cm^{-1}: analytical'],'Marker','none','MarkerSize',12,'MarkerFaceColor',TheColor,'MarkerEdgeColor','k','LineStyle',':'), hold on
    end
end
clearvars i_a i_s i_z

figure(1)
xlabel('z (cm)'), set(gca,'xtick',[set_of_Lzs]), xlim([0 max(set_of_Lzs)])
ylabel('T (%)'),  set(gca,'ytick',0:0.1:1),      ylim([0 1])
title('transmittance for multiple pairs of \mu_a & \mu_s vs. depth (z)'),
set(gca,'fontsize',24), axis square,
legend('show','Location','northeast')
save('mus_5_25.mat','set_of_mua','set_of_mus','g','n','set_of_Ts','set_of_Lzs')
end
