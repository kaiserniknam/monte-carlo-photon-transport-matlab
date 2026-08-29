function [] = Photon_33_6 ()
% Repository group: 07_dpf_parameter_studies
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: w & s vs. d over mu_s & mu_a: more photons, more combinations
% the analysis: population stat - to plot p and q for mu_a's and mu_s's
% Arnab dye experiment

clc
close all
set(0,'DefaultFigureWindowStyle','docked')

% optical & geometry properties
% z_air = 0.0; % the thickness of air layer
Lx = 29.1; Ly = 29.1; Lz = 05.0; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
% beam_phi = 0;     % polar angle of beam in cm
% beam_tht = 0;     % azimuthal angle of beam in cm
% beam_X = 0.0;     % X deviation of beam in cm
% beam_Y = 0.0;     % Y deviation of beam in cm

n = 1.4; g = 0.95;
set_of_mua = (.00:.05:0.5);
set_of_mus = (000:050:500);
N_bins = 200;
d_edges = linspace(0,sqrt((Lx/2).^2+(Ly/2).^2),N_bins+1); clearvars Lx Ly Lz
d_cntrs = 1/2*(d_edges(1:end-1)+d_edges(2:end-0));
d_detector = 0.8;
set_of_Ds = [2,3,4.5];
set_of_Ds = 1.5:0.5:4.5;
set_of_ODs = nan(length(set_of_mua),length(set_of_mus),length(set_of_Ds));

clearvars N_bins

for i_a = 1:length(set_of_mua)
    for i_s = 1:length(set_of_mus)
        % read dbase
        mua = set_of_mua(i_a);
        mus = set_of_mus(i_s);
        the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_33_mua_',sprintf('%.2f',mua),'_mus_',sprintf('%.2f',mus),'.mat'];
        t_db = load(the_filename); clearvars the_filename
        disp(['mua = ',num2str(mua,'%.2f'),', mus = ',num2str(mus,'%.0f')])

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
        disp(num2str([length(u) sum(u==u_unique(i_fx)) sum(u==u_unique(i_fx))/length(u)*100]))
        t_db.no_of_photons = t_db.no_of_photons-sum(u~=u_unique(i_fx));
        clearvars u u_unique freq i_f i_fx

        % 1-D sorting
        [~,~,ind_diff] = histcounts(sqrt(t_db.x(t_db.c==0).^2+t_db.y(t_db.c==0).^2+t_db.z(t_db.c==0).^2),d_edges); % d_diffuse bins
        clearvars d_trns_edges d_diff_edges

        % I vs. d
        fun_x = @mean; fun_y = @sum; TheOutFun = @(x)(-log(x));
        index_in = ind_diff; TheCode = 0;
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2);
        y_temp = t_db.w(t_db.c==TheCode)./t_db.no_of_photons;
        x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
        y_bind = accumarray(index_in,y_temp,[],fun_y,nan);
        for pnt = 1:length(set_of_Ds)
            idx_involved = abs(x_bind-set_of_Ds(pnt))<=d_detector/2;
            set_of_ODs(i_a,i_s,pnt) = TheOutFun(sum(y_bind(idx_involved))*d_detector/2/pi/set_of_Ds(pnt));
            clearvars idx_involved
        end
        clearvars pnt fun_x fun_y TheOutFun index_in TheCode x_temp y_temp x_bind y_bind
        disp(num2str(set_of_ODs(i_a,i_s,:)))
        disp('----------------------------')

        TheColor = get_color(1-mus/max(set_of_mus),mus/max(set_of_mus));
        subplot(3,4,i_a)
        plot(set_of_Ds,squeeze(set_of_ODs(i_a,i_s,:)),'Color',TheColor,'DisplayName',['\mu_s = ',num2str(mus,'%.0f')],'Marker','none','LineStyle','-','MarkerEdgeColor','w','MarkerFaceColor',TheColor,'MarkerSize',8,'LineWidth',1.5), hold on
        xlabel('source separation (cm)'), ylabel('OD (a.u.)'), title(['\mu_a = ',num2str(mua,'%.2f')])
        set(gca,'fontsize',12), axis tight, grid on, hold on, axis([1 5.5 4 25])

        clearvars mdl p blnFit blnShow fun_x fun_y fun_s i_fig index_in TheCode TheOutFun x_temp y_temp s_temp x_bind y_bind s_bind
        clearvars t_db the_filename mua mus ind_diff ind_trns TheColor
    end
end
clearvars i_a i_s
end

function [out] = get_color(mua_norm,mus_norm)
col_mua_0_mus_0 = [1 1 0];
col_mua_0_mus_1 = [1 0 0];
col_mua_1_mus_0 = [0 0 1];
col_mua_1_mus_1 = [1 0.0 1];
% out = [mua_norm 0 mus_norm];
out = ...
    (1-mua_norm)*(1-mus_norm)*col_mua_0_mus_0 + ...
    (  mua_norm)*(1-mus_norm)*col_mua_1_mus_0 + ...
    (1-mua_norm)*(  mus_norm)*col_mua_0_mus_1 + ...
    (  mua_norm)*(  mus_norm)*col_mua_1_mus_1;
end
