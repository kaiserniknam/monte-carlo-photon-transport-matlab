function [] = Photon_49_1 ()
% Repository group: 06_laboratory_phantom_models
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: Replication of Arnab's Experiment: Dissolving x grams of TiO₂ in a solution of 10 g Sodium Polyacrylate (PAS) in 440 mL of water.
% the same as Photon_40, but between 0 & 3
% general analysis

clc
close all

z_air = 0.0; % the thickness of air layer
nPhotonsReq = 1e5; % number of requested photon
nExamplePaths = 1000000; % number of example photon paths
% Optical & size properties
Lx = 29.1; Ly = 29.1; Lz = 10.0; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
beam_phi = 0;     % polar angle of beam in cm
beam_tht = 0;     % azimuthal angle of beam in cm
beam_X = 0.0;     % X deviation of beam in cm
beam_Y = 0.0;     % Y deviation of beam in cm

N_bins = 151;
set_of_g = (0:0.1:3)+eps;
set_of_data     = nan(length(set_of_g),1+2+2+1);
set_of_cnst_DPF = nan(length(set_of_g),1);
set_of_idst_DPF = nan(length(set_of_g),2);
set_of_savg_DPF = nan(length(set_of_g),N_bins);
set_of_true_DPF = nan(length(set_of_g),N_bins);

for i_g = 1:length(set_of_g)
    % read dbase
    the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_49_g_',sprintf('%.2f',set_of_g(i_g)),'.mat'];
    t_db = load(the_filename); clearvars the_filename
    mua = t_db.mua;
    mus = t_db.mus;
    set_of_data(i_g,1) = set_of_g(i_g);
    set_of_data(i_g,2) = mua;
    set_of_data(i_g,3) = mus;
    % disp(['mua = ',num2str(mua,'%.2f'),', mus = ',num2str(mus,'%.0f')])
    disp([num2str(set_of_g(i_g)),',',num2str(mua,'%.2f'),',',num2str(mus,'%.0f'),',',num2str(t_db.g),',',num2str(t_db.n)])

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
    d_diff_edges = linspace(0,sqrt((Lx/2).^2+(Ly/2).^2),N_bins+1); [~,~,ind_diff] = histcounts(sqrt(t_db.x(t_db.c==0).^2+t_db.y(t_db.c==0).^2+t_db.z(t_db.c==0).^2),d_diff_edges); % d_diffuse bins
    clearvars d_trns_edges d_diff_edges



    % I vs. s (scatterplot)
    i_fig = 1; TheOutFun = @(x)(-log(x)); TheLineStyle = 'none'; TheMarker = 'o'; x_label = 's (cm)'; y_label = 'OD (a.u.)';
    TheLegend = ['TiO_2 = ',num2str(set_of_g(i_g)),' g'];
    % (diffuse)
    i_subplot = 1; index_in = ind_diff; TheCode = 0; TheTitle = 'I vs. s'; TheColor = [(i_g-1)/(length(set_of_g)-1) 0 1-(i_g-1)/(length(set_of_g)-1)];
    x_temp = t_db.s(t_db.c==TheCode); y_temp = t_db.w(t_db.c==TheCode);
    x_temp = x_temp; y_temp = TheOutFun(y_temp);
    if get_plot(mua,mus)
        figure(i_fig),
        plot(x_temp,y_temp,'Color',TheColor,'DisplayName',TheLegend,'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','w','MarkerFaceColor',TheColor,'MarkerSize',8), hold on
        xlim([0 150]), xlabel(x_label), ylabel(y_label), title(TheTitle), axis([0 50 0 5])
        set(gca,'fontsize',16), axis square, grid on, hold on, legend('show','Location','northeast')
    end
    clearvars i_subplot index_in x_temp y_temp TheCode TheTitle mdl
    clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label TheLegend TheColor



    % s vs. d
    fun_x = @mean; fun_y = @mean;
    index_in = ind_diff; TheCode = 0;
    x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2);
    y_temp = t_db.s(t_db.c==TheCode);
    x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
    y_bind = accumarray(index_in,y_temp,[],fun_y,nan);
    set_of_savg_DPF(i_g,1:length(y_bind)) = y_bind./x_bind;
    clearvars fun_x fun_y index_in TheCode x_temp y_temp x_bind y_bind



    % I vs. d
    i_fig = 2; fun_x = @mean; fun_y = @sum; TheOutFun = @(x)(-log(x)); fun_s = @(x)(sum(exp(-mua.*x))./t_db.no_of_photons); TheLineStyle = '-'; TheMarker = 'none'; x_label = 'd (cm)'; y_label = 'OD (a.u.)'; Therloess = 0.0; blnFit = false;
    TheLegend = ['TiO_2 = ',num2str(set_of_g(i_g)),' g'];
    % (diffuse)
    i_subplot = 1; index_in = ind_diff; TheCode = 0; TheTitle = 'I vs. d'; TheColor = [(i_g-1)/(length(set_of_g)-1) 0 1-(i_g-1)/(length(set_of_g)-1)];
    x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2);
    y_temp = t_db.w(t_db.c==TheCode)./t_db.no_of_photons;
    s_temp = t_db.s(t_db.c==TheCode);

    x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
    y_bind = accumarray(index_in,y_temp,[],fun_y,nan); y_bind = TheOutFun(y_bind);
    s_bind = accumarray(index_in,s_temp,[],fun_s,nan); s_bind = TheOutFun(s_bind)./mua./x_bind;

    mdl = fitlm(x_bind(~isnan(x_bind)),y_bind(~isnan(x_bind)));
    p = mdl.Coefficients.pValue(2);
    if p<=0.05
        set_of_cnst_DPF(i_g,1) = (x_bind(~isnan(x_bind))\y_bind(~isnan(x_bind)))/mua;
        set_of_idst_DPF(i_g,1) = mdl.Coefficients.Estimate(2)/mua;
        set_of_idst_DPF(i_g,2) = mdl.Coefficients.Estimate(1)/mua;
    else
        set_of_cnst_DPF(i_g,1) = nan;
        set_of_idst_DPF(i_g,1) = nan;
        set_of_idst_DPF(i_g,2) = nan;
    end
    set_of_true_DPF(i_g,1:length(s_bind)) = s_bind;

    [m,b,p] = make_a_subplot(i_fig,i_subplot,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,TheColor,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit,TheLegend);
    set_of_data(i_g,4) = m;
    set_of_data(i_g,5) = b;
    set_of_data(i_g,6) = p;
    clearvars m b p



    if get_plot(mua,mus)
        figure(i_fig), xlabel(x_label), ylabel(y_label), title(TheTitle), axis([0 20 0 20])
        set(gca,'fontsize',16), axis square, grid on, hold on, legend('show','Location','southeast')
    end
    clearvars i_subplot index_in x_temp y_temp s_bind s_temp TheCode TheTitle mdl
    clearvars i_fig fun_x fun_y fun_s TheOutFun TheLineStyle TheMarker x_label y_label Therloess blnFit



    % corrected DPF vs. d
    i_fig = 3; fun_x = @mean; fun_y = @(x)(sum(exp(-mua.*x))./t_db.no_of_photons); TheOutFun = @(x)(log(x)./mua); TheLineStyle = '-'; TheMarker = 'none'; x_label = 'd (cm)'; y_label = 'corrected DPF (a.u.)'; Therloess = 0.5; blnFit = false;
    % (diffuse)
    i_subplot = 1; index_in = ind_diff; TheCode = 0; TheTitle = 'diffuse';
    x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); y_temp = t_db.s(t_db.c==TheCode);
    x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
    y_bind = accumarray(index_in,y_temp,[],fun_y,nan); y_bind = -TheOutFun(y_bind)./x_bind;
    TheColor = [(i_g-1)/(length(set_of_g)-1) 0 1-(i_g-1)/(length(set_of_g)-1)];
    TheLegend = ['TiO_2 = ',num2str(set_of_g(i_g)),' g'];
    if get_plot(mua,mus)
        figure(i_fig), xlabel(x_label), ylabel(y_label), yscale("log"), title(TheTitle), % axis([0 10 0 150])
        plot(x_bind,y_bind,'Color',TheColor,'DisplayName',TheLegend,'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','w','MarkerFaceColor',TheColor,'MarkerSize',8,'LineWidth',1.5), hold on
        set(gca,'fontsize',16), axis square, grid on, hold on, legend('show','Location','northeast')
    end
    clearvars i_subplot index_in x_temp y_temp x_bind y_bind TheCode TheTitle TheColor TheLegend
    clearvars i_fig fun_x fun_y TheOutFun TheLineStyle TheMarker x_label y_label Therloess blnFit



    clearvars t_db the_filename mua mus ind_diff ind_trns
end
save('Photon_49_1.mat','set_of_data','set_of_g','set_of_cnst_DPF','set_of_idst_DPF','set_of_savg_DPF','set_of_true_DPF')
end

function [m,b,p] = make_a_subplot(i_fig,i_subplot,index_in,x_temp,y_temp,fun_x,fun_y,mua,mus,TheColor,TheOutFun,TheLineStyle,TheMarker,Therloess,blnFit,TheLegend)
x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
y_bind = accumarray(index_in,y_temp,[],fun_y,nan);
mdl = fitlm(x_bind,TheOutFun(y_bind));
p = mdl.Coefficients.pValue(2);
if p<=0.05
    b = mdl.Coefficients.Estimate(1);
    m = mdl.Coefficients.Estimate(2);
else
    b = nan;
    m = nan;
end
if get_plot(mua,mus)
    figure(i_fig)
    if Therloess<0 & length(x_bind)>3
        y_bind = smooth(x_bind,y_bind,Therloess,'rloess');
    end
    plot(x_bind,TheOutFun(y_bind),'Color',TheColor,'DisplayName',TheLegend,'Marker',TheMarker,'LineStyle',TheLineStyle,'MarkerEdgeColor','w','MarkerFaceColor',TheColor,'MarkerSize',8,'LineWidth',1.5), hold on
    if blnFit
        plot(x_bind,m*x_bind+b,'Color',TheColor,'HandleVisibility','off'), hold on
    end
end
end
function [out] = get_plot(mua,mus)
out = true;
end
