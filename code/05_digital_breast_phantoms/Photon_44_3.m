function [] = Photon_44_3 ()
% Repository group: 05_digital_breast_phantoms
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: simulating phantoms prepared by Diego
% Similar to Photon_42: comparisons between Tumor and No Tumor cases, with the tumor composed of blood. The source is placed at four different locations. Tumor radius is 1.5/2 cm. Simulations are run for different hematocrit (HCT) levels.
% analyzing simulated data: I's for no-Tumor cases

clc
close all
db33 = load("DB/Photon_33_5.mat");

% Optical properties of breast and tumor [1-4]
% mua_glnd = 0.55; mus_glnd = 332.7; g_glnd = 0.965; n_glnd = 1.4;
% mua_adps = 0.82; mus_adps = 313.8; g_adps = 0.976; n_adps = 1.4;
set_of_percnt = [25,50,75]; % percent of adipose
set_of_versns = [0,1,2]; % version of each density
set_of_beam_X = [-6,-3,0,+3]; % beam_X location (in cm)
set_of_Tumr_Z = [0.75,1.25,1.75]; % Tumor_Z location (in cm)
set_of_HCT_Cn = [0,10,20,30,40]; % HCT (in %)
% N_bins = 150; Lx = 20.4; Ly = 6.48; Lz = 5.16;
% d_edges = linspace(0,sqrt((Lx/2).^2+(Ly/2).^2+(Lz/2).^2),N_bins+1);
% db33.d_cntrs = 1/2.*(d_edges(1:end-1)+d_edges(2:end-0));
N1 = 17; N2 = 37;
d_edges = db33.d_edges;
set_of_Is = nan(length(set_of_percnt),length(set_of_versns),length(set_of_beam_X),length(set_of_Tumr_Z),length(db33.d_cntrs));
set_of_Es = nan(length(set_of_percnt),length(set_of_versns),length(set_of_beam_X),length(set_of_Tumr_Z),length(db33.d_cntrs),3);
set_of_1s = nan(length(set_of_percnt),length(set_of_versns),length(set_of_beam_X),length(set_of_Tumr_Z),2);
set_of_2s = nan(length(set_of_percnt),length(set_of_versns),length(set_of_beam_X),length(set_of_Tumr_Z),2);
set_of_Ts = nan(length(set_of_percnt),length(set_of_versns),length(set_of_beam_X),length(set_of_Tumr_Z),2);

for i_percent = 1:length(set_of_percnt)
    for i_versns = 1:length(set_of_versns)
        for i_beam_X = 1:length(set_of_beam_X)
            for i_Tumr_Z = 1:length(set_of_Tumr_Z)
                for i_HCT_Cn = 1 %:length(set_of_HCT_Cn)
                    the_percnt = set_of_percnt(i_percent);
                    the_verson = set_of_versns(i_versns);
                    the_beam_X = set_of_beam_X(i_beam_X);
                    the_Tumr_Z = set_of_Tumr_Z(i_Tumr_Z);
                    the_HCT_Cn = set_of_HCT_Cn(i_HCT_Cn);
                    data = load(['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_44_Pd',num2str(the_percnt),'_',num2str(the_verson),'_',num2str(the_HCT_Cn),'_',num2str(the_beam_X),'_',num2str(the_Tumr_Z),'_y_-z_x.mat']);
                    disp(['Pd = ',num2str(the_percnt,'%2.0f'),', v = ',num2str(the_verson,'%1.0f'),', HCT = ',num2str(the_HCT_Cn,'%1.0f'),', beam X = ',num2str(the_beam_X,'%2.0f'),', Tumor Z = ',num2str(the_Tumr_Z,'%1.2f')])

                    fun_x = @mean; fun_y = @sum; TheOutFun = @(x)(-log(x));
                    x_temp = sqrt(...
                        (data.p_in(:,1)-data.p_ot(:,1)).^2+ ...
                        (data.p_in(:,2)-data.p_ot(:,2)).^2+ ...
                        (data.p_in(:,3)-data.p_ot(:,3)).^2); y_temp = data.w_ot/data.no_of_photons;
                    [~,~,index_in] = histcounts(x_temp,d_edges);
                    x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
                    y_bind = accumarray(index_in,y_temp,[],fun_y,nan);
                    set_of_Is(i_percent,i_versns,i_beam_X,i_Tumr_Z,1:length(x_bind)) = TheOutFun(y_bind);



                    % whole signal
                    A = nan(length(db33.set_of_mua),length(db33.set_of_mus));
                    for i_a = 1:length(db33.set_of_mua)
                        for i_s = 1:length(db33.set_of_mus)
                            A(i_a,i_s) = r_squared(squeeze(db33.set_of________I(i_a,i_s,:)),squeeze(set_of_Is(i_percent,i_versns,i_beam_X,i_Tumr_Z,:)));
                        end
                    end
                    [~,idx] = min(1-A(:));
                    [row,col] = ind2sub(size(A),idx);
                    set_of_Ts(i_percent,i_versns,i_beam_X,i_Tumr_Z,1) = db33.set_of_mua(row);
                    set_of_Ts(i_percent,i_versns,i_beam_X,i_Tumr_Z,2) = db33.set_of_mus(col);
                    % figure
                    % plot(db33.d_cntrs,squeeze(db33.set_of________I(row,col,:))), hold on; plot(db33.d_cntrs,squeeze(set_of_Is(i_percent,i_versns,i_beam_X,i_Tumr_Z,:))), hold off
                    % title(['r^2 = ',num2str(A(row,col)),', \mu_a = ',num2str(set_of_Ts(i_percent,i_versns,i_beam_X,i_Tumr_Z,1)),', \mu_s = ',num2str(set_of_Ts(i_percent,i_versns,i_beam_X,i_Tumr_Z,2))])
                    clearvars A idx row col i_a i_s

                    % first part
                    A = nan(length(db33.set_of_mua),length(db33.set_of_mus));
                    for i_a = 1:length(db33.set_of_mua)
                        for i_s = 1:length(db33.set_of_mus)
                            A(i_a,i_s) = r_squared(squeeze(db33.set_of________I(i_a,i_s,1:N1)),squeeze(set_of_Is(i_percent,i_versns,i_beam_X,i_Tumr_Z,1:N1)));
                        end
                    end
                    [~,idx] = min(1-A(:));
                    [row,col] = ind2sub(size(A),idx);
                    set_of_1s(i_percent,i_versns,i_beam_X,i_Tumr_Z,1) = db33.set_of_mua(row);
                    set_of_1s(i_percent,i_versns,i_beam_X,i_Tumr_Z,2) = db33.set_of_mus(col);
                    % figure
                    % plot(db33.d_cntrs,squeeze(db33.set_of________I(row,col,:))), hold on; plot(db33.d_cntrs,squeeze(set_of_Is(i_percent,i_versns,i_beam_X,i_Tumr_Z,:))), hold off
                    % title(['r^2 = ',num2str(A(row,col)),', \mu_a = ',num2str(set_of_1s(i_percent,i_versns,i_beam_X,i_Tumr_Z,1)),', \mu_s = ',num2str(set_of_1s(i_percent,i_versns,i_beam_X,i_Tumr_Z,2))])
                    clearvars A idx row col i_a i_s

                    % second part
                    A = nan(length(db33.set_of_mua),length(db33.set_of_mus));
                    for i_a = 1:length(db33.set_of_mua)
                        for i_s = 1:length(db33.set_of_mus)
                            A(i_a,i_s) = r_squared(squeeze(db33.set_of________I(i_a,i_s,N1+1:N2)),squeeze(set_of_Is(i_percent,i_versns,i_beam_X,i_Tumr_Z,N1+1:N2)));
                        end
                    end
                    [~,idx] = min(1-A(:));
                    [row,col] = ind2sub(size(A),idx);
                    set_of_2s(i_percent,i_versns,i_beam_X,i_Tumr_Z,1) = db33.set_of_mua(row);
                    set_of_2s(i_percent,i_versns,i_beam_X,i_Tumr_Z,2) = db33.set_of_mus(col);
                    % figure
                    % plot(db33.d_cntrs,squeeze(db33.set_of________I(row,col,:))), hold on; plot(db33.d_cntrs,squeeze(set_of_Is(i_percent,i_versns,i_beam_X,i_Tumr_Z,:))), hold off
                    % title(['r^2 = ',num2str(A(row,col)),', \mu_a = ',num2str(set_of_2s(i_percent,i_versns,i_beam_X,i_Tumr_Z,1)),', \mu_s = ',num2str(set_of_2s(i_percent,i_versns,i_beam_X,i_Tumr_Z,2))])
                    clearvars A idx row col i_a i_s

                    clearvars fun_x fun_y TheOutFun x_temp y_temp index_in x_bind y_bind
                    clearvars data the_percnt the_verson the_beam_X the_Tumr_Z the_HCT_Cn

                    for i_d = 1:size(set_of_Es,5)
                        A = squeeze(db33.set_of________I(:,:,i_d));
                        [~,idx] = min(abs(A(:)-set_of_Is(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_d)));
                        TheClosestValue = A(idx);
                        [row,col] = ind2sub(size(A),idx);
                        % disp(num2str([set_of_Is(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_d) db33.set_of________I(row,col,i_d) TheClosestValue]))
                        set_of_Es(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_d,1) = db33.set_of_mua(row);
                        set_of_Es(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_d,2) = db33.set_of_mus(col);
                        set_of_Es(i_percent,i_versns,i_beam_X,i_Tumr_Z,i_d,3) = TheClosestValue;
                        clearvars A idx TheClosestValue row col
                    end
                    clear i_d
                end
            end
        end
    end
end
clearvars i_tissue i_percent i_versns i_beam_X i_Tumr_Z i_HCT_Cn d_edges N_bins



for i_percent = 1:length(set_of_percnt)
    for i_versns = 1:length(set_of_versns)
        the_percnt = set_of_percnt(i_percent);
        the_verson = set_of_versns(i_versns);

        figure(1), subplot(3,3,(i_percent-1)*3+i_versns)
        A = [reshape(set_of_1s(i_percent,i_versns,:,:,1),1,[]);reshape(set_of_1s(i_percent,i_versns,:,:,2),1,[])]; [c,ic,ia] = unique(A.','rows');
        plot(c(:,1),c(:,2),'Marker','o','LineStyle','none','MarkerFaceColor','b','markersize',16,'DisplayName',['d < ',num2str(db33.d_cntrs(N1))],'Color','b'), hold on
        clearvars A c ia ic
        figure(1), subplot(3,3,(i_percent-1)*3+i_versns)
        A = [reshape(set_of_2s(i_percent,i_versns,:,:,1),1,[]);reshape(set_of_2s(i_percent,i_versns,:,:,2),1,[])]; [c,ic,ia] = unique(A.','rows');
        plot(c(:,1),c(:,2),'Marker','^','LineStyle','none','MarkerFaceColor','r','markersize',16,'DisplayName',[num2str(db33.d_cntrs(N1)),' < d < ',num2str(db33.d_cntrs(N2))],'Color','r'), hold on
        clearvars A c ia ic
        figure(1), subplot(3,3,(i_percent-1)*3+i_versns)
        A = [reshape(set_of_Ts(i_percent,i_versns,:,:,1),1,[]);reshape(set_of_Ts(i_percent,i_versns,:,:,2),1,[])]; [c,ic,ia] = unique(A.','rows');
        plot(c(:,1),c(:,2),'Marker','v','LineStyle','none','MarkerFaceColor','g','markersize',16,'DisplayName',['as a whole'],'Color','g'), hold on
        clearvars A c ia ic

        xlabel('\mu_a'), ylabel('\mu_s'), title(['Pd = ',num2str(the_percnt),', v = ',num2str(the_verson)])
        set(gca,'fontsize',12), axis tight, grid on; axis([0.0 0.5 0 500])
        legend('show','Location','northeast','NumColumns',1)
        clearvars i_fig index_in x_temp y_temp x_bind y_bind fun_x fun_y TheOutFun x_label y_label
        clearvars data the_percnt the_verson the_HCT_Cn the_beam_X the_Tumr_Z
    end
end
clearvars i_tissue i_percent i_versns i_beam_X i_Tumr_Z i_HCT_Cn d_edges N_bins



for i_percent = 1:length(set_of_percnt)
    for i_versns = 1:length(set_of_versns)
        the_percnt = set_of_percnt(i_percent);
        the_verson = set_of_versns(i_versns);

        figure(2)
        A = [reshape(set_of_1s(i_percent,i_versns,:,:,1),1,[]);reshape(set_of_1s(i_percent,i_versns,:,:,2),1,[])]; [c,ic,ia] = unique(A.','rows');
        plot(c(:,1),c(:,2),'Marker','o','LineStyle','none','MarkerFaceColor','b','markersize',16,'DisplayName',['d < ',num2str(db33.d_cntrs(N1))],'Color','b'), hold on
        clearvars A c ia ic
        figure(2)
        A = [reshape(set_of_2s(i_percent,i_versns,:,:,1),1,[]);reshape(set_of_2s(i_percent,i_versns,:,:,2),1,[])]; [c,ic,ia] = unique(A.','rows');
        plot(c(:,1),c(:,2),'Marker','^','LineStyle','none','MarkerFaceColor','r','markersize',16,'DisplayName',[num2str(db33.d_cntrs(N1)),' < d < ',num2str(db33.d_cntrs(N2))],'Color','r'), hold on
        clearvars A c ia ic
        figure(2)
        A = [reshape(set_of_Ts(i_percent,i_versns,:,:,1),1,[]);reshape(set_of_Ts(i_percent,i_versns,:,:,2),1,[])]; [c,ic,ia] = unique(A.','rows');
        plot(c(:,1),c(:,2),'Marker','v','LineStyle','none','MarkerFaceColor','g','markersize',16,'DisplayName',['as a whole'],'Color','g'), hold on
        clearvars A c ia ic

        xlabel('\mu_a'), ylabel('\mu_s'), title(['All PD''s, versions'])
        set(gca,'fontsize',12), axis tight, grid on; axis([0.0 0.5 0 500])
        legend('show','Location','southwest','NumColumns',3)
        clearvars i_fig index_in x_temp y_temp x_bind y_bind fun_x fun_y TheOutFun x_label y_label
        clearvars data the_percnt the_verson the_HCT_Cn the_beam_X the_Tumr_Z
    end
end
clearvars i_tissue i_percent i_versns i_beam_X i_Tumr_Z i_HCT_Cn d_edges N_bins



fit_data = [];
for i_percent = 1:length(set_of_percnt)
    for i_versns = 1:length(set_of_versns)
        for i_beam_X = 1:length(set_of_beam_X)
            the_percnt = set_of_percnt(i_percent);
            the_verson = set_of_versns(i_versns);
            the_beam_X = set_of_beam_X(i_beam_X);

            figure(3)
            A = [reshape(set_of_1s(i_percent,i_versns,i_beam_X,:,1),1,[]);reshape(set_of_1s(i_percent,i_versns,i_beam_X,:,2),1,[])]; [c,ic,ia] = unique(A.','rows'); fit_data = [ fit_data ; c];
            for icx = 1:size(c,1)
                plot(c(icx,1),c(icx,2),'Marker',get_marker(i_beam_X),'LineStyle','none','MarkerFaceColor','b','markersize',ic(icx)*7,'DisplayName',['d < ',num2str(db33.d_cntrs(N1)),', beamX = ',num2str(the_beam_X)],'Color','b'), hold on
            end
            clearvars A c ia ic icx

            figure(3)
            A = [reshape(set_of_2s(i_percent,i_versns,i_beam_X,:,1),1,[]);reshape(set_of_2s(i_percent,i_versns,i_beam_X,:,2),1,[])]; [c,ic,ia] = unique(A.','rows'); fit_data = [ fit_data ; c];
            for icx = 1:size(c,1)
                plot(c(icx,1),c(icx,2),'Marker',get_marker(i_beam_X),'LineStyle','none','MarkerFaceColor','r','markersize',ic(icx)*7,'DisplayName',[num2str(db33.d_cntrs(N1)),' < d < ',num2str(db33.d_cntrs(N2)),', beamX = ',num2str(the_beam_X)],'Color','r'), hold on
            end
            clearvars A c ia ic

            figure(3)
            A = [reshape(set_of_Ts(i_percent,i_versns,i_beam_X,:,1),1,[]);reshape(set_of_Ts(i_percent,i_versns,i_beam_X,:,2),1,[])]; [c,ic,ia] = unique(A.','rows'); fit_data = [ fit_data ; c];
            for icx = 1:size(c,1)
                plot(c(icx,1),c(icx,2),'Marker',get_marker(i_beam_X),'LineStyle','none','MarkerFaceColor','g','markersize',ic(icx)*7,'DisplayName',['as a whole',', beamX = ',num2str(the_beam_X)],'Color','g'), hold on
            end
            clearvars A c ia ic icx

            xlabel('\mu_a'), ylabel('\mu_s'), title(['All PD''s, versions'])
            set(gca,'fontsize',14), axis tight, grid on; axis([0.0 0.5 0 500])
            legend('show','Location','southwest','NumColumns',5,'fontsize',10)
            clearvars i_fig index_in x_temp y_temp x_bind y_bind fun_x fun_y TheOutFun x_label y_label
            clearvars data the_percnt the_verson the_HCT_Cn the_beam_X the_Tumr_Z
        end
    end
end
[params] = fit_rational_model(fit_data(:,1), fit_data(:,2));
plot(db33.set_of_mua,(params(1) + params(2)*db33.set_of_mua)./(params(3) + params(4)*db33.set_of_mua),'k')
clearvars i_tissue i_percent i_versns i_beam_X i_Tumr_Z i_HCT_Cn d_edges N_bins
clearvars fit_data


for i_percent = 1:length(set_of_percnt)
    for i_versns = 1:length(set_of_versns)
        for i_beam_X = 1:length(set_of_beam_X)
            figure(3+i_beam_X)
            for i_Tumr_Z = 1:length(set_of_Tumr_Z)
                the_percnt = set_of_percnt(i_percent);
                the_verson = set_of_versns(i_versns);
                the_beam_X = set_of_beam_X(i_beam_X);
                the_Tumr_Z = set_of_beam_X(i_Tumr_Z);

                subplot(3,3,(i_percent-1)*3+i_versns)
                i_a = db33.set_of_mua == set_of_1s(i_percent,i_versns,i_beam_X,i_Tumr_Z,1);
                i_s = db33.set_of_mus == set_of_1s(i_percent,i_versns,i_beam_X,i_Tumr_Z,2);
                plot(db33.d_cntrs,squeeze(db33.set_of________I(i_a,i_s,:)),                'LineWidth',1,'Color','b','LineStyle','-','DisplayName',['\mu_a = ',num2str(db33.set_of_mua(i_a)),', \mu_s = ',num2str(db33.set_of_mus(i_s))]),  hold on
                plot(db33.d_cntrs,squeeze(set_of_Is(i_percent,i_versns,i_beam_X,i_Tumr_Z,:)),'LineWidth',1,'Color','r','LineStyle','-.','HandleVisibility','off'), hold on

                xlabel('\mu_a'), ylabel('\mu_s'), title(['1^{st} part: Pd = ',num2str(the_percnt),', v = ',num2str(the_verson),', X = ',num2str(the_beam_X)])
                set(gca,'fontsize',12), axis tight, grid on; axis([0 5 0 20])
                legend('show','Location','southeast','NumColumns',1)
                clearvars i_fig index_in x_temp y_temp x_bind y_bind fun_x fun_y TheOutFun x_label y_label i_a i_s
                clearvars data the_percnt the_verson the_HCT_Cn the_beam_X the_Tumr_Z
            end
        end
    end
end
clearvars i_tissue i_percent i_versns i_beam_X i_Tumr_Z i_HCT_Cn d_edges N_bins



for i_percent = 1:length(set_of_percnt)
    for i_versns = 1:length(set_of_versns)
        for i_beam_X = 1:length(set_of_beam_X)
            figure(7+i_beam_X)
            for i_Tumr_Z = 1:length(set_of_Tumr_Z)
                the_percnt = set_of_percnt(i_percent);
                the_verson = set_of_versns(i_versns);
                the_beam_X = set_of_beam_X(i_beam_X);
                the_Tumr_Z = set_of_beam_X(i_Tumr_Z);

                subplot(3,3,(i_percent-1)*3+i_versns)
                i_a = db33.set_of_mua == set_of_2s(i_percent,i_versns,i_beam_X,i_Tumr_Z,1);
                i_s = db33.set_of_mus == set_of_2s(i_percent,i_versns,i_beam_X,i_Tumr_Z,2);
                plot(db33.d_cntrs,squeeze(db33.set_of________I(i_a,i_s,:)),                'LineWidth',1,'Color','b','LineStyle','-','DisplayName',['\mu_a = ',num2str(db33.set_of_mua(i_a)),', \mu_s = ',num2str(db33.set_of_mus(i_s))]),  hold on
                plot(db33.d_cntrs,squeeze(set_of_Is(i_percent,i_versns,i_beam_X,i_Tumr_Z,:)),'LineWidth',1,'Color','r','LineStyle','-.','HandleVisibility','off'), hold on

                xlabel('\mu_a'), ylabel('\mu_s'), title(['2^{nd} part: Pd = ',num2str(the_percnt),', v = ',num2str(the_verson),', X = ',num2str(the_beam_X)])
                set(gca,'fontsize',12), axis tight, grid on; axis([0 5 0 20])
                legend('show','Location','southeast','NumColumns',1)
                clearvars i_fig index_in x_temp y_temp x_bind y_bind fun_x fun_y TheOutFun x_label y_label i_a i_s
                clearvars data the_percnt the_verson the_HCT_Cn the_beam_X the_Tumr_Z
            end
        end
    end
end
clearvars i_tissue i_percent i_versns i_beam_X i_Tumr_Z i_HCT_Cn d_edges N_bins



for i_percent = 1:length(set_of_percnt)
    for i_versns = 1:length(set_of_versns)
        for i_beam_X = 1:length(set_of_beam_X)
            figure(11+i_beam_X)
            for i_Tumr_Z = 1:length(set_of_Tumr_Z)
                the_percnt = set_of_percnt(i_percent);
                the_verson = set_of_versns(i_versns);
                the_beam_X = set_of_beam_X(i_beam_X);
                the_Tumr_Z = set_of_beam_X(i_Tumr_Z);

                subplot(3,3,(i_percent-1)*3+i_versns)
                i_a = db33.set_of_mua == set_of_Ts(i_percent,i_versns,i_beam_X,i_Tumr_Z,1);
                i_s = db33.set_of_mus == set_of_Ts(i_percent,i_versns,i_beam_X,i_Tumr_Z,2);
                plot(db33.d_cntrs,squeeze(db33.set_of________I(i_a,i_s,:)),                'LineWidth',1,'Color','b','LineStyle','-','DisplayName',['\mu_a = ',num2str(db33.set_of_mua(i_a)),', \mu_s = ',num2str(db33.set_of_mus(i_s))]),  hold on
                plot(db33.d_cntrs,squeeze(set_of_Is(i_percent,i_versns,i_beam_X,i_Tumr_Z,:)),'LineWidth',1,'Color','r','LineStyle','-.','HandleVisibility','off'), hold on

                xlabel('\mu_a'), ylabel('\mu_s'), title(['whole: Pd = ',num2str(the_percnt),', v = ',num2str(the_verson),', X = ',num2str(the_beam_X)])
                set(gca,'fontsize',12), axis tight, grid on; axis([0 5 0 20])
                legend('show','Location','southeast','NumColumns',1)
                clearvars i_fig index_in x_temp y_temp x_bind y_bind fun_x fun_y TheOutFun x_label y_label i_a i_s
                clearvars data the_percnt the_verson the_HCT_Cn the_beam_X the_Tumr_Z
            end
        end
    end
end
clearvars i_tissue i_percent i_versns i_beam_X i_Tumr_Z i_HCT_Cn d_edges N_bins



for i_percent = 1:length(set_of_percnt)
    for i_versns = 1:length(set_of_versns)
        for i_beam_X = 1:length(set_of_beam_X)
            figure(15+i_beam_X)
            for i_Tumr_Z = 1:length(set_of_Tumr_Z)
                the_percnt = set_of_percnt(i_percent);
                the_verson = set_of_versns(i_versns);
                the_beam_X = set_of_beam_X(i_beam_X);
                the_Tumr_Z = set_of_beam_X(i_Tumr_Z);

                subplot(3,3,(i_percent-1)*3+i_versns)
                plot(db33.d_cntrs,squeeze(set_of_Is(i_percent,i_versns,i_beam_X,i_Tumr_Z,:)),  'LineWidth',1,'Color','b','LineStyle','-','DisplayName','original'),  hold on
                plot(db33.d_cntrs,squeeze(set_of_Es(i_percent,i_versns,i_beam_X,i_Tumr_Z,:,3)),'LineWidth',1,'Color','r','LineStyle','-.','HandleVisibility','off'), hold on

                xlabel('\mu_a'), ylabel('\mu_s'), title(['whole: Pd = ',num2str(the_percnt),', v = ',num2str(the_verson),', X = ',num2str(the_beam_X)])
                set(gca,'fontsize',12), axis tight, grid on; axis([0 5 0 20])
                legend('show','Location','southeast','NumColumns',1)
                clearvars i_fig index_in x_temp y_temp x_bind y_bind fun_x fun_y TheOutFun x_label y_label i_a i_s
                clearvars data the_percnt the_verson the_HCT_Cn the_beam_X the_Tumr_Z
            end
        end
    end
end
clearvars i_tissue i_percent i_versns i_beam_X i_Tumr_Z i_HCT_Cn d_edges N_bins


set_of_Is_analytic = nan(size(db33.set_of________I));
for i_d = 1:length(db33.d_cntrs)
    subplot(1,2,1)
    figure(20), contourf(db33.set_of_mua,db33.set_of_mus,squeeze(db33.set_of________I(:,:,i_d)),0:19), clim([0 19]), shading interp
    colormap jet, colorbar
    xlabel('\mu_a'), xlim([min(db33.set_of_mua) max(db33.set_of_mua)])
    ylabel('\mu_s'), ylim([min(db33.set_of_mus) max(db33.set_of_mus)])
    set(gca,'fontsize',18), title(['simulation: d = ',num2str(db33.d_cntrs(i_d))])
    axis square

    n = 1.4; g = 0.95;
    Set_I = nan(length(db33.set_of_mua),length(db33.set_of_mus));
    for ia = 1:length(db33.set_of_mua)
        for is = 1:length(db33.set_of_mus)
            Set_I(ia,is) = R_diffuse(db33.d_cntrs(i_d), db33.set_of_mua(ia), db33.set_of_mus(is)*(1-g), n);
            set_of_Is_analytic(ia,is,i_d) = Set_I(ia,is);
        end
    end
    figure(20), subplot(1,2,2)
    contourf(db33.set_of_mua,db33.set_of_mus,-log(Set_I).',0:19); clim([0 19]), shading interp
    colormap jet, colorbar
    xlabel('\mu_a'), xlim([min(db33.set_of_mua) max(db33.set_of_mua)])
    ylabel('\mu_s'), ylim([min(db33.set_of_mus) max(db33.set_of_mus)])
    set(gca,'fontsize',18), title(['analytical: d = ',num2str(db33.d_cntrs(i_d))])
    axis square

    clearvars n g Set_I ia is
    pause(0.1)
end


for ia = 1:2:length(db33.set_of_mua)
    for is = 1:2:length(db33.set_of_mus)
        figure
        plot(db33.d_cntrs,squeeze(db33.set_of________I(ia,is,:)),'b','DisplayName','simulation'), hold on
        plot(db33.d_cntrs,-log(squeeze(set_of_Is_analytic(ia,is,:))),'r','DisplayName','analytical'), hold on
        xlabel('d'), xlim([min(db33.d_edges) max(db33.d_edges)])
        ylabel('OD (a.u.)'), ylim([0 19])
        set(gca,'fontsize',18), title(['\mu_a = ',num2str(db33.set_of_mua(ia)),', \mu_s = ',num2str(db33.set_of_mus(is))])
        axis square

    end
end
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
function [R2] = r_squared (y_true,y_pred)
y_true = y_true(:);
y_pred = y_pred(:);
idx_in = ~(isnan(y_true)|isnan(y_pred));
y_true = y_true(idx_in);
y_pred = y_pred(idx_in);
ss_tot = sum((y_true - mean(y_true)).^2);
ss_res = sum((y_true - y_pred).^2);
R2 = 1 - (ss_res / ss_tot);
end
function [mrkr] = get_marker(i_beamX)
if     i_beamX==1
    mrkr = 'o';
elseif i_beamX==2
    mrkr = 's';
elseif i_beamX==3
    mrkr = '^';
elseif i_beamX==4
    mrkr = 'v';
else
    mrkr = 'd';
end
end
function [params] = fit_rational_model(x, y)
% Initial guess
init_params = [1, 1, 1, 1];

% Define the model function
model_fun = @(p, x) (p(1) + p(2)*x) ./ (p(3) + p(4)*x);

% Fit using lsqcurvefit (requires Optimization Toolbox)
options = optimoptions('lsqcurvefit','Display','off');
params = lsqcurvefit(model_fun, init_params, x, y, [], [], options);
end
function R = R_diffuse(rho, mua, musp, n)
% R_diffuse: Computes diffuse reflectance R(rho) at tissue surface
% under diffusion approximation for a pencil beam source.
%
% Inputs:
%   rho  - radial distance(s) [cm] (scalar or vector)
%   mua  - absorption coefficient [cm^-1]
%   musp - reduced scattering coefficient [cm^-1]
%   n    - refractive index of the medium (e.g., 1.4 for tissue)
%
% Output:
%   R    - diffuse reflectance at each rho [1/cm^2]

% Diffusion coefficient
D = 1./(3.*(mua+musp));

% Effective attenuation coefficient
mueff = sqrt(3.*mua.*(mua+musp));

% Source depth in diffusion approximation
z0 = 1./(mua+musp);

% Effective reflection coefficient (Farrell 1992 approx)
Reff = -1.440./n.^2 + 0.710./n + 0.668 + 0.0636.*n;

% Extrapolated boundary parameter A(n)
A = (1 + Reff)./(1 - Reff);

% Extrapolated boundary distance
zb = 2.*A.*D;

% Distance from source + boundary (in denominator and exponential)
r1 = sqrt(rho.^2 + (z0 + zb).^2);
disp(num2str([r1 rho]))


% Full expression for diffuse reflectance
R = (1 ./ (4 * pi * D)) .* (z0 .* (z0 + 2 * zb)) ./ (r1.^3) .* exp(-mueff .* r1);
end
