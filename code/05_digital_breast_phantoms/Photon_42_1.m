function [] = Photon_42_1 ()
% Repository group: 05_digital_breast_phantoms
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: simulating phantoms prepared by Diego
% analyzing simulated data

clc
close all
code_num = 42;

% Optical properties of breast and tumor [1-4]
% mua_glnd = 0.55; mus_glnd = 332.7; g_glnd = 0.965; n_glnd = 1.4;
% mua_adps = 0.82; mus_adps = 313.8; g_adps = 0.976; n_adps = 1.4;
set_of_percnt = [25,50,75]; % percent of adipose
set_of_versns = [0,1,2]; % version of each density
set_of_endmat = nan(length(set_of_percnt),length(set_of_versns),3);
set_of_endmap = nan(length(set_of_percnt),length(set_of_versns),510,162);

for i_percent = 1:length(set_of_percnt)
    for i_versns = 1:length(set_of_versns)
        the_version = set_of_versns(i_versns);
        the_percent = set_of_percnt(i_percent);
        data = load(['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_42_Pd',num2str(the_percent),'_',num2str(the_version),'_y_-z_x.mat']);

        % % removing noisy points
        % idx = -10<=data.p_ot(:,1)&data.p_ot(:,1)<=+10&-3.0<=data.p_ot(:,2)&data.p_ot(:,2)<=+3.0;
        % data.p_in = data.p_in(idx,:);
        % data.i_in = data.i_in(idx,:);
        % data.m_in = data.m_in(idx,:);
        % data.w_in = data.w_in(idx,:);
        % data.p_ot = data.p_ot(idx,:);
        % data.i_ot = data.i_ot(idx,:);
        % data.m_ot = data.m_ot(idx,:);
        % data.w_ot = data.w_ot(idx,:);
        % data.s    = data.s(idx,:)   ;
        % disp(['passed photons = ',num2str(sum(idx)/data.no_of_photons*100,'%.2f')])
        % data.no_of_photons = sum(idx);
        % clearvars idx

        % Display information on the percentage of adipose tissue (1), and the percentage of photons that ended in each tissue type: air, adipose, and fibrogranular.
        set_of_endmat(i_percent,i_versns,1) = sum(data.m_ot==1)/length(data.m_ot);
        set_of_endmat(i_percent,i_versns,2) = sum(data.m_ot==2)/length(data.m_ot);
        set_of_endmat(i_percent,i_versns,3) = sum(data.m_ot==3)/length(data.m_ot);
        disp(['Pd = ',num2str(the_percent),' ~ ',num2str(100*sum(reshape(data.M_raw==2,[],1))/(sum(reshape(data.M_raw==2,[],1))+sum(reshape(data.M_raw==3,[],1))),'%.2f'),', version = ',num2str(the_version),': ( ',num2str(set_of_endmat(i_percent,i_versns,1),'%.2f'),' , ',num2str(set_of_endmat(i_percent,i_versns,2),'%.2f'),' , ',num2str(set_of_endmat(i_percent,i_versns,3),'%.2f'),' )'])

        % Distribution of the z-positions of photons that ended in each tissue type: air, adipose, and fibrogranular / boxplot
        figure(1), subplot(3,3,(i_percent-1)*3+i_versns)
        boxplot(data.p_ot(:,3),data.m_ot,'Symbol','o')
        set(gca,'fontsize',12), axis tight, grid on
        if i_percent==1, title(['version = ',num2str(the_version)]); else, title("");end
        if i_percent==length(set_of_percnt), set(gca,'xtick',1:3); set(gca,'xticklabel',{'air','fat','fibre'}); else, set(gca,'xtick',1:3); set(gca,'xticklabel',{'','',''}); end
        set(gca,'ytick',0:1:4), ylim([0 4])
        if i_versns==1, ylabel([{['Pd = ',num2str(the_percent)]},{'ending z (cm)'}]); else, ylabel(""); end
        xlim([0 +4])

        % Distribution of the z-positions of photons that ended in each tissue type: air, adipose, and fibrogranular / histogram
        figure(2), subplot(3,3,(i_percent-1)*3+i_versns)
        % ot = get_point(data.p_ot,data.m_ot,1); histogram(ot(:,3),'FaceColor','b','EdgeColor','k'), hold on, clearvars ot
        ot = get_point(data.p_ot,data.m_ot,2); histogram(ot(:,3),'FaceColor','y','EdgeColor','k'), hold on, clearvars ot
        ot = get_point(data.p_ot,data.m_ot,3); histogram(ot(:,3),'FaceColor','r','EdgeColor','k'), hold on, clearvars ot
        set(gca,'fontsize',12), axis tight, grid on
        if i_percent==1, title(['version = ',num2str(the_version)]); else, title("");end
        if i_percent==length(set_of_percnt), xlabel('ending z (cm)'); else, xlabel(''); end
        if i_versns==1, ylabel([{['Pd = ',num2str(the_percent)]},{'freq'}]); else, ylabel(""); end

        % Distribution of the z-positions of photons that ended in each tissue type: air, adipose, and fibrogranular
        figure(3), subplot(3,3,(i_percent-1)*3+i_versns)
        ot = data.s(data.m_ot==1); histogram(ot,'FaceColor','b','EdgeColor','none'), hold on, clearvars ot
        ot = data.s(data.m_ot==2); histogram(ot,'FaceColor','y','EdgeColor','none'), hold on, clearvars ot
        ot = data.s(data.m_ot==3); histogram(ot,'FaceColor','r','EdgeColor','none'), hold on, clearvars ot
        set(gca,'fontsize',12), axis tight, grid on
        if i_percent==1, title(['version = ',num2str(the_version)]); else, title("");end
        if i_percent==length(set_of_percnt), xlabel('s (cm)'); else, xlabel(''); end
        if i_versns==1, ylabel([{['Pd = ',num2str(the_percent)]},{'freq'}]); else, ylabel(""); end

        % Distribution of OD
        i_m = 1;
        N_lines = 38; [set_of_lines,my_colormap] = make_colormap(0,18,N_lines);
        x_sg = data.p_ot(data.m_ot==i_m,1);
        y_sg = data.p_ot(data.m_ot==i_m,2);
        u_sg = data.w_ot(data.m_ot==i_m);
        n_sg = data.no_of_photons; Nx = size(data.M_raw,1); Ny = size(data.M_raw,2);
        [t_db,x_c,y_c] = do_sum (x_sg,y_sg,u_sg,Nx*data.dl,Ny*data.dl,Nx,Ny,@(w,d,n)(-log(sum(w)./n)),n_sg,0,0);
        set_of_endmap(i_percent,i_versns,:,:) = t_db;
        clearvars i_m x_sg y_sg n_sg u_sg Nx Ny
        figure(4), subplot(3,3,(i_percent-1)*3+i_versns)
        contourf(x_c, y_c, t_db.', set_of_lines,'LineColor', 'none', 'LineWidth', 0.1), clearvars t_db
        set(gca,'fontsize',12), axis equal, axis tight, grid on; h = colorbar; % ylabel(h,'OD/d-OD/d_{not-comp} (cm^{-1})')
        colormap('jet');  clim([min(set_of_lines) max(set_of_lines)])
        set(gca,'fontsize',12), axis tight, grid on
        if i_percent==1, title(['version = ',num2str(the_version)]); else, title("");end
        if i_percent==length(set_of_percnt), xlabel('x (cm)'); else, xlabel(''); end
        if i_versns==1, ylabel([{['Pd = ',num2str(the_percent)]},{'y (cm)'}]); else, ylabel(""); end
        axis([-2 +2 -2 +2])
        clearvars my_colormap h N_lines set_of_lines

        % Distribution of the z-positions of photons that ended in each tissue type: air, adipose, and fibrogranular
        figure(5), subplot(3,3,(i_percent-1)*3+i_versns)
        ot = get_point(data.p_ot,data.m_ot,1); plot3(ot(:,1),ot(:,2),-ot(:,3),'b.','MarkerSize',0.1), hold on, clearvars ot
        ot = get_point(data.p_ot,data.m_ot,2); plot3(ot(:,1),ot(:,2),-ot(:,3),'y.','MarkerSize',0.1), hold on, clearvars ot
        ot = get_point(data.p_ot,data.m_ot,3); plot3(ot(:,1),ot(:,2),-ot(:,3),'r.','MarkerSize',0.1), hold on, clearvars ot
        set(gca,'fontsize',12), axis tight, grid on
        if i_percent==1, title(['version = ',num2str(the_version)]); else, title("");end
        xlabel('x (cm)'), ylabel('y (cm)')
        if i_versns==1, ylabel([{['Pd = ',num2str(the_percent)]},{'z (cm)'}]); else, ylabel(""); end
        xlim([-3 +3]), set(gca,'xtick',-3:+1:+3)
        ylim([-3 +3]), set(gca,'ytick',-3:+1:+3)
        zlim([-3 +0]), set(gca,'ztick',-3:+1:+3)


        % Distribution of the z-positions of photons that ended in each tissue type: air, adipose, and fibrogranular
        figure(6), subplot(3,3,(i_percent-1)*3+i_versns)
        ot = get_point(data.p_ot,data.m_ot,1); plot3(ot(:,1),ot(:,2),-ot(:,3),'b.','MarkerSize',0.1), hold on, clearvars ot
        ot = get_point(data.p_ot,data.m_ot,2); plot3(ot(:,1),ot(:,2),-ot(:,3),'y.','MarkerSize',0.1), hold on, clearvars ot
        ot = get_point(data.p_ot,data.m_ot,3); plot3(ot(:,1),ot(:,2),-ot(:,3),'r.','MarkerSize',0.1), hold on, clearvars ot
        for i_x  = 1:size(data.M_raw,1)
            for i_y  = 1:size(data.M_raw,2)
                if length(unique(squeeze(data.M_raw(i_x,i_y,:)))) > 1
                    i_z = find(squeeze(data.M_raw(i_x,i_y,:))~=1,1,"first");
                    plot3((i_x-size(data.M_raw,1)/2)*data.dl,(i_y-size(data.M_raw,2)/2)*data.dl,-i_z*data.dl,'b.','MarkerSize',0.1), hold on, clearvars i_z
                end
            end
        end
        clearvars i_x i_y
        set(gca,'fontsize',12), axis tight, grid on
        if i_percent==1, title(['version = ',num2str(the_version)]); else, title("");end
        xlabel('x (cm)'), ylabel('y (cm)'), view([-180 90])
        if i_versns==1, ylabel([{['Pd = ',num2str(the_percent)]},{'z (cm)'}]); else, ylabel(""); end
        % xlim([-3 +3]), set(gca,'xtick',-3:+1:+3)
        % ylim([-3 +3]), set(gca,'ytick',-3:+1:+3)
        % zlim([-3 +0]), set(gca,'ztick',-3:+1:+3)

        clearvars the_percent the_version data
    end
end
clearvars i_tissue i_percenti_versns

% Absorption rate versus tissue type
figure(11)
the_colors = [0 0 1;1 1 0;1 0 0];
for i_percent = 1:length(set_of_percnt)
    boxplot(squeeze(set_of_endmat(i_percent,:,:)),'BoxStyle','outline','Colors',the_colors), hold on
    h = findobj(gca,'Tag','Box');
    for j = 1:3
        patch(get(h(j),'XData'),get(h(j),'YData'),get(h(j),'Color'),'FaceAlpha',set_of_percnt(i_percent)/100,'DisplayName',[get_material(j),', Pd = ',num2str(set_of_percnt(i_percent)),'%']);
    end
    clearvars h j
end
set(gca,'xtick',1:3), set(gca,'xticklabel',{'air','fat','fibre'})
ylabel('absorption rate (%)'), set(gca,'ytick',0:0.1:1), ylim([0 0.8])
legend('show','Location','northeast')
set(gca,'fontsize',18)
clearvars the_colors
clearvars i_tissue i_percenti_versns

% Absorption rate versus Pd
figure(12)
the_colors = [0 0 1;1 1 0;1 0 0];
for i_tissue = 1:3
    boxplot(squeeze(set_of_endmat(:,:,i_tissue).'),'BoxStyle','outline','Colors',[the_colors(i_tissue,:);the_colors(i_tissue,:);the_colors(i_tissue,:)]), hold on
    h = findobj(gca,'Tag','Box');
    for j = 1:length(set_of_percnt)
        patch(get(h(j),'XData'),get(h(j),'YData'),get(h(j),'Color'),'FaceAlpha',1-set_of_percnt(j)/100,'DisplayName',[get_material(i_tissue),', Pd = ',num2str(set_of_percnt(j)),'%']);
    end
    clearvars h j
end
set(gca,'xtick',1:3), set(gca,'xticklabel',{'Pd = 25','Pd = 50','Pd = 75'})
ylabel('absorption rate (%)'), set(gca,'ytick',0:0.1:1), ylim([0 0.8])
legend('show','Location','west')
set(gca,'fontsize',18)
clearvars the_colors i_tissue

% D_OD: OD/d-OD/d_{Pd=1}
figure(13)
set_of_lines = -1:0.1:+1;
for i_percent = 1:length(set_of_percnt)
    for i_versns = 1:length(set_of_versns)
        the_version = set_of_versns(i_versns);
        the_percent = set_of_percnt(i_percent);
        subplot(3,3,(i_percent-1)*3+i_versns)
        contourf(x_c, y_c, squeeze(set_of_endmap(i_percent,i_versns,:,:)).'-squeeze(set_of_endmap(1,i_versns,:,:)).', set_of_lines, 'LineColor', 'none', 'LineWidth', 0.1)
        set(gca,'fontsize',12), axis equal, axis tight, grid on; h = colorbar; ylabel(h,'OD/d-OD/d_{Pd=25} (cm^{-1})')
        colormap('jet');  clim([min(set_of_lines) max(set_of_lines)])
        if i_percent==1, title(['version = ',num2str(the_version)]); else, title("");end
        if i_percent==length(set_of_percnt), xlabel('x (cm)'); else, xlabel(''); end
        if i_versns==1, ylabel([{['Pd = ',num2str(the_percent)]},{'y (cm)'}]); else, ylabel(""); end
        axis([-2 +2 -2 +2])
        clearvars the_percent the_version
    end
end

% D_OD: OD/d-OD/d_{dv}
figure(14)
set_of_lines = -1:0.1:+1;
for i_percent = 1:length(set_of_percnt)
    for i_versns = 1:length(set_of_versns)
        i1 = mod(i_versns+0,3)+1;
        i2 = mod(i_versns+1,3)+1;
        the_version = set_of_versns(i_versns);
        the_percent = set_of_percnt(i_percent);
        subplot(3,3,(i_percent-1)*3+i_versns)
        contourf(x_c, y_c, squeeze(set_of_endmap(i_percent,i1,:,:)).'-squeeze(set_of_endmap(i_percent,i2,:,:)).', set_of_lines, 'LineColor', 'none', 'LineWidth', 0.1)
        set(gca,'fontsize',12), axis equal, axis tight, grid on; h = colorbar; ylabel(h,['OD/d_{v=',num2str(i1),'}-OD/d_{v=',num2str(i2),'} (cm^{-1})'])
        colormap('jet');  clim([min(set_of_lines) max(set_of_lines)])
        if i_percent==1, title(['version = ',num2str(the_version)]); else, title("");end
        if i_percent==length(set_of_percnt), xlabel('x (cm)'); else, xlabel(''); end
        if i_versns==1, ylabel([{['Pd = ',num2str(the_percent)]},{'y (cm)'}]); else, ylabel(""); end
        axis([-2 +2 -2 +2])
        clearvars the_percent the_version
    end
end

for i_fig = 14:-1:1, figure(i_fig), saveas(gcf,['Photon_',num2str(code_num),'_',num2str(1),'_fig_',sprintf('%2.0f',i_fig),'.jpg' ],'jpeg'); end
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
function [out] = get_point(in,c,c_chose)
out = in(c==c_chose,:);
end
function [out] = get_material(in)
if in==1
    out = 'air';
elseif in==2
    out = 'fat';
elseif in==3
    out = 'fibre';
else
    out = 'none';
end
end
function [u_out,x_c,y_c] = do_sum  (x_sg,y_sg,u_sg,Lx,Ly,Nx_bin,Ny_bin,TheFun,n_sg,beam_X,beam_Y)
x_edges = linspace(-Lx/2,+Lx/2,Nx_bin+1); % # of x bins
y_edges = linspace(-Ly/2,+Ly/2,Ny_bin+1); % # of y bins
[~,~,~,ind_x_sg,ind_y_sg] = histcounts2(x_sg,y_sg,x_edges,y_edges);
x_c = 1/2*(x_edges(1:end-1)+x_edges(2:end-0)).';
y_c = 1/2*(y_edges(1:end-1)+y_edges(2:end-0)).';
clearvars x_edges y_edges

u_out =  nan(Nx_bin,Ny_bin);
clearvars Nx_bin Ny_bin
U = [ind_x_sg,ind_y_sg];
[unique_U,~,~] = unique(U,'rows'); clearvars U
for i_U = 1:size(unique_U,1)
    idx_from_sg = ind_x_sg==unique_U(i_U,1)&ind_y_sg==unique_U(i_U,2);
    u_sg_loop = u_sg(idx_from_sg); x_sg_loop = x_sg(idx_from_sg); y_sg_loop = y_sg(idx_from_sg); d_sg_loop = sqrt((x_sg_loop-beam_X).^2+(y_sg_loop-beam_Y).^2);
    u_out(unique_U(i_U,1),unique_U(i_U,2)) = TheFun(u_sg_loop,mean(d_sg_loop),n_sg);
    clearvars u_sg_loop x_sg_loop y_sg_loop d_sg_loop i_counter N t_sg idx_from_sg idx_from_sg
end
end
function [set_of_lines,my_colormap] = make_colormap(l_start,l_end,N_lines)
set_of_lines = linspace(l_start,l_end,2*N_lines+1);
my_colormap = [1, 1, 1];
if     l_end<=0
    for idx = 1:2*N_lines
        my_colormap = [...
            [1-idx/2/N_lines, 1-idx/2/N_lines, 1];...
            my_colormap];
    end
elseif l_start>=0
    for idx = 1:2*N_lines
        my_colormap = [...
            my_colormap;...
            [1, 1-idx/2/N_lines, 1-idx/2/N_lines]];
    end
else
    Np = sum(set_of_lines>0);
    Nn = sum(set_of_lines<0);
    for idx = 1:Np
        my_colormap = [...
            my_colormap;...
            [1, 1-idx/Np, 1-idx/Np]];
    end
    for idx = 1:Nn
        my_colormap = [...
            [1-idx/Nn, 1-idx/Nn, 1];...
            my_colormap];
    end
end
end
