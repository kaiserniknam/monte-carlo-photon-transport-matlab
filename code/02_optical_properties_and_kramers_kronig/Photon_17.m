function [] = Photon_17 ()
% Repository group: 02_optical_properties_and_kramers_kronig
% Version role: primary simulation or study entry point
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: plot mua & mus vs lambda
clc
close all

db = load('DB/breast_db.mat');
set_of_refs = unique(db.lst(:,1));
set_of_case = unique(db.lst(:,2));
for i_ref = 1:length(set_of_refs)
    for i_case = 1:length(set_of_case)
        lmd = db.lst(db.lst(:,1)==set_of_refs(i_ref)&db.lst(:,2)==set_of_case(i_case),3);
        mua = db.lst(db.lst(:,1)==set_of_refs(i_ref)&db.lst(:,2)==set_of_case(i_case),4);
        mus = db.lst(db.lst(:,1)==set_of_refs(i_ref)&db.lst(:,2)==set_of_case(i_case),5).*(1-db.lst(db.lst(:,1)==i_ref&db.lst(:,2)==set_of_case(i_case),7));
        mup = db.lst(db.lst(:,1)==set_of_refs(i_ref)&db.lst(:,2)==set_of_case(i_case),6);
        mus(isnan(mus)) = mup(isnan(mus)); clearvars mup
        if isnan(sum(mus)), keyboard; end

        subplot(1,2,1)
        semilogy(lmd,mua,'Marker',get_marker(i_ref),'MarkerFaceColor',get_color(i_case),'MarkerEdgeColor','k','DisplayName',['ref #',num2str(i_ref),', ',get_title(i_case),' breast'],'MarkerSize',16,'LineStyle','none'), hold on
        xlabel('wavelength (\lambda, nm)'), ylabel('\mu_a (cm^{-1})'), title('\mu_a vs. wavelength')
        legend('show','Location','best')
        axis tight, axis square, set(gca,'fontsize',18)
        subplot(1,2,2)
        semilogy(lmd,mus,'Marker',get_marker(i_ref),'MarkerFaceColor',get_color(i_case),'MarkerEdgeColor','k','DisplayName',['ref #',num2str(i_ref),', ',get_title(i_case),' breast'],'MarkerSize',16,'LineStyle','none'), hold on
        xlabel('wavelength (\lambda, nm)'), ylabel('\mu_s'' (cm^{-1})'), title('\mu_s'' vs. wavelength')
        legend('show','Location','best')
        axis tight, axis square, set(gca,'fontsize',18)
        clearvars mua mus lmd
    end
end
end
function [out] = get_marker(idx)
if     idx == 1
    out = 'o';
elseif idx == 2
    out = 's';
elseif idx == 3
    out = 'h';
elseif idx == 4
    out = 'd';
elseif idx == 5
    out = '^';
elseif idx == 6
    out = 'v';
else
    keyboard
end
end
function [out] = get_color(idx)
if     idx == 1
    out = 'b';
elseif idx == 2
    out = 'r';
else
    keyboard
end
end
function [out] = get_title(idx)
if     idx == 1
    out = 'healthy';
elseif idx == 2
    out = 'cancer';
else
    keyboard
end
end
