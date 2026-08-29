function [] = Photon_46 ()
% Repository group: 06_laboratory_phantom_models
% Version role: primary simulation or study entry point
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: Reading Arnab data - April 18, 2025

clc
close all

TheData = nan(4,3,3,3); % Data matrix: [4 Tartazine (background dye) concentrations X 3 Indian Ink (signal) concentrations × 3 TiO2 concentrations × 3 measurements]
set_of_Trzn = [0.0,   0.1,  0.3,  0.6];
set_of_Iink = [0.01, 0.05, 0.10];
set_of_TiO2 = [0.5, 1.0, 2.0];

TheData(1,1,:,:) = [0.570800    0.56980     0.57120
                    0.196550    0.20105	    0.19404
                    0.159945	0.16201	    0.16056]; % Tartazine = 0.0, Indian Ink = 0.01;
TheData(1,2,:,:) = [0.530400	0.52014	    0.54160
                    0.174820	0.18540 	0.17894
                    0.140599	0.14821	    0.13270]; % Tartazine = 0.0, Indian Ink = 0.05;
TheData(1,3,:,:) = [0.441200 	0.39810	    0.43100
                    0.141500    0.14121	    0.14960
                    0.110295	0.10867	    0.11592]; % Tartazine = 0.0, Indian Ink = 0.10;

TheData(2,1,:,:) = [0.54280	    0.53946 	0.54198
                    0.18283	    0.17956	    0.18898
                    0.14862	    0.14629	    0.14967]; % Tartazine = 0.1, Indian Ink = 0.01;
TheData(2,2,:,:) = [0.49193	    0.48921 	0.50243
                    0.15153	    0.15287	    0.15058
                    0.12498	    0.12895	    0.11998]; % Tartazine = 0.1, Indian Ink = 0.05;
TheData(2,3,:,:) = [0.35350	    0.36052 	0.34935
                    0.12433	    0.12860	    0.11934
                    0.08511	    0.09016 	0.08546]; % Tartazine = 0.1, Indian Ink = 0.10;

TheData(3,1,:,:) = [0.52150	    0.49864	    0.51460
                    0.17543	    0.17881	    0.16957
                    0.13844 	0.12967	    0.13597]; % Tartazine = 0.3, Indian Ink = 0.01;
TheData(3,2,:,:) = [0.44080	    0.42680     0.44590
                    0.14153	    0.14762	    0.15160
                    0.11649	    0.10970 	0.11297]; % Tartazine = 0.3, Indian Ink = 0.05;
TheData(3,3,:,:) = [0.31540	    0.29460 	0.31496
                    0.11433	    0.12160 	0.11584
                    0.07551	    0.07895 	0.06989]; % Tartazine = 0.3, Indian Ink = 0.10;

TheData(4,1,:,:) = [0.57231	    0.56897	    0.574600
                    0.19567	    0.19642	    0.192841
                    0.16269	    0.16569	    0.159670]; % Tartazine = 0.6, Indian Ink = 0.01;
TheData(4,2,:,:) = [0.53408	    0.52980 	0.531590
                    0.174153	0.17476 	0.175160
                    0.140379	0.13997	    0.138297]; % Tartazine = 0.6, Indian Ink = 0.05;
TheData(4,3,:,:) = [0.421540	0.42946 	0.421496
                    0.150330	0.14916	    0.145840
                    0.107551	0.107895	0.116989]; % Tartazine = 0.6, Indian Ink = 0.10;


TheBase = nan(3,1); % Data matrix: [4 Tartazine (background dye) concentrations X 3 Indian Ink (signal) concentrations × 3 TiO2 concentrations × 3 measurements]
TheBase(:,:) = [0.58670
                0.20883
                0.16476];


figure(1)
for i_TiO2 = 1:3
    plot(set_of_Iink,mean(TheData(1,:,i_TiO2,:),4),'LineWidth',2,'LineStyle','-', 'Color',get_color(i_TiO2),'DisplayName',['TiO_2 ',num2str(set_of_TiO2(i_TiO2)),' g ',num2str(set_of_Trzn(1)),' dye'],'Marker','s','MarkerSize',16), hold on
    plot(set_of_Iink,mean(TheData(2,:,i_TiO2,:),4),'LineWidth',2,'LineStyle','--','Color',get_color(i_TiO2),'DisplayName',['TiO_2 ',num2str(set_of_TiO2(i_TiO2)),' g ',num2str(set_of_Trzn(4)),' dye'],'Marker','s','MarkerSize',16), hold on
end
xlabel('Ink Concentration (%)')
ylabel('Signal Intensity (mV)')
set(gca,'FontSize',24)
axis([0.01 0.1 0.0 0.6 ]), grid on, legend('show','Location','northeast','NumColumns',1)
clearvars i_TiO2

figure(2)
idx = 1; subplot(1,3,idx)
bar(set_of_Iink,TheBase(idx)-squeeze(mean(TheData([1,3],:,idx,:),4)).');
xlabel('Ink Concentration'), set(gca,'xtick',[0.01,0.05,0.10]), grid on
ylabel('Signal Strength (mV)'), ylim([0.00 0.32]), title([num2str(set_of_TiO2(idx)),'g TiO_2'])
set(gca,'fontsize',20)
clearvars idx
idx = 2; subplot(1,3,idx)
bar(set_of_Iink,TheBase(idx)-squeeze(mean(TheData([1,3],:,idx,:),4)).');
xlabel('Ink Concentration'), set(gca,'xtick',[0.01,0.05,0.10]), grid on
ylabel('Signal Strength (mV)'), ylim([0.00 0.12]), title([num2str(set_of_TiO2(idx)),'g TiO_2'])
set(gca,'fontsize',20)
clearvars idx
idx = 3; subplot(1,3,idx)
bar(set_of_Iink,TheBase(idx)-squeeze(mean(TheData([1,3],:,idx,:),4)).');
xlabel('Ink Concentration'), set(gca,'xtick',[0.01,0.05,0.10]), grid on
ylabel('Signal Strength (mV)'), ylim([0.00 0.14]), title([num2str(set_of_TiO2(idx)),'g TiO_2'])
set(gca,'fontsize',20)
clearvars idx

figure(3)
for i_TiO2 = 1:3
    for i_dye = 2:4
        y_temp = squeeze(mean(TheData([1,i_dye],:,i_TiO2,:),4));
        y_temp(1,:) = TheBase(i_TiO2) - y_temp(1,:);
        y_temp(2,:) = TheBase(i_TiO2) - y_temp(2,:);
        subplot(3,3,i_TiO2+3*(i_dye-2))
        bar(set_of_Iink,y_temp)
        if i_dye==4,  xlabel('Ink Concentration');    else; xlabel(''); end
        if i_TiO2==1, ylabel('signal strength (mV)'); else; ylabel(''); end
        title([num2str(set_of_TiO2(i_TiO2)),' TiO_2 - ',num2str(set_of_Trzn(i_dye)),' M'])
        set(gca,'FontSize',12)
    end
end
end

function [out] = get_color(i_color)
if     i_color==1
    out = 'b';
elseif i_color==2
    out = 'r';
elseif i_color==3
    out = 'g';
else
    out = 'k';
end
end
