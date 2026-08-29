function [] = Photon_45_3 ()
% Repository group: 05_digital_breast_phantoms
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: simulating phantoms prepared by Diego
% Similar to Photon_44: but in better intervals
% analyzing simulated data: Testing an idea

clc
close all
load('DB/Photon_45_2.mat');

% i_fig = 1; idx_d = 2:4:55; Y_Limit = [0 20]; i_beam_X = 2; idx_HCT_Cn = 2:length(set_of_HCT_Cn);
% for i_d = 1:length(idx_d)
%     TheColor = [(i_d-1)/(length(idx_d)-1),0,1-(i_d-1)/(length(idx_d)-1)];
%     for i_Tumr_Z = 1:length(set_of_Tumr_Z)
%         figure(i_fig), subplot(2,2,i_Tumr_Z)
%         plot(set_of_HCT_Cn(idx_HCT_Cn),squeeze(mean(mean(-log(set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn,idx_d(i_d))),1),2)),'Color',TheColor,'DisplayName',['d = ',num2str(d_edges(idx_d(i_d)))],'LineWidth',2), hold on
%         xlabel('HCT (%)'), xlim([min(set_of_HCT_Cn) max(set_of_HCT_Cn)])
%         ylabel('OD (a.u.)'), ylim([min(Y_Limit) max(Y_Limit)])
%         set(gca,'fontsize',12), grid on; % legend('show','numColumns',1,'Orientation','horizontal'),
%         title(['beamX = ',num2str(set_of_beam_X(i_beam_X)),' cm, Tumor depth = ',num2str(set_of_Tumr_Z(i_Tumr_Z)),' cm'])
%     end
%     clearvars TheColor
% end
% clearvars i_fig idx_d Y_Limit i_beam_X idx_HCT_Cn i_d i_Tumr_Z
%
% i_fig = 2; idx_d = 2:4:55; Y_Limit = [0 6]; i_beam_X = 2; idx_HCT_Cn = 2:length(set_of_HCT_Cn);
% for i_d = 1:length(idx_d)
%     TheColor = [(i_d-1)/(length(idx_d)-1),0,1-(i_d-1)/(length(idx_d)-1)];
%     for i_Tumr_Z = 1:length(set_of_Tumr_Z)
%         figure(i_fig), subplot(2,2,i_Tumr_Z)
%         plot(set_of_HCT_Cn(idx_HCT_Cn),squeeze(mean(mean(-log(set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn,idx_d(i_d))),1),2))./set_of_HCT_Cn(idx_HCT_Cn).','Color',TheColor,'DisplayName',['d = ',num2str(d_edges(idx_d(i_d)))],'LineWidth',2), hold on
%         xlabel('HCT (%)'), xlim([min(set_of_HCT_Cn) max(set_of_HCT_Cn)])
%         ylabel('OD/HCT (a.u.)'), ylim([min(Y_Limit) max(Y_Limit)])
%         set(gca,'fontsize',12), grid on; % legend('show','numColumns',1,'Orientation','horizontal'),
%         title(['beamX = ',num2str(set_of_beam_X(i_beam_X)),' cm, Tumor depth = ',num2str(set_of_Tumr_Z(i_Tumr_Z)),' cm'])
%     end
%     clearvars TheColor
% end
% clearvars i_fig idx_d Y_Limit i_beam_X idx_HCT_Cn i_d i_Tumr_Z
%
% i_fig = 3; idx_d = 2:4:55; Y_Limit = [0 35]; i_beam_X = 2; idx_HCT_Cn = 2:length(set_of_HCT_Cn);
% for i_d = 1:length(idx_d)
%     TheColor = [(i_d-1)/(length(idx_d)-1),0,1-(i_d-1)/(length(idx_d)-1)];
%     for i_Tumr_Z = 1:length(set_of_Tumr_Z)
%         figure(i_fig), subplot(2,2,i_Tumr_Z)
%         plot(set_of_HCT_Cn(idx_HCT_Cn),squeeze(mean(mean(-log(set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn,idx_d(i_d))),1),2))./d_edges(idx_d(i_d)).','Color',TheColor,'DisplayName',['d = ',num2str(d_edges(idx_d(i_d)))],'LineWidth',2), hold on
%         xlabel('HCT (%)'), xlim([min(set_of_HCT_Cn) max(set_of_HCT_Cn)])
%         ylabel('OD/d (cm^{-1})'), ylim([min(Y_Limit) max(Y_Limit)])
%         set(gca,'fontsize',12), grid on; % legend('show','numColumns',1,'Orientation','horizontal'),
%         title(['beamX = ',num2str(set_of_beam_X(i_beam_X)),' cm, Tumor depth = ',num2str(set_of_Tumr_Z(i_Tumr_Z)),' cm'])
%     end
%     clearvars TheColor
% end
% clearvars i_fig idx_d Y_Limit i_beam_X idx_HCT_Cn i_d i_Tumr_Z
%
% i_fig = 4; idx_d = 2:4:55; Y_Limit = [0 15]; i_beam_X = 1; idx_HCT_Cn = 2:length(set_of_HCT_Cn);
% for i_d = 1:length(idx_d)
%     TheColor = [(i_d-1)/(length(idx_d)-1),0,1-(i_d-1)/(length(idx_d)-1)];
%     for i_Tumr_Z = 1:length(set_of_Tumr_Z)
%         figure(i_fig), subplot(2,2,i_Tumr_Z)
%         plot(set_of_HCT_Cn(idx_HCT_Cn),squeeze(mean(mean(-log(set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn,idx_d(i_d))),1),2))./d_edges(idx_d(i_d))./set_of_HCT_Cn(idx_HCT_Cn).','Color',TheColor,'DisplayName',['d = ',num2str(d_edges(idx_d(i_d)))],'LineWidth',2), hold on
%         xlabel('HCT (%)'), xlim([min(set_of_HCT_Cn) max(set_of_HCT_Cn)])
%         ylabel('OD/d/HCT (cm^{-1})'), ylim([min(Y_Limit) max(Y_Limit)])
%         set(gca,'fontsize',12), grid on; % legend('show','numColumns',1,'Orientation','horizontal'),
%         title(['beamX = ',num2str(set_of_beam_X(i_beam_X)),' cm, Tumor depth = ',num2str(set_of_Tumr_Z(i_Tumr_Z)),' cm'])
%     end
%     clearvars TheColor
% end
% clearvars i_fig idx_d Y_Limit i_beam_X idx_HCT_Cn i_d i_Tumr_Z
%
% i_fig = 5; idx_d = 2:4:55; Y_Limit = [0 15]; i_beam_X = 1; idx_HCT_Cn = 2:length(set_of_HCT_Cn);
% for i_d = 1:length(idx_d)
%     for i_Tumr_Z = 1:length(set_of_Tumr_Z)
%         figure(i_fig)
%         TheColor = [(i_d-1)/(length(idx_d)-1),(i_Tumr_Z-1)./(length(set_of_Tumr_Z)-1),1-(i_d-1)/(length(idx_d)-1)];
%         plot(set_of_HCT_Cn(idx_HCT_Cn),squeeze(mean(mean(-log(set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn,idx_d(i_d))),1),2))./d_edges(idx_d(i_d))./set_of_HCT_Cn(idx_HCT_Cn).','Color',TheColor,'DisplayName',['d = ',num2str(d_edges(idx_d(i_d)))],'LineWidth',2), hold on
%         xlabel('HCT (%)'), xlim(9[min(set_of_HCT_Cn) max(set_of_HCT_Cn)])
%         ylabel('OD/d/HCT (cm^{-1})'), ylim([min(Y_Limit) max(Y_Limit)])
%         set(gca,'fontsize',12), grid on; % legend('show','numColumns',1,'Orientation','horizontal'),
%         title(['beamX = ',num2str(set_of_beam_X(i_beam_X)),' cm, Tumor depth = ',num2str(set_of_Tumr_Z(i_Tumr_Z)),' cm'])
%     end
%     clearvars TheColor
% end
% clearvars i_fig idx_d Y_Limit i_beam_X idx_HCT_Cn i_d i_Tumr_Z

i_fig = 6; idx_d = 2:4:35; Y_Limit = [0 0.8]; i_beam_X = 2; idx_HCT_Cn = 1:length(set_of_HCT_Cn);
for i_d = 1:length(idx_d)
    TheColor = [(i_d-1)/(length(idx_d)-1),0,1-(i_d-1)/(length(idx_d)-1)];
    for i_Tumr_Z = 1:length(set_of_Tumr_Z)
        figure(i_fig), subplot(2,2,i_Tumr_Z)
        plot(set_of_HCT_Cn(idx_HCT_Cn),(squeeze(mean(mean(-log(set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn,idx_d(i_d))),1),2))),'Color',TheColor,'DisplayName',['d = ',num2str(d_edges(idx_d(i_d)))],'LineWidth',2), hold on
        xlabel('HCT (%)'), xlim([min(set_of_HCT_Cn) max(set_of_HCT_Cn)])
        ylabel('OD (a.u.)'), % ylim([min(Y_Limit) max(Y_Limit)])
        set(gca,'fontsize',12), grid on; % legend('show','numColumns',1,'Orientation','horizontal'),
        title(['beamX = ',num2str(set_of_beam_X(i_beam_X)),' cm, Tumor depth = ',num2str(set_of_Tumr_Z(i_Tumr_Z)),' cm'])
    end
    clearvars TheColor
end
clearvars i_fig idx_d Y_Limit i_beam_X idx_HCT_Cn i_d i_Tumr_Z

i_fig = 7; idx_d = 2:2:35; Y_Limit = [0 0.3]; i_beam_X = 2; idx_HCT_Cn = 2:length(set_of_HCT_Cn);
for i_d = 1:length(idx_d)
    TheColor = [(i_d-1)/(length(idx_d)-1),0,1-(i_d-1)/(length(idx_d)-1)];
    for i_Tumr_Z = 1:length(set_of_Tumr_Z)
        figure(i_fig), subplot(2,2,i_Tumr_Z)
        plot(set_of_HCT_Cn(idx_HCT_Cn),gradient(squeeze(mean(mean(-log(set_of_I_s(:,:,i_beam_X,i_Tumr_Z,idx_HCT_Cn,idx_d(i_d))),1),2)),mean(diff(set_of_HCT_Cn))),'Color',TheColor,'DisplayName',['d = ',num2str(d_edges(idx_d(i_d)))],'LineWidth',2), hold on
        xlabel('HCT (%)'), xlim([min(set_of_HCT_Cn) max(set_of_HCT_Cn)])
        ylabel('OD/HCT (a.u.)'), ylim([min(Y_Limit) max(Y_Limit)])
        set(gca,'fontsize',12), grid on; legend('show','numColumns',1,'Orientation','horizontal'),
        title(['beamX = ',num2str(set_of_beam_X(i_beam_X)),' cm, Tumor depth = ',num2str(set_of_Tumr_Z(i_Tumr_Z)),' cm'])
    end
    clearvars TheColor
end
clearvars i_fig idx_d Y_Limit i_beam_X idx_HCT_Cn i_d i_Tumr_Z
end
