function [] = Photon_39_9 ()
% Repository group: 04_compression_beam_and_dpf
% Version role: analysis or follow-up variant
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: s vs. d for breast tissue (compression-induced hemodynamic changes in breast tissue: 2-D measurement -> beam at several distances)
% the visualization: for lab meeting / wrong code

clc
close all
set(0,'DefaultFigureWindowStyle','docked')

rx = 1.00; ry = 1.00; rz = 1.00; % Tumor radius within the breast (in cm)(1.0 to 2.5 cm).
set_of_beam_X = [0.0,0.5,1.0,1.5]; % X deviation of beam in cm
set_of_beam_Y = [0.0,0.5,1.0,1.5]; % Y deviation of beam in cm
set_of_scal   = [0,1:0.5:3];   % set of scales
code_num = 39;
load(['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_',num2str(code_num),'.mat'])

% plot DPF
for i_s = 1:length(set_of_scal)
    for i_X = 1:length(set_of_beam_X)
        for i_Y = 1:length(set_of_beam_Y)
            beam_X = set_of_beam_X(i_X);
            beam_Y = set_of_beam_Y(i_Y);
            beam_d = sqrt(beam_X.^2+beam_Y.^2);
            scal   = set_of_scal(i_s);

            [d_avg,dpf_avg] = circumferenceAverageMatrix(squeeze(set_of_DPFx(i_s,i_X,i_Y,:)), squeeze(set_of_DPFy(i_s,i_X,i_Y,:)), squeeze(set_of_DPFs(i_s,i_X,i_Y,:,:)).',beam_X, beam_Y);
            figure(1), subplot(4,4,(i_Y-1).*length(set_of_beam_Y)+i_X),
            plot(d_avg,dpf_avg,'DisplayName',get_title(scal),'LineWidth',2), hold on
            set(gca,'fontsize',12), axis tight, grid on
            if i_Y==1, title(get_X_label(beam_X)); else, title("");end
            if i_Y==length(set_of_beam_Y), xlabel('d (cm)'); else, xlabel(''); end
            if i_X==1, ylabel([{get_Y_label(beam_Y)},{'DPF (a.u.)'}]); else, ylabel(""); end
            axis([0 +7.5 0 25])
            clearvars d_avg dpf_avg

            [d_avg,OD_avg] = circumferenceAverageMatrix(squeeze(set_of_OD_x(i_s,i_X,i_Y,:)), squeeze(set_of_OD_y(i_s,i_X,i_Y,:)), squeeze(set_of_OD_s(i_s,i_X,i_Y,:,:)).',beam_X, beam_Y);
            figure(2), subplot(4,4,(i_Y-1).*length(set_of_beam_Y)+i_X),
            plot(d_avg,-log(OD_avg),'DisplayName',get_title(scal),'LineWidth',2), hold on
            set(gca,'fontsize',12), axis tight, grid on
            if i_Y==1, title(get_X_label(beam_X)); else, title("");end
            if i_Y==length(set_of_beam_Y), xlabel('d (cm)'); else, xlabel(''); end
            if i_X==1, ylabel([{get_Y_label(beam_Y)},{'OD (a.u.)'}]); else, ylabel(""); end
            axis([0 +7.5 0 18])
            clearvars d_avg OD_avg

            clearvars beam_X beam_Y beam_d scal
        end
    end
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
function [out] = get_title(scal)
if scal == 0
    out = 'no tumor';
else
    out = ['press. ',sprintf('%.1f',scal),' times'];
end
end
function [out] = get_X_label(beam_X)
out = ['X = ',num2str(beam_X)];
end
function [out] = get_Y_label(beam_Y)
out = ['Y = ',num2str(beam_Y)];
end
function [r_avg,circ_avg] = circumferenceAverageMatrix(x_c,y_c,Xin,beam_X, beam_Y)
[X,Y] = ndgrid(x_c,y_c);
r_c = unique(sqrt((x_c-0).^2+(y_c-0).^2));
circ_avg = nan(length(r_c)-1,1);
r_avg = 1/2*(r_c(1:end-1)+r_c(2:end-0));
for i_r = 1:length(r_avg)
    mask = r_c(i_r)<=sqrt((X-beam_X).^2+(Y-beam_Y).^2)...
                    &sqrt((X-beam_X).^2+(Y-beam_Y).^2)<=r_c(i_r+1);
    circ_avg(i_r) = mean(Xin(mask));
end
end
