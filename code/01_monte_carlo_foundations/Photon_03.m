function [] = Photon_03 ()
% Repository group: 01_monte_carlo_foundations
% Version role: primary simulation or study entry point
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% the study of Beer-Lambert law using mcxlab (MCXStudio)
% implementing Monte Carlo Simulation - MCXLab

mua = 0.2; g = 0.9; n = 1.4;
mus = 0.05:0.05:0.95;
r_w = nan(length(mus),17);
r_n = nan(length(mus),17);
for i_s = 1:length(mus)
    [r_w(i_s,:),r_n(i_s,:)] = get_measuremetns (mua, mus(i_s), g, n);
end
plot(mus,log10(r_w),'linewidth',2)
% plot(mus,10*log10(r_n),'linewidth',2)

end
function [detw,detc] = get_measuremetns (mua, mus, g, n)
% Photon Packet Transport using Monte Carlo method
clc
format long

% Setup Definition
M = zeros(101, 101, 100);   % AIR
M(:,:, 1:95) = 1;           % Material
M = uint8(M);
cfg.vol = M;
cfg.nphoton = 1e6;

% Properties: an N by 4 array, each row specifies [mua, mus, g, n] in order.
cfg.prop = [0,     0, 1.0, 1.0;   % Air
            mua, mus,   g, n ];   % Material
% time it runs: Run for a large time duration to capture all photons
cfg.tstart = 0;
cfg.tend   = 1e-2;
cfg.tstep  = 1e-2;

% Source: Pencil Source placed on top of the block at the center
cfg.srctype = 'pencil';
cfg.srcpos = [51, 51,  96];  % Placed on top
cfg.srcdir = [0,   0,  -1];  % Pointing straight downwards
cfg.isspecular = 0;          % No specular reflection
cfg.bc = 'aaaaaa';           % Boundary Conditions to make all light absorb at all the edges

% Detectors: Placed at different distances and directions on the surface of the block
cfg.detpos = [51, 51, 95, 2;
              51, 71, 95, 2;
              51, 65, 95, 2;
              51, 59, 95, 2;
              51, 53, 95, 2;
              51, 49, 95, 2;
              51, 43, 95, 2;
              51, 37, 95, 2;
              51, 31, 95, 2;
              71, 51, 95, 2;
              65, 51, 95, 2;
              59, 51, 95, 2;
              53, 51, 95, 2;
              49, 51, 95, 2;
              43, 51, 95, 2;
              37, 51, 95, 2;
              31, 51, 95, 2;
              ];

cfg.isreflect = 1;  % consider refractive index mismatch
cfg.unitinmm = 1.0; % 1 mm unit
% save trajectory
cfg.maxjumpdebug = 1;
cfg.issaveref = 1;


% Preview the model
% mcxpreview(cfg);


cfg.seed = 42;
[~, detpt, ~, ~] = mcxlab(cfg);
clc

% Count the number of photons at each detector individually
photonw = mcxdetweight(detpt,detpt.prop,detpt.unitinmm);                        % returns the detected weight of each photon
detc = accumarray(detpt.data(1,:)', ones(size(detpt.data(1,:)')))./cfg.nphoton; % sum all photons numbers captured by each detector
detw = accumarray(detpt.data(1,:)', photonw                     )./cfg.nphoton; % sum all photons weights captured by each detector


% detector_number = length(cfg.detpos);
% count = zeros(1, detector_number);
% for i=1:detector_number
%     count(i) = sum(detpt.data(1, :) == i);
% end
%
% detw    = accumarray(detpt.data(1,:)', photonw); % sum all photons captured by each detector
end
