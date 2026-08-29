function [z,t,C] = Tartrazine_01(zm,tm)
% Repository group: 10_transport_models
% Version role: primary simulation or study entry point
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Numerical solution of the tartrazine advection-diffusion PDE
%
% Inputs:
%   zm : maximum tissue depth in cm
%   tm : maximum simulation time in seconds
%
% Outputs:
%   z  : depth vector in meters
%   t  : time vector in seconds
%   C  : concentration matrix
%        rows    = time
%        columns = depth
%
% Example:
%   [z,t,C] = Tartrazine_02(6,24*3600);

clc;
close all;

%% Default inputs
if nargin < 1
    zm = 0.2;       % cm
end
if nargin < 2
    tm = 2*3600; % seconds
end

%% Physical parameters
% Effective diffusion coefficient in chicken breast
D = 6.6e-11;                % m^2/s
% Effective downward advection velocity
u = 0;                      % m/s

%% Surface concentration schedule
% Example:
% 1 mol/L from t = 0 to 30 minutes
% 2 mol/L after 30 minutes

Cs_1 = 1.0;                 % mol/L
Cs_2 = 2.0;                 % mol/L

t_switch = 30*60;           % seconds

% For constant 0.6 mol/L, use:
% Cs_1 = 0.6;
% Cs_2 = 0.6;

%% Initial concentration parameters

% Uniform initial background concentration
Ci_background = 0.00;       % mol/L

% Optional nonuniform initial concentration
initial_peak_amplitude = 0.00;  % mol/L
initial_peak_depth     = 10e-3; % m
initial_peak_width     = 3e-3;  % m

%% Numerical grids

Nz = 401;
Nt = 600;

z = linspace(0,zm*1e-2,Nz);
t = linspace(0,tm,Nt);

%% Solve the PDE numerically

% m = 0 means Cartesian one-dimensional geometry
m = 0;

sol = pdepe( ...
    m, ...
    @pde_model, ...
    @initial_condition, ...
    @boundary_conditions, ...
    z, ...
    t);

% Scalar solution
C = sol(:,:,1);

% pdepe may return fewer time points if integration stops early
t_solution = t(1:size(C,1));

%% Plot 1: numerical concentration profiles

time_min = [10 30 60 180 360 1440];

% Keep only times inside the simulation
time_min = time_min(time_min*60 <= tm);

figure('Color','w');
hold on;

for i = 1:length(time_min)

    requested_time = time_min(i)*60;

    % Interpolate numerical solution at the requested time
    C_now = interp1( ...
        t_solution, ...
        C, ...
        requested_time, ...
        'pchip');

    plot( ...
        z*1e3, ...
        C_now, ...
        'LineWidth',2, ...
        'DisplayName',sprintf('%g min',time_min(i)));
end

xlabel('Depth (mm)');
ylabel('Tartrazine concentration (mol/L)');
title('Numerical solution of tartrazine transport');

legend('Location','best');
grid on;
box on;

%% Plot 2: time-depth concentration heatmap

figure('Color','w');

contourf( ...
    t_solution/3600, ...
    z*1e3, ...
    C', ...
    60, ...
    'LineColor','k');

% Depth zero appears at the top
set(gca,'YDir','reverse');

xlabel('Time (hours)');
ylabel('Depth (mm)');
title('Numerical tartrazine concentration');

cb = colorbar;
ylabel(cb,'Concentration (mol/L)');

% White-to-orange colormap
ncol = 256;

cmap = [ ...
    ones(ncol,1), ...
    linspace(1,0.55,ncol)', ...
    linspace(1,0,ncol)' ];

colormap(cmap);

caxis([min(C(:)),max(C(:))]);

hold on;

if t_switch <= tm
    xline( ...
        t_switch/3600, ...
        'k--', ...
        'Second application', ...
        'LineWidth',1.5, ...
        'LabelVerticalAlignment','bottom');
end

box on;

%% Plot 3: surface concentration schedule

Cs_plot = arrayfun(@surface_concentration,t_solution);

figure('Color','w');

stairs(t_solution/3600,Cs_plot,'LineWidth',2);

xlabel('Time (hours)');
ylabel('Surface concentration (mol/L)');
title('Applied surface concentration');

grid on;
box on;

%% Nested PDE functions

    function [c,f,s] = pde_model(x,t_local,C_local,dCdz)
        %#ok<INUSD>

        % MATLAB pdepe form:
        %
        % c*dC/dt = d(f)/dz + s
        %
        % Here:
        %
        % f = D*dC/dz - u*C
        %
        % Therefore:
        %
        % dC/dt = D*d2C/dz2 - u*dC/dz

        c = 1;

        f = D*dCdz - u*C_local;

        s = 0;
    end

    function C0 = initial_condition(x)
        % Initial concentration C(z,0)
        %
        % The following example contains:
        %   1. a homogeneous nonzero background
        %   2. an additional Gaussian concentration region

        C0 = Ci_background + ...
             initial_peak_amplitude .* ...
             exp( ...
             -0.5 .* ...
             ((x-initial_peak_depth) ./ initial_peak_width).^2);

        % For a homogeneous nonzero initial concentration, use:
        %
        % C0 = Ci_background;

        % For zero initial concentration, use:
        %
        % C0 = 0;
    end

    function [pl,ql,pr,qr] = ...
            boundary_conditions(xl,Cl,xr,Cr,t_local)
        %#ok<INUSD>

        %% Left boundary: prescribed surface concentration

        Cs_now = surface_concentration(t_local);

        % pl + ql*f = 0
        %
        % Cl - Cs_now = 0
        pl = Cl-Cs_now;
        ql = 0;

        %% Right boundary: zero concentration gradient

        % We want:
        %
        % dC/dz = 0
        %
        % Since f = D*dC/dz-u*C, use:
        %
        % u*C + f = D*dC/dz = 0

        pr = u*Cr;
        qr = 1;
    end

    function Cs_now = surface_concentration(t_local)

        % Smooth transition scale (seconds)
        transition = 10;

        % Important times
        t1 = 10*60;   % end of first application
        t2 = 25*60;   % start of second application
        t3 = 35*60;   % end of second application

        % 0 -> 2 mol/L
        rise1 = 0.5*(1+tanh((t_local-0)/transition));

        % 2 -> 0 mol/L
        fall1 = 0.5*(1+tanh((t_local-t1)/transition));

        % 0 -> 2 mol/L
        rise2 = 0.5*(1+tanh((t_local-t2)/transition));

        % 2 -> 0 mol/L
        fall2 = 0.5*(1+tanh((t_local-t3)/transition));

        Cs_now = ...
            2*(rise1-fall1 + rise2-fall2);

    end
end
