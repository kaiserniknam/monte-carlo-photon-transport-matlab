function [] = Photon_02 ()
% Repository group: 01_monte_carlo_foundations
% Version role: primary simulation or study entry point
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
clear, clc, close all
% implementing Monte Carlo Simulation - my cdoe (2): trying to have a tissue

Ny=400; Nx=Ny; Nz=Ny;
dx = 0.010; % mm, 10-um bin size.
dz = dx;
x = ([1:Nx]'-Nx/2)*dx;
z = [1:Nz]'*dx;
T = double(zeros(Nx,Ny,Nz)); % number of voxels in tissue cube
T = T + 4; % fill background with skin (4 = dermis)
zsurf = 1; % position of air/skin surface
for iz=1:Nz % for every depth z(iz)
    % air
    if iz<=round(zsurf/dz) % dz = voxel size
    T(:,:,iz) = 1; % 1 = air
    end
    % epidermis (100 um thick)
    if iz>round(zsurf/dz) & iz<=round((zsurf+0.100)/dz)
        T(:,:,iz) = 5; % 5 = epidermis
    end
    % blood vessel @ xc, zc, radius, oriented along y axis
    xc = 0; % [mm], center of blood vessel
    zc = Nz/2*dz; % [mm], center of blood vessel
    vesselradius = 0.5; % blood vessel radius [mm]
    for ix=1:Nx
        xd = x(ix) - xc; % vessel, x distance from vessel center
        zd = z(iz) - zc; % vessel, z distance from vessel center
        r = sqrt(xd^2 + zd^2); % r from vessel center
        if (r<=vesselradius) % if r is within vessel
            T(ix,:,iz) = 3; % 3 = blood. Extend pattern over all y.
        end
    end %ix
end % iz

imagesc(squeeze(T(:,1,:))), axis equal, axis tight


nm = 460; % requested wavelength
nmLIB = 300;
MU2 = interp1(nmLIB,MU,nm); % horizontal vector
j=5; % pointer to tissue type
tissue(j).name = 'epidermis';
B = 0; % blood volume fraction
S = 0; % saturation of hemoglobin
W = 0.75; % water volume fraction
M = 0.03; % melanosome volume fraction
F = 0; % fat volume fraction
a = 2; % mm^-1, scattering strength = reduced scattering coeff. at lambda_ref=500nm
fray = 0.0; % fraction of reduced scattering at 500nm due to Rayleigh scattering
% 1-fray = fraction of reduced scattering at 500nm due to Mie scattering
bmie = 1.0; % Mie scattering power
g1 = 0.90; % anisotropy of scatter
musp = a*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie); % reduced scattering
X = [B*S B*(1-S) W M F]'; % vertical vector
tissue(j).mua = MU2*X; % matrix multiplication --> absorption coeff.
tissue(j).mus = a/(1-g1); % scattering coeff.
tissue(j).g1 = g1; % anisotropy of scatter

end
