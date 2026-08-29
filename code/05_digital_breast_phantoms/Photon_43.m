function [] = Photon_43 ()
% Repository group: 05_digital_breast_phantoms
% Version role: primary simulation or study entry point
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: Generate a tumor to insert into the simulation phantom (code by Diego)

clc
close all
% Created on Mon Oct 28 12:55:02 2024
% Author: Diego (converted to MATLAB)

% Set lesion parameters
lesion_size = 25;    % Diameter of the lesion core
spicule_count = 100; % Number of spicules radiating from the core

% Create a 3D grid
grid_size = 100;  % Size of the cube that will contain the lesion
lesion_grid = zeros(grid_size, grid_size, grid_size);

% Define the lesion core
center = floor(grid_size / 2) + 1;
core_radius = floor(lesion_size / 2);

% Populate the core of the lesion
for x = center - core_radius : center + core_radius
    for y = center - core_radius : center + core_radius
        for z = center - core_radius : center + core_radius
            if sqrt((x - center)^2 + (y - center)^2 + (z - center)^2) <= core_radius
                lesion_grid(x, y, z) = 20;  % Higher density for the core
            end
        end
    end
end

% Set random seed for reproducibility
rng(40);

% Add spicules
for i = 1:spicule_count
    theta = rand() * 2 * pi;
    phi = rand() * pi;

    spicule_start_x = round(center + core_radius * sin(phi) * cos(theta));
    spicule_start_y = round(center + core_radius * sin(phi) * sin(theta));
    spicule_start_z = round(center + core_radius * cos(phi));

    num_points = randi([3, 19]);

    for j = 1:num_points
        x = round(spicule_start_x + j * sin(phi) * cos(theta));
        y = round(spicule_start_y + j * sin(phi) * sin(theta));
        z = round(spicule_start_z + j * cos(phi));

        if x >= 1 && x <= grid_size && y >= 1 && y <= grid_size && z >= 1 && z <= grid_size
            lesion_grid(x, y, z) = 0.5;  % Lower density for spicules
        end
    end
end

% Visualization
[x, y, z] = ind2sub(size(lesion_grid), find(lesion_grid > 0));
values = lesion_grid(lesion_grid > 0);

figure;
scatter3(x, y, z, 10, values, 'filled');
xlabel('X axis');
ylabel('Y axis');
zlabel('Z axis');
title('3D Spiculated Lesion');
axis equal;
view(3);
colormap turbo;
colorbar;
end
