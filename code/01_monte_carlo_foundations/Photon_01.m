function [] = Photon_01 ()
% Repository group: 01_monte_carlo_foundations
% Version role: primary simulation or study entry point
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Photon Packet Transport using Monte Carlo method
% implementing Monte Carlo Simulation - my cdoe (1) + ChatGPT

clc
close all
set(0,'DefaultFigureWindowStyle','docked')

% Simulation parameters
nPhotons = 1000; % Number of photon packets
nSteps = 100;     % Maximum number of steps per photon packet
mu_a = 0.01;      % Absorption coefficient
mu_s = 0.1;       % Scattering coefficient
g = 0.9;          % Anisotropy factor for Henyey-Greenstein phase function

% Initialize arrays to store results
absorption = zeros(1, nPhotons); % Record absorption events
the_ppaths = nan(nPhotons,nSteps,3); % Record absorption events
final_locs = nan(3, nPhotons); % Record absorption events

% Loop over each photon packet
for photon = 1:nPhotons
    % Initialize photon position and direction
    position = [0, 0, 0];  % Start at the origin
    direction = [0, 0, 1]; % Initially moving along z-axis

    for step = 1:nSteps
        % Compute step size based on scattering and absorption coefficients
        s = -log(rand) / (mu_a + mu_s);

        % Update photon position
        the_ppaths(photon,step,:) = position;
        position = position + s * direction;

        % Check for absorption
        if rand < (mu_a / (mu_a + mu_s))
            absorption(photon) = 1; % Photon absorbed
            final_locs(:,photon) = position; % final position
            break;
        end

        % Update photon direction using the Henyey-Greenstein phase function
        cos_theta = (1 / (2*g)) * (1 + g^2 - ((1 - g^2) / (1 - g + 2*g*rand))^2);
        sin_theta = sqrt(1 - cos_theta^2);
        phi = 2 * pi * rand;

        % Update direction cosines
        new_direction = [sin_theta * cos(phi), sin_theta * sin(phi), cos_theta];

        % Rotate new direction to the current photon's frame
        direction = new_direction;
    end
end

% Calculate percentage of absorbed photons
percentAbsorbed = sum(absorption) / nPhotons * 100;
fprintf('Percentage of photons absorbed: %.2f%%\n', percentAbsorbed);

% Plotting absorbed photon positions (optional, for visualization)
figure(2),
for photon = 1:nPhotons
    plot3(-the_ppaths(photon,:,1),-the_ppaths(photon,:,2),-the_ppaths(photon,:,3),'linewidth',1), hold on;
end
title('Photon paths');
xlabel('x');
ylabel('y');
zlabel('z');
grid on; view([90 0]), axis tight

figure(3), scatter3(-final_locs(1,:),-final_locs(2,:),-final_locs(3,:),'filled','color','k'), hold on
title('Photon abs. positions');
xlabel('x');
ylabel('y');
zlabel('z');
grid on; view([90 0]), axis tight

end
