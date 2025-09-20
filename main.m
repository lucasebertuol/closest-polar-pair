% Algorithm to find the two closest polar points and indicate the resolution 
% that preserves all information in the Cartesian-polar mapping.
clear; clc; close all;

%N = [101, 201, 501, 1001, 1501]; % Image Size
%N = [101, 201, 501, 1001]; % Image Size
N = [201, 301, 501]; % Image Size
% [n, x1, y1, x2, y2, Deltarho, Deltatheta]
results = [];
% [n, x1, y1, x2, y2, Deltarho, Deltatheta, Rho*Theta]
results_deg = [];
% [n, dmin, resolution, Rho, Theta, Rho*Theta]
results_dmin = [];
% [n, dmin_deg, resolution, Rho, Theta, Rho*Theta]
results_dmin_deg = [];
tic
for n = N
    disp("Calculando para n = " + num2str(n) + "...");
    center = (n - 1) / 2 + 1;
    % [rho, theta, theta_deg, x, y;]
    polar_points = [];

    % Cartesian-polar
    for x = 0:center-1
        for y = 0:center-1
            rho = sqrt(x^2 + y^2);
            theta = atan2(y, x);
            % Check if theta is less than 45 degrees
            if theta <= pi/4
                polar_points = [polar_points; [rho, theta, rad2deg(theta), x, y]];
            end
        end
    end

    R = size(polar_points, 1);
    dmin = inf;
    polar_points_min = [];
    dmin_deg = inf;
    polar_points_min_deg = [];

    % Search the minimum distance
    for i=1:R-1
        for ii=i+1:R
            d = sqrt((polar_points(ii, 1) - polar_points(i, 1))^2 + ...
            (polar_points(ii, 2) - polar_points(i, 2))^2);

            if d < dmin
                dmin = d;
                polar_points_min = [polar_points(i, :); polar_points(ii, :)];
            end

            d = sqrt((polar_points(ii, 1) - polar_points(i, 1))^2 + ...
            (polar_points(ii, 3) - polar_points(i, 3))^2);

            if d < dmin_deg
                dmin_deg = d;
                polar_points_min_deg = [polar_points(i, :); polar_points(ii, :)];
            end
        end
    end
    
    resolution = dmin/sqrt(2);
    rhomax = center-1;
    Rho = ceil(rhomax / resolution + 1);
    Theta = ceil(pi/resolution);
    
    results_dmin = [results_dmin; [n, round(dmin*1000)/1000, ...
        round(resolution*1000)/1000, Rho, Theta, Rho*Theta]];

    resolution = dmin_deg/sqrt(2);
    Rho = ceil(rhomax / resolution + 1);
    Theta = ceil(360/resolution + 1);
    
    results_dmin_deg = [results_dmin_deg; [n, round(dmin_deg*1000)/1000, ...
        round(resolution*1000)/1000, Rho, Theta, Rho*Theta]];

    Deltarho = abs(polar_points_min(2, 1) - polar_points_min(1, 1));
    Deltatheta = abs(polar_points_min(2, 2) - polar_points_min(1, 2));

    Deltarho_deg = abs(polar_points_min_deg(2, 1) - polar_points_min_deg(1, 1));
    Deltatheta_deg = abs(polar_points_min_deg(2, 3) - polar_points_min_deg(1, 3));

    results = [results; [n, polar_points_min(1, 4), polar_points_min(1, 5), ...
    polar_points_min(2, 4), polar_points_min(2, 5), Deltarho, Deltatheta]];

    Rho = ceil(rhomax / Deltarho_deg + 1);
    Theta = ceil(360/Deltatheta_deg + 1);

    results_deg = [results_deg; [n, polar_points_min_deg(1, 4), ...
    polar_points_min_deg(1, 5), polar_points_min_deg(2, 4), ...
    polar_points_min_deg(2, 5), Deltarho_deg, Deltatheta_deg, Rho*Theta]];
end
toc

% disp('For mapping in rad')
% table(results(:,1), results(:,2), results(:,3), results(:,4), results(:,5), ...
% results(:,6), results(:,7), 'VariableNames', {'n', 'x1', 'y1', 'x2', 'y2', ...
% 'Deltarho', 'Deltatheta'})

% disp('For mapping in degrees')
% table(results_deg(:,1), results_deg(:,2), results_deg(:,3), results_deg(:,4), ...
% results_deg(:,5), results_deg(:,6), results_deg(:,7), results_deg(:,8), ...
% 'VariableNames', {'n', 'x1', 'y1', 'x2', 'y2', 'Deltarho', 'Deltatheta', ...
% 'Rho*Theta'})

disp('For mapping in rad')
table(results_dmin(:,1), results_dmin(:,2), results_dmin(:,3), ...
    results_dmin(:,4), results_dmin(:,5), results_dmin(:,6), 'VariableNames', ...
    {'n', 'dmin', 'resolution', 'Rho', 'Theta', 'Rho*Theta'})

% disp('For mapping in degrees')
% table(results_dmin_deg(:,1), results_dmin_deg(:,2), results_dmin_deg(:,3), ...
%     results_dmin_deg(:,4), results_dmin_deg(:,5), results_dmin_deg(:,6), ...
%     'VariableNames', {'n', 'dmin', 'resolution', 'Rho', 'Theta', 'Rho*Theta'})

% disp('For mapping in degrees')
% table(results_deg(:,1), results_deg(:,2), results_deg(:,3), results_deg(:,4), ...
% results_deg(:,5), round(results_deg(:,6)*1000)/1000, ...
% round(results_deg(:,7)*1000)/1000, results_deg(:,8), 'VariableNames', {'n', ...
% 'x1', 'y1', 'x2', 'y2', 'Deltarho', 'Deltatheta', ...
% 'Rho*Theta'})