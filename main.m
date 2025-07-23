% Algorithm to find the two closest polar points and indicate the resolution 
% that preserves all information in the Cartesian-polar mapping.
clear; clc; close all;

N = [11, 51, 101, 501, 1001, 5001, 10001]; % Image Size
% [n, x1, y1, x2, y2, Deltarho, Deltatheta]
results = [];
% [n, x1, y1, x2, y2, Deltarho, Deltatheta]
results_deg = [];
tic
for n = N
    disp("Calculando para n = " + num2str(n));
    center = round(n / 2);
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

    Deltarho = abs(polar_points_min(2, 1) - polar_points_min(1, 1));
    Deltatheta = abs(polar_points_min(2, 2) - polar_points_min(1, 2));

    Deltarho_deg = abs(polar_points_min_deg(2, 1) - polar_points_min_deg(1, 1));
    Deltatheta_deg = abs(polar_points_min_deg(2, 3) - polar_points_min_deg(1, 3));

    results = [results; [n, polar_points_min(1, 4), polar_points_min(1, 5), ...
    polar_points_min(2, 4), polar_points_min(2, 5), Deltarho, Deltatheta]];

    results_deg = [results_deg; [n, polar_points_min_deg(1, 4), ...
    polar_points_min_deg(1, 5), polar_points_min_deg(2, 4), ...
    polar_points_min_deg(2, 5), Deltarho_deg, Deltatheta_deg]];
end
toc

disp('For mapping in rad')
table(results(:,1), results(:,2), results(:,3), results(:,4), results(:,5), ...
results(:,6), results(:,7), 'VariableNames', {'n', 'x1', 'y1', 'x2', 'y2', ...
'Deltarho', 'Deltatheta'})

disp('For mapping in degrees')
table(results_deg(:,1), results_deg(:,2), results_deg(:,3), results_deg(:,4), ...
results_deg(:,5), results_deg(:,6), results_deg(:,7), 'VariableNames', {'n', ...
'x1', 'y1', 'x2', 'y2', 'Deltarho', 'Deltatheta'})