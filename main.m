clear; clc; close all;

N = 11; % Image Size
while N <= 100000
    center = round(N / 2);
    % [rho, theta, theta_deg;]
    polar_points = [];

    % Cartesian-polar
    for x = 0:center-1
        for y = 0:center-1
            rho = sqrt(x^2 + y^2);
            theta = atan2(y, x);
            % Check if theta is less than 45 degrees
            if theta <= pi/4
                polar_points = [polar_points; [rho, theta, rad2deg(theta)]];
            end
        end
    end

    R = size(polar_points, 1);

    for i=1:R-1
        for ii=i+1:R
        end
    end

    N = (N-1)*10+1;
end