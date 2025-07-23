% Direct mapping algorithm to verify whether all information has been preserved
% for a given resolution.

clear; clc; close all;

% ---------------------------------- Settings ---------------------------------
N = 101; % Image Size
Deltarho = ;
Deltatheta = ; % In degrees
% -----------------------------------------------------------------------------

center = round(N/2);
rho = sqrt(center^2 + center^2);
Rho = ceil(rho / Deltarho);
Theta = ceil(180/Deltatheta + 1);

inim = ones(N, N);
figure
imshow(inim);

outim = zeros(Rho, Theta);

% Direct mapping
for r = 1:center
    for c = 1:N
        x = c - center;
        y = center - r;

        rho = round(sqrt(x^2+y^2) / Deltarho) + 1;
        theta = round(rad2deg(atan2(y, x)) / Deltatheta) + 1;

        if rho <= Rho && theta <= Theta
            outim(rho, theta) = inim(r, c);
        end
    end
end
figure
imshow(outim)

preservation_percentage = sum(outim(:) == 1) / (N * ceil(N / 2)) * 100