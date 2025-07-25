% Direct mapping algorithm to verify whether all information has been preserved
% for a given resolution.

clear; clc; close all;

% ---------------------------------- Settings ---------------------------------
N = 201; % Image Size
Deltarho = 0.0049999;
Deltatheta = 0.57294; % In degrees
dorig = sqrt((0.0049999)^2 + (0.57294)^2)
% -----------------------------------------------------------------------------

center = round(N/2);
rho = sqrt(center^2 + center^2);
Rho = ceil(rho / Deltarho + 1);
Theta = ceil(180/Deltatheta + 1);

Deltarho = rho / (Rho - 1);
Deltatheta = 180 / (Theta - 1);
dactual = sqrt(Deltarho^2 + Deltatheta^2)

% inim = double(rgb2gray(imread('eagle_512x512.png')));
% inim = inim(1:511, 1:511);
inim = ones(N, N);
% figure
% imshow(inim);

outim = zeros(Rho, Theta);

% Direct mapping
for r = 1:center
    for c = 1:N
        x = c - center;
        y = center - r;

        rho = round(sqrt(x^2+y^2) / Deltarho) + 1;
        theta = round(atan2d(y, x) / Deltatheta) + 1;

        if rho <= Rho && theta <= Theta
            outim(rho, theta) = inim(r, c);
        end
    end
end
% figure
% imshow(mat2gray(outim))

preservation_percentage = sum(outim(:) == 1) / ((N-1)/2 * N + N) * 100