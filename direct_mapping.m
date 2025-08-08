% Direct mapping algorithm to verify whether all information has been preserved
% for a given resolution.

clear; clc; close all;

% ---------------------------------- Settings ---------------------------------
N = 101; % Image Size
dorig = sqrt((0.7035)^2 + (0.5787)^2)
Deltarho = 1.15*(dorig/sqrt(2));
Deltatheta = Deltarho; % In degrees
% -----------------------------------------------------------------------------

center = round(N/2);
rho = sqrt((center-1)^2 + (center-1)^2);
Rho = ceil(rho / Deltarho+1); %Tudo indica que não deveria se somar 1
Theta = ceil(360/Deltatheta+1); %Tudo indica que não deveria se somar 1

% Deltarho = rho / (Rho - 1)
% Deltatheta = 360 / (Theta - 1)
dactual = sqrt(Deltarho^2 + Deltatheta^2)

% inim = double(rgb2gray(imread('eagle_512x512.png')));
% inim = inim(1:511, 1:511);
inim = ones(N, N);
% figure
% imshow(inim);

outim = zeros(Rho, Theta);
count = 0;

% Direct mapping
for r = 1:N
    for c = 1:N
        x = c - center;
        y =  center-r;

        rho = round(sqrt(x^2+y^2) / Deltarho) + 1;
        theta = round((atan2d(y, x)+180) / Deltatheta) + 1;

        if rho <= Rho && theta <= Theta
            if outim(rho, theta) ~= 0
                count = count+1;
            end
            if (rho == 86 && theta ==313)
                [r, c]
                [sqrt(x^2+y^2) / Deltarho + 1, (atan2d(y, x)+180) / Deltatheta + 1]
            end
            outim(rho, theta) = inim(r, c);
        end
    end
end
if count > 0
  disp('Não salvou todos os pontos...')
  count
end
% figure
% imshow(mat2gray(outim))

preservation_percentage = sum(outim(:)==1) / (N^2) * 100

%d para os pontos próximos a 45 graus
% rho1 = (sqrt((N - center)^2+(N - center)^2) / Deltarho)+1;
% rho2 = (sqrt((N-1 - center)^2+(N - center)^2) / Deltarho)+1;
% theta1 = ((atan2d(N-center, N - center)+ 180)/ Deltatheta)+1;
% theta2 = ((atan2d(N-1-center, N - center)+ 180)/ Deltatheta)+1;
% d = sqrt((rho1 - rho2)^2 + (theta1 - theta2)^2)
% 
% rho1 = (sqrt((2 - center)^2+(90 - center)^2) / Deltarho) + 1;
% rho2 = (sqrt((2 - center)^2+(91 - center)^2) / Deltarho) + 1;
% theta1 = ((atan2d(center-2, 90 - center) + 180)/ Deltatheta) + 1;
% theta2 = ((atan2d(center - 2, 91 - center) + 180)/ Deltatheta) + 1;
% d = sqrt((rho1 - rho2)^2 + (theta1 - theta2)^2)