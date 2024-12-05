%% monochromatization of image with saturation factor
close all;
clear;
clc;

%% loading the original image
Iorig = im2double(imread('AGN_D.M.png'));

figure;
imshow(Iorig);
title('Original image');

%% turning into monochrome: mean color from hsv_mix.m or certain chosen hue
% converting to hsv
Iorig_hsv = rgb2hsv(Iorig);

% converting hsv image to hsv color space coordinates (cartesian coordinates)
Iorig_hsv3D = zeros(size(Iorig, 1), size(Iorig,2), 3);
Iorig_hsv3D(:,:,3) = Iorig_hsv(:,:,3);
Iorig_hsv3D(:,:,1) = Iorig_hsv(:,:,2) .* cos(Iorig_hsv(:,:,1) * 2*pi);
Iorig_hsv3D(:,:,2) = Iorig_hsv(:,:,2) .* sin(Iorig_hsv(:,:,1) * 2*pi);

% mix of image in hsv color space system (mix = mean = average)
mean_Iorig_mun_v = squeeze(mean(Iorig_hsv3D, [1 2]));

% turning the mix, a single color, from hsv color space coordinates to hsv
mean_Iorig_h =  0; % place holder
if mean_Iorig_mun_v(1) == 0 && mean_Iorig_mun_v(2) >= 0 % hue
    mean_Iorig_h = (pi/2) / (2*pi);
elseif mean_Iorig_mun_v(1) == 0 && mean_Iorig_mun_v(2) < 0
    mean_Iorig_h = - (pi/2) / (2*pi);
elseif mean_Iorig_mun_v(1) >= 0 && mean_Iorig_mun_v(2) == 0
    mean_Iorig_h = (0) / (2*pi);
elseif mean_Iorig_mun_v(1) < 0 && mean_Iorig_mun_v(2) == 0
    mean_Iorig_h = - (pi) / (2*pi);
elseif mean_Iorig_mun_v(1) > 0 && mean_Iorig_mun_v(2) > 0
    mean_Iorig_h = atan(mean_Iorig_mun_v(2)/mean_Iorig_mun_v(1));
    mean_Iorig_h = mean_Iorig_h / (2*pi);
elseif mean_Iorig_mun_v(1) > 0 && mean_Iorig_mun_v(2) < 0
    mean_Iorig_h = atan(mean_Iorig_mun_v(2)/mean_Iorig_mun_v(1));
    mean_Iorig_h = mean_Iorig_h / (2*pi);
elseif mean_Iorig_mun_v(1) < 0 && mean_Iorig_mun_v(2) > 0
    mean_Iorig_h = atan(mean_Iorig_mun_v(2)/mean_Iorig_mun_v(1));
    mean_Iorig_h = (pi + mean_Iorig_h) / (2*pi);
elseif mean_Iorig_mun_v(1) < 0 && mean_Iorig_mun_v(2) < 0
    mean_Iorig_h = atan(mean_Iorig_mun_v(2)/mean_Iorig_mun_v(1));
    mean_Iorig_h = (pi + mean_Iorig_h) / (2*pi);
end

% turning into mono chrome: changing the h of the original image hsv
Imono_hsv = Iorig_hsv;
Imono_hsv(:,:,1) = mean_Iorig_h;

% certain hue instead of mix, uncomment <------
% Imono_hsv(:,:,1) = 300/360;
% Imono_hsv(:,:,1) = 0.099091918411174;
% Imono_hsv(:,:,1) = 0.147482992563468;

% hvs2rgb
Imono_rgb = hsv2rgb(Imono_hsv);

% showing the results
figure;
imshow(cat(2, Iorig, Imono_rgb));
title('Original image vs. monochrome');

%% scaling the saturation
final_max_sat_01 = 0.5; % maximum saturation
Iorig_max_sat = max(Iorig_hsv(:,:,2), [], 'all');
Imono_hsv_satscl = Imono_hsv;
if Iorig_max_sat ~= 0
    Imono_hsv_satscl(:,:,2) = Imono_hsv(:,:,2) / Iorig_max_sat * final_max_sat_01;
end
Imono_rgb_satscl = hsv2rgb(Imono_hsv_satscl);
figure;
imshow(cat(2, Iorig, Imono_rgb_satscl));
title('Original image vs. saturation-scaled monochrome');
