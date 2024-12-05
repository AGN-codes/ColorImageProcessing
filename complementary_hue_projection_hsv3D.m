%% complementary-hue projection in hsv3D
close all;
clear;
clc;

%% loading the original image
Iorig = im2double(imread('lonetree.jpg'));
% Iorig = imresize(Iorig, 0.4);

figure;
imshow(Iorig);
title('Original image');

%% 
% certain hue <----------------------- insert hue here
% between 0 and 1
hue_of_proj = 0.099091918411174;
hue_of_proj = 0.147482992563468;
hue_of_proj = mod(hue_of_proj, 1);

% converting to hsv
Iorig_hsv = rgb2hsv(Iorig);

% converting hsv image to hsv color space coordinates (cartesian coordinates)
Iorig_hsv3D = zeros(size(Iorig, 1), size(Iorig,2), 3);
Iorig_hsv3D(:,:,3) = Iorig_hsv(:,:,3);
Iorig_hsv3D(:,:,1) = Iorig_hsv(:,:,2) .* cos(Iorig_hsv(:,:,1) * 2*pi);
Iorig_hsv3D(:,:,2) = Iorig_hsv(:,:,2) .* sin(Iorig_hsv(:,:,1) * 2*pi);

%
Iproj_hsv3D = Iorig_hsv3D;
Iproj_hsv3D(:,:,1) = sqrt(Iorig_hsv3D(:,:,1).^2 + Iorig_hsv3D(:,:,2).^2) .* ...
    cos(Iorig_hsv(:,:,1) * 2*pi - hue_of_proj * 2*pi) .* ...
    cos(hue_of_proj * 2*pi);
Iproj_hsv3D(:,:,2) = sqrt(Iorig_hsv3D(:,:,1).^2 + Iorig_hsv3D(:,:,2).^2) .* ...
    cos(Iorig_hsv(:,:,1) * 2*pi - hue_of_proj * 2*pi) .* ...
    sin(hue_of_proj * 2*pi);

%
Iproj_hsv3D = reshape(Iproj_hsv3D, [], 3);

%
Iproj_hsv = zeros(size(Iorig, 1), size(Iorig,2), 3);
Iproj_hsv = reshape(Iproj_hsv, [], 3);

%
i = 1;
Iproj_hsv_p = zeros(1,3); % place holder
for i = 1:size(Iproj_hsv3D,1)
    this_hsv3D = Iproj_hsv3D(i, :);
    Iproj_hsv_p(2) = sqrt(this_hsv3D(1)^2 + this_hsv3D(2)^2); % saturation
    Iproj_hsv_p(3) = this_hsv3D(3); % value
    if this_hsv3D(1) == 0 && this_hsv3D(2) >= 0 % hue
        Iproj_hsv_p(1) = (pi/2) / (2*pi);
    elseif this_hsv3D(1) == 0 && this_hsv3D(2) < 0
        Iproj_hsv_p(1) = (3*pi/2) / (2*pi);
    elseif this_hsv3D(1) >= 0 && this_hsv3D(2) == 0
        Iproj_hsv_p(1) = (0) / (2*pi);
    elseif this_hsv3D(1) < 0 && this_hsv3D(2) == 0
        Iproj_hsv_p(1) = (pi) / (2*pi);
    elseif this_hsv3D(1) > 0 && this_hsv3D(2) > 0
        Iproj_hsv_p(1) = atan(this_hsv3D(2)/this_hsv3D(1));
        Iproj_hsv_p(1) = Iproj_hsv_p(1) / (2*pi);
    elseif this_hsv3D(1) > 0 && this_hsv3D(2) < 0
        Iproj_hsv_p(1) = atan(this_hsv3D(2)/this_hsv3D(1));
        Iproj_hsv_p(1) = (2*pi+Iproj_hsv_p(1)) / (2*pi);
    elseif this_hsv3D(1) < 0 && this_hsv3D(2) > 0
        Iproj_hsv_p(1) = atan(this_hsv3D(2)/this_hsv3D(1));
        Iproj_hsv_p(1) = (pi + Iproj_hsv_p(1)) / (2*pi);
    elseif this_hsv3D(1) < 0 && this_hsv3D(2) < 0
        Iproj_hsv_p(1) = atan(this_hsv3D(2)/this_hsv3D(1));
        Iproj_hsv_p(1) = (pi + Iproj_hsv_p(1)) / (2*pi);
    end

    Iproj_hsv (i, :) = Iproj_hsv_p;
    i = i + 1;
end

%
Iproj_hsv3D = reshape(Iproj_hsv3D, size(Iorig, 1), size(Iorig, 2), 3);
Iproj_hsv = reshape(Iproj_hsv, size(Iorig, 1), size(Iorig, 2), 3);
Iproj_rgb = hsv2rgb(Iproj_hsv);

%
figure;
imshow(cat(2,Iorig,Iproj_rgb));
title('Original image vs complementary-hue projection');

%% saturation histogram before and after
figure;
imhist(Iorig_hsv(:,:,2));
[ih_counts, ih_binLocations] = imhist(Iorig_hsv(:,:,2));
ylim([0 1.05*max(ih_counts,[],'all')]);
xlim([-0.05 1.05]);

figure;
imhist(Iproj_hsv(:,:,2));
[ih_counts, ih_binLocations] = imhist(Iproj_hsv(:,:,2));
ylim([0 1.05*max(ih_counts,[],'all')]);
xlim([-0.05 1.05]);

%% hue histogram before and after
figure;
imhist(Iorig_hsv(:,:,1));
[ih_counts, ih_binLocations] = imhist(Iorig_hsv(:,:,1));
ylim([0 1.05*max(ih_counts,[],'all')]);
xlim([-0.05 1.05]);

figure;
imhist(Iproj_hsv(:,:,1));
[ih_counts, ih_binLocations] = imhist(Iproj_hsv(:,:,1));
ylim([0 1.05*max(ih_counts,[],'all')]);
xlim([-0.05 1.05]);

%% for single hue projection
%
same_hue = (Iproj_hsv(:,:,1)  == hue_of_proj);
figure;
imhist(same_hue);
[ih_counts, ih_binLocations] = imhist(same_hue);
ylim([0 1.05*max(ih_counts,[],'all')]);
xlim([-0.05 1.05]);

%
diff_hue = Iproj_hsv(:,:,1) - hue_of_proj;
figure;
imhist(diff_hue);
[ih_counts, ih_binLocations] = imhist(diff_hue);
ylim([0 1.05*max(ih_counts,[],'all')]);
xlim([-0.05 1.05]);

%% an image of different hues
hue_spec_inhsv = zeros(1,256,3);
hue_spec_inhsv(:,:,1) = (0:255)/255;
hue_spec_inhsv(:,:,[2,3]) = 1;
hue_spec_inhsv = repmat(hue_spec_inhsv, 256, 1, 1);
hue_spec_inhsv = cat(2, hue_spec_inhsv, hue_spec_inhsv);
hue_spec_inrgb = hsv2rgb(hue_spec_inhsv);
Iorig = hue_spec_inrgb;