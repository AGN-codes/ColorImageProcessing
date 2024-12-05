%% hsv color space mix
% mixes colors in an image in hsv color space
close all;
clear;
clc;

%% loading the original image
Iorig = im2double(imread('lonetree.jpg'));
% Iorig = Iorig(1:50, 1:30, :);
% Iorig = imresize(Iorig, 10);

figure;
imshow(Iorig);
title('Original image');

%% to test for single color, uncomment the following lines
% rgb_all = [255, 200, 20];
% rgb_all = rgb_all/255;
% Iorig(:,:,1) = rgb_all(1);
% Iorig(:,:,2) = rgb_all(2);
% Iorig(:,:,3) = rgb_all(3);

%% mean color
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
mean_Iorig_hsv = zeros(1,3); % place holder
mean_Iorig_hsv(2) = sqrt(mean_Iorig_mun_v(1)^2 + mean_Iorig_mun_v(2)^2); % saturation
mean_Iorig_hsv(3) = mean_Iorig_mun_v(3); % value
if mean_Iorig_mun_v(1) == 0 && mean_Iorig_mun_v(2) >= 0 % hue
    mean_Iorig_hsv(1) = (pi/2) / (2*pi);
elseif mean_Iorig_mun_v(1) == 0 && mean_Iorig_mun_v(2) < 0
    mean_Iorig_hsv(1) = (3*pi/2) / (2*pi);
elseif mean_Iorig_mun_v(1) >= 0 && mean_Iorig_mun_v(2) == 0
    mean_Iorig_hsv(1) = (0) / (2*pi);
elseif mean_Iorig_mun_v(1) < 0 && mean_Iorig_mun_v(2) == 0
    mean_Iorig_hsv(1) = (pi) / (2*pi);
elseif mean_Iorig_mun_v(1) > 0 && mean_Iorig_mun_v(2) > 0
    mean_Iorig_hsv(1) = atan(mean_Iorig_mun_v(2)/mean_Iorig_mun_v(1));
    mean_Iorig_hsv(1) = mean_Iorig_hsv(1) / (2*pi);
elseif mean_Iorig_mun_v(1) > 0 && mean_Iorig_mun_v(2) < 0
    mean_Iorig_hsv(1) = atan(mean_Iorig_mun_v(2)/mean_Iorig_mun_v(1));
    mean_Iorig_hsv(1) = (2*pi+mean_Iorig_hsv(1)) / (2*pi);
elseif mean_Iorig_mun_v(1) < 0 && mean_Iorig_mun_v(2) > 0
    mean_Iorig_hsv(1) = atan(mean_Iorig_mun_v(2)/mean_Iorig_mun_v(1));
    mean_Iorig_hsv(1) = (pi + mean_Iorig_hsv(1)) / (2*pi);
elseif mean_Iorig_mun_v(1) < 0 && mean_Iorig_mun_v(2) < 0
    mean_Iorig_hsv(1) = atan(mean_Iorig_mun_v(2)/mean_Iorig_mun_v(1));
    mean_Iorig_hsv(1) = (pi + mean_Iorig_hsv(1)) / (2*pi);
end
mean_Iorig_hsv

% convering the mix from hsv to rgb
mean_Iorig_rgb = hsv2rgb(mean_Iorig_hsv); 

Imix = Iorig;
Imix(:,:,1) = mean_Iorig_rgb(1);
Imix(:,:,2) = mean_Iorig_rgb(2);
Imix(:,:,3) = mean_Iorig_rgb(3);

figure;
imshow(cat(2, Iorig, Imix));
title('Original image vs. mix');

mean_Iorig_hsv_fulcro = mean_Iorig_hsv;
mean_Iorig_hsv_fulcro(2) = 1;
mean_Iorig_hsv_fulcro(3) = 1;
mean_Iorig_rgb_fulcro = hsv2rgb(mean_Iorig_hsv_fulcro);
Imix_fulcro = Iorig;
Imix_fulcro(:,:,1) = mean_Iorig_rgb_fulcro(1);
Imix_fulcro(:,:,2) = mean_Iorig_rgb_fulcro(2);
Imix_fulcro(:,:,3) = mean_Iorig_rgb_fulcro(3);
figure;
imshow(cat(2, Iorig, Imix, Imix_fulcro));
title('Original image vs. mix vs. mix full chroma');




%% when hue in hsv is calculated from hsv color space / hsv3d
% it can be outside [0 1]
% so i checked it here
% -> the prolem got solved (wrong inverse transform)
grad = -1:0.01:1;
hsv3D_XY = zeros(size(grad,2), size(grad,2), 2);
repmat_grad = repmat(grad, size(grad,2), 1);
hsv3D_XY(:,:,1) = repmat_grad;
hsv3D_XY(:,:,2) = repmat_grad';
hsv3D_XY_vec = reshape(hsv3D_XY, [], 2);
hsv3D_XY_vec_2hue = zeros(1, size(hsv3D_XY_vec, 2));
for i = 1:size(hsv3D_XY_vec, 1)
    this_hsv3D = hsv3D_XY_vec(i,:);
    this_h = 0; % place holder
    if this_hsv3D(1) == 0 && this_hsv3D(2) >= 0 % hue
        this_h = (pi/2) / (2*pi);
    elseif this_hsv3D(1) == 0 && this_hsv3D(2) < 0
        this_h = (3*pi/2) / (2*pi);
    elseif this_hsv3D(1) >= 0 && this_hsv3D(2) == 0
        this_h = (0) / (2*pi);
    elseif this_hsv3D(1) < 0 && this_hsv3D(2) == 0
        this_h = (pi) / (2*pi);
    elseif this_hsv3D(1) > 0 && this_hsv3D(2) > 0
        this_h = atan(this_hsv3D(2)/this_hsv3D(1));
        this_h = this_h / (2*pi);
    elseif this_hsv3D(1) > 0 && this_hsv3D(2) < 0
        this_h = atan(this_hsv3D(2)/this_hsv3D(1));
        this_h = (2*pi+this_h) / (2*pi);
    elseif this_hsv3D(1) < 0 && this_hsv3D(2) > 0
        this_h = atan(this_hsv3D(2)/this_hsv3D(1));
        this_h = (pi + this_h) / (2*pi);
    elseif this_hsv3D(1) < 0 && this_hsv3D(2) < 0
        this_h = atan(this_hsv3D(2)/this_hsv3D(1));
        this_h = (pi + this_h) / (2*pi);
    end

    hsv3D_XY_vec_2hue(i) = this_h; 
end
more1_hsv3D_XY_vec_2hue = hsv3D_XY_vec_2hue > 1;
less0_hsv3D_XY_vec_2hue = hsv3D_XY_vec_2hue < 0;
tot_more1_hsv3D_XY_vec_2hue = sum(more1_hsv3D_XY_vec_2hue, "all");
tot_less0_hsv3D_XY_vec_2hue = sum(less0_hsv3D_XY_vec_2hue, "all");

%% case
% for the problem mentioned above, which was solved
this_hsv3D = [-1 -1];
if this_hsv3D(1) == 0 && this_hsv3D(2) >= 0 % hue
    this_h = (pi/2) / (2*pi);
elseif this_hsv3D(1) == 0 && this_hsv3D(2) < 0
    this_h = (3*pi/2) / (2*pi);
elseif this_hsv3D(1) >= 0 && this_hsv3D(2) == 0
    this_h = (0) / (2*pi);
elseif this_hsv3D(1) < 0 && this_hsv3D(2) == 0
    this_h = (pi) / (2*pi);
elseif this_hsv3D(1) > 0 && this_hsv3D(2) > 0
    this_h = atan(this_hsv3D(2)/this_hsv3D(1));
    this_h = this_h / (2*pi);
elseif this_hsv3D(1) > 0 && this_hsv3D(2) < 0
    this_h = atan(this_hsv3D(2)/this_hsv3D(1));
    this_h = (2*pi+this_h) / (2*pi);
elseif this_hsv3D(1) < 0 && this_hsv3D(2) > 0
    this_h = atan(this_hsv3D(2)/this_hsv3D(1));
    this_h = (pi + this_h) / (2*pi);
elseif this_hsv3D(1) < 0 && this_hsv3D(2) < 0
    this_h = atan(this_hsv3D(2)/this_hsv3D(1));
    this_h = (pi + this_h) / (2*pi);
end
this_h;
