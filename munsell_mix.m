%% munsell color space mix - needs more work
% mixes colors in an image
% THIS WAS A FAILURE
% REFER TO hsv_mix.m FOR MORE INFO
close all;
clear;
clc;

%% loading the original image
Iorig = im2double(imread('goya.jpeg'));
% Iorig = Iorig(1:50, 1:30, :);
% Iorig = imresize(Iorig, 10);

figure;
imshow(Iorig);
title('Original image');

%% to test for single color, uncomment the following lines
% rgb_all = [255, 0, 0];
% Iorig(:,:,1) = rgb_all(1)/255;
% Iorig(:,:,2) = rgb_all(2)/255;
% Iorig(:,:,3) = rgb_all(3)/255;

%% mean color
% converting to hsv & ycbcr
Iorig_hsv = rgb2hsv(Iorig);
Iorig_ycbcr = rgb2ycbcr(Iorig);

% selecting the hue applied to the whole image
Iorig_hcyv = cat(3, Iorig_hsv(:,:,1), Iorig_hsv(:,:,2).*Iorig_hsv(:,:,3), Iorig_ycbcr(:,:,1), Iorig_hsv(:,:,3));
Iorig_mun_v = zeros(size(Iorig, 1), size(Iorig,2), 4);
Iorig_mun_v(:,:,3) = Iorig_hcyv(:,:,3);
Iorig_mun_v(:,:,4) = Iorig_hcyv(:,:,4);
Iorig_mun_v(:,:,1) = Iorig_hcyv(:,:,2) .* cos(Iorig_hcyv(:,:,1) * 2*pi);
Iorig_mun_v(:,:,2) = Iorig_hcyv(:,:,2) .* sin(Iorig_hcyv(:,:,1) * 2*pi);
mean_Iorig_mun_v = squeeze(mean(Iorig_mun_v, [1 2]));
mean_Iorig_hsv =  zeros(1,3);
mean_Iorig_chroma = sqrt(mean_Iorig_mun_v(1)^2 + mean_Iorig_mun_v(2)^2);
mean_Iorig_hsv(2) = mean_Iorig_chroma * mean_Iorig_mun_v(4);
mean_Iorig_hsv(3) = mean_Iorig_mun_v(4);
if mean_Iorig_mun_v(1) == 0 && mean_Iorig_mun_v(2) >= 0
    mean_Iorig_hsv(1) = (pi/2) / (2*pi);
elseif mean_Iorig_mun_v(1) == 0 && mean_Iorig_mun_v(2) < 0
    mean_Iorig_hsv(1) = - (pi/2) / (2*pi);
elseif mean_Iorig_mun_v(1) >= 0 && mean_Iorig_mun_v(2) == 0
    mean_Iorig_hsv(1) = (0) / (2*pi);
elseif mean_Iorig_mun_v(1) < 0 && mean_Iorig_mun_v(2) == 0
    mean_Iorig_hsv(1) = - (pi) / (2*pi);
elseif mean_Iorig_mun_v(1) > 0 && mean_Iorig_mun_v(2) > 0
    mean_Iorig_hsv(1) = atan(mean_Iorig_mun_v(2)/mean_Iorig_mun_v(1));
    mean_Iorig_hsv(1) = mean_Iorig_hsv(1) / (2*pi);
elseif mean_Iorig_mun_v(1) > 0 && mean_Iorig_mun_v(2) < 0
    mean_Iorig_hsv(1) = atan(mean_Iorig_mun_v(2)/mean_Iorig_mun_v(1));
    mean_Iorig_hsv(1) = mean_Iorig_hsv(1) / (2*pi);
elseif mean_Iorig_mun_v(1) < 0 && mean_Iorig_mun_v(2) > 0
    mean_Iorig_hsv(1) = atan(mean_Iorig_mun_v(2)/mean_Iorig_mun_v(1));
    mean_Iorig_hsv(1) = (pi + mean_Iorig_hsv(1)) / (2*pi);
elseif mean_Iorig_mun_v(1) < 0 && mean_Iorig_mun_v(2) < 0
    mean_Iorig_hsv(1) = atan(mean_Iorig_mun_v(2)/mean_Iorig_mun_v(1));
    mean_Iorig_hsv(1) = (pi + mean_Iorig_hsv(1)) / (2*pi);
end
mean_Iorig_rgb = hsv2rgb(mean_Iorig_hsv);

Imono = Iorig;
Imono(:,:,1) = mean_Iorig_rgb(1);
Imono(:,:,2) = mean_Iorig_rgb(2);
Imono(:,:,3) = mean_Iorig_rgb(3);

figure;
imshow(cat(2, Iorig, Imono));
title('Original image vs. mixed image');

%% difference between gray levels in different color spaces
Iorig_hsv = rgb2hsv(Iorig);
Iorig_ycbcr = rgb2ycbcr(Iorig);
Iorig_gray = rgb2gray(Iorig);
Iorig_lightness = rgb2lightness(Iorig)/100;

figure;
imshow(    cat(  1, cat(2, Iorig_gray, Iorig_ycbcr(:,:,1)), ...
    cat(2, Iorig_hsv(:,:,3), Iorig_lightness)  )    );

figure;
tiledlayout(2,2);
nexttile;
imhist(Iorig_gray);
[ih_counts, ~] = imhist(Iorig_gray);
ylim([0 max(ih_counts,[],'all')]);

nexttile;
imhist(Iorig_ycbcr(:,:,1));
[ih_counts, ~] = imhist(Iorig_ycbcr(:,:,1));
ylim([0 max(ih_counts,[],'all')]);

nexttile;
imhist(Iorig_hsv(:,:,3));
[ih_counts, ~] = imhist(Iorig_hsv(:,:,3));
ylim([0 max(ih_counts,[],'all')]);

nexttile;
imhist(Iorig_lightness);
[ih_counts, ~] = imhist(Iorig_lightness);
ylim([0 max(ih_counts,[],'all')]);

%% replacing hsv's value with true value
% unfinished
hue_spec_inhsv = zeros(1,256,3);
hue_spec_inhsv(:,:,1) = (0:255)/255;
hue_spec_inhsv(:,:,[2,3]) = 1;
hue_spec_inhsv = repmat(hue_spec_inhsv, 256, 1, 1);
hue_spec_inhsv = cat(2, hue_spec_inhsv, hue_spec_inhsv);
hue_spec_inrgb = hsv2rgb(hue_spec_inhsv);
figure;
imshow(hue_spec_inrgb);
hue_spec_ntsc = rgb2ntsc(hue_spec_inrgb);
hue_spec_ingray = hue_spec_ntsc(:,:,1);
hue_spec_ingray_inrgb = ntsc2rgb(    cat(  3,hue_spec_ingray, ...
    zeros(size(hue_spec_inhsv,1), size(hue_spec_inhsv,2)), ...
    zeros(size(hue_spec_inhsv,1), size(hue_spec_inhsv,2))  )    );
figure;
imshow(cat(1, hue_spec_ingray, hue_spec_inhsv(:,:,3)));
figure;
imshow(cat(1, hue_spec_inrgb, hue_spec_ingray_inrgb));

hue_spec_hs_inhsv = hue_spec_inhsv;
hue_spec_hs_inhsv(:,:,2) = 0.5;
hue_spec_hv_inhsv = hue_spec_inhsv;
hue_spec_hs_inhsv(:,:,3) = 0.5;


%%
A = hue_spec_inrgb - hue_spec_ingray_inrgb;
B = squeeze(A(1, :, :));
figure;
hold on;
plot(B(:,1));
plot(B(:,2));
plot(B(:,3));
hold off;

