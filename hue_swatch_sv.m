%% Hue Swatch
% Saturation and Value Swatch with Fixed Hue
close all;
clear;
clc;

%% hsv plane with fixed hue
% place holder
hue_swatch_hsv = zeros(256, 256, 3);
% certain hue <----------------------- insert hue here
hue_swatch_hsv(:,:,1) = 100/360;
% hue_swatch_hsv(:,:,1) = 0.099091918411174;
% gradient of saturation
grad_255 = (0:255)/255;
mat_grad = repmat(grad_255, 256, 1); % rows of grad_255 stacked 256 times
hue_swatch_hsv(:,:,2) = mat_grad;
% gradient of value
grad_255 = flip(grad_255);
mat_grad = repmat(grad_255, 256, 1); % rows of grad_255 stacked 256 times
mat_grad = mat_grad'; % columns of grad_255 stacked 256 times
hue_swatch_hsv(:,:,3) = mat_grad;
% turn from hsv to rgb
hue_swatch_rgb = hsv2rgb(hue_swatch_hsv);
% plot
figure;
imshow(hue_swatch_rgb);

%% hue spectrum
% spectrum_hsv = zeros(1,256,3);
% spectrum_hsv(:,:,1) = (0:255)/255;
% spectrum_hsv(:,:,[2,3]) = 1;
% spectrum_hsv = repmat(spectrum_hsv, 256, 1, 1);
% figure;
% imshow(hsv2rgb(spectrum_hsv));