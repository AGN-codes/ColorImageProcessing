%% sharpening and blurring in frequency domain
close all;
clear;
clc;

%% Loading the image and seeing it
Iorig = im2double(rgb2gray(imread('thomas-curryer-unsplash-small.jpg')));
figure;
imshow(Iorig);
xlabel("original image")

%% blurring image in spatial domain
blFil = fspecial('average', 5);
Ibl = imfilter(Iorig, blFil);

figure;
imshow(cat(1, Iorig, Ibl));
xlabel('blurring result in spatial domain');

%% blurring image in frequency domain
Iorig_f = fft2(Iorig);
blFil_f = fft2(blFil, size(Iorig,1), size(Iorig,2));


Ibl_f = Iorig_f .* blFil_f;
Ibl_if = ifft2(Ibl_f);

figure;
imshow(cat(1, Iorig, Ibl_if));
xlabel('blurring result in frequency domain');

%% sharpening image in spatial domain
% implemented with laplacian
shFil_1 = fspecial('laplacian', 1);
shFil_05 = fspecial('laplacian', 0.5);
shFil_0 = fspecial('laplacian', 0);

Ish_1 = Iorig - imfilter(Ibl, shFil_1);
Ish_05 = Iorig - imfilter(Ibl, shFil_05);
Ish_0 = Iorig - imfilter(Ibl, shFil_0);

figure;
imshow(cat(2,cat(1, Iorig, Ish_1),cat(1, Ish_05, Ish_0)));
xlabel('sharpening result in spatial domain, for alpha = [orig, 1; 0.5, 0]');

%% sharpening image in frequency domain
% implemented with laplacian
Iorig_f = fft2(Iorig);
shFil_1_f = fft2(shFil_1, size(Iorig,1), size(Iorig,2));
shFil_05_f = fft2(shFil_05, size(Iorig,1), size(Iorig,2));
shFil_0_f = fft2(shFil_0, size(Iorig,1), size(Iorig,2));

Ish_1_f = Iorig_f - (Iorig_f .* shFil_1_f);
Ish_05_f = Iorig_f - (Iorig_f .* shFil_05_f);
Ish_0_f = Iorig_f - (Iorig_f .* shFil_0_f);

Ish_1_if = ifft2(Ish_1_f);
Ish_05_if = ifft2(Ish_05_f);
Ish_0_if = ifft2(Ish_0_f);

figure;
imshow(cat(2,cat(1, Iorig, Ish_1_if),cat(1, Ish_05_if, Ish_0_if)));
xlabel('sharpening result in frequency domain, for alpha = [orig, 1; 0.5, 0]');