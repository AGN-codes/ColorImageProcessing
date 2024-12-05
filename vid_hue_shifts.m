%% video of hue shift
% run each section manually
close all;
clear;
clc;

%% loading the original image
Iorig = im2double(imread('goya.jpeg'));

figure;
imshow(Iorig);
title('Original image');

%% creating the 4D image matrix
hue_num_steps = 100;
Vrgb = zeros(hue_num_steps, size(Iorig, 1), size(Iorig, 2), 3);

Ihsv = rgb2hsv(Iorig);
for i = 1:hue_num_steps
    Vrgb(i,:,:,:) = Ihsv;
    Vrgb(i,:,:,1) = mod(Vrgb(i,:,:,1) + (i-1)/hue_num_steps, 1);
    Vrgb(i,:,:,:) = hsv2rgb(squeeze(Vrgb(i,:,:,:)));
    i
end

%% showing the shifts
FrameRate = 60;
figure;
ax = axes;

vidFrameSeq = 1:hue_num_steps;
vidFrameSeq = cat(2, vidFrameSeq, vidFrameSeq);
for i = vidFrameSeq
    imshow(squeeze(Vrgb(i,:,:,:)), 'Parent', ax);
    pause(1/FrameRate);

end
close;