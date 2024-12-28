%% SVD demonstration
close all;
clear;
clc;

%% loading the original image and SVD
I = im2double(rgb2gray(imread('thomas-curryer-unsplash.jpg')));
I = imresize(I, 0.4);

figure;
imshow(I);
title('Original Image');

[U,S,V] = svd(I);

%% Most significant elements
No_SinV = 30; 

U_c = U;
U_c(:, No_SinV+1:end) = 0;
S_c = S;
S_c(No_SinV+1:end, No_SinV+1:end) = 0;
V_c = V;
V_c(:, No_SinV+1:end) = 0;

Ic = U_c *  S_c * V_c';

total_error = sum(abs(I-Ic), 'all');

figure;
% imshow(cat(2,I,Ic, rescale(I-Ic), Ic + bwmorph(edge(I), 'thicken', 2)));
imshow(cat(2,I,Ic, rescale(I-Ic)));

title('Original Image, the Reconstructed, and the difference');

%% showing the video
No_SinV_Seq = [1, 20, 50, 90, 150, 200, 250, 300];
No_SinV_Seq = [1:2:50, 55:5:150, 160:20:300, 400];

FrameRate = 60;
figure;
ax = axes;

videoFrameSequence = cat(2, No_SinV_Seq, No_SinV_Seq);
for i = videoFrameSequence
    U_c = U;
    U_c(i+1:end, i+1:end) = 0;
    S_c = S;
    S_c(i+1:end, i:end) = 0;
    V_c = V;
    V_c(i+1:end, i+1:end) = 0;
    
    Ic = U_c *  S_c * V_c';

    imshow(cat(2,I,Ic, rescale(I-Ic)), 'Parent', ax);

    %pause(1/FrameRate);
    pause(0);
    
    i
end
close;