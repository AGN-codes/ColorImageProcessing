%% turning 
close all;
clear;
clc;

%%
I = im2double(imread('backside.jpeg'));
Iycbcr = rgb2ycbcr(I);
Iycbcr(:,:,[2,3]) = rand(size(I,1),size(I,2),2);
Iycbcr = ycbcr2rgb(Iycbcr);

figure;
imshow(cat(2,I,Iycbcr));

%%
J = rgb2hsv(Iycbcr);
figure;
imhist(J(:,:,1));
figure;
imhist(J(:,:,2));
figure;
imhist(J(:,:,3));