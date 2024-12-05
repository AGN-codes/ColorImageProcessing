%% 
close all;
clear;
clc;

%% 
I = imread('goya.jpeg');
II = rgb2ycbcr(I);
II(:,:,[2 3]) = 255 - II(:,:,[2 3]);
II = ycbcr2rgb(II);
figure;
imshow(cat(2,I,II,0.6*I+0.4*II));

I = im2double(I);
II =  rgb2hsv(I);
figure;
imshow(cat(2,I, II, 0.6*I+0.4*II));

I2 = II;
I2(:,:,1) = ones(size(I,1),size(I,2))*mean(I2(:,:,1), "all");
I3 = II;
I3(:,:,2) = ones(size(I,1),size(I,2))*mean(I3(:,:,2), "all");
I4 = II;
I4(:,:,3) = ones(size(I,1),size(I,2))*mean(I4(:,:,3), "all");

I2 = hsv2rgb(I2);
I3 = hsv2rgb(I3);
I4 = hsv2rgb(I4);

figure;
imshow(cat(1,cat(2,I,I2),cat(2,I3,I4)));

II2 = II;
II2(:,:,1) = histeq(II2(:,:,1));
II3 = II;
II3(:,:,2) = histeq(II3(:,:,2));
II4 = II;
II4(:,:,3) = histeq(II4(:,:,3));

II2 = hsv2rgb(II2);
II3 = hsv2rgb(II3);
II4 = hsv2rgb(II4);

figure;
imshow(cat(1,cat(2,I,II2),cat(2,II3,II4)));

figure;
x = rgb2ycbcr(I4);
x = I4(:,:,1);
imshow(cat(2,x,histeq(x)));

y = II(:,:,3);
yeq = histeq(y);
z = rgb2ycbcr(I);
z = z(:,:,1);
zeq = histeq(z);
figure;
imshow(cat(1,cat(2,y,yeq),cat(2,z,zeq)));


% figure;
% imhist(II(:,:,1));
% figure;
% imhist(II(:,:,2));
% figure;
% imhist(II(:,:,3));