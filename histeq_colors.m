%% histeq on color planes
close all;
clear;
clc;

%%
Iorig = imread('lonetree.jpg');

Iblah = rgb2ycbcr(Iorig); 

I1 = Iblah;
I2 = Iblah;
I3 = Iblah;

I1(:,:,1) = histeq(Iblah(:,:,1));
I2(:,:,2) = histeq(Iblah(:,:,2));
I3(:,:,3) = histeq(Iblah(:,:,3));

I1 = ycbcr2rgb(I1);
I2 = ycbcr2rgb(I2);
I3 = ycbcr2rgb(I3);

Ir = Iorig;
Ig = Iorig;
Ib = Iorig;

Ir(:,:,1) = histeq(Iorig(:,:,1));
Ig(:,:,2) = histeq(Iorig(:,:,2));
Ib(:,:,3) = histeq(Iorig(:,:,3));

figure;
imshowpair(Iorig, I1, 'montage');
figure;
imshowpair(I2, I3, 'montage');

figure;
imshowpair(Iorig, Ir, 'montage');
figure;
imshowpair(Ig, Ib, 'montage');