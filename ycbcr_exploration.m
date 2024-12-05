%% ycbcr exploration
close all;
clear;
clc;

%% loading the image
Iorig = imread("backside.jpeg");
Iorig = im2double(Iorig);

%% varying
I_ycbcr = Iorig; % line needed for ycbcr/rgb demonstration
I_ycbcr = rgb2ycbcr(Iorig); % ! -> comment for rgb demonstration

A = I_ycbcr(:, :, 1);
B = I_ycbcr(:, :, 2);
C = I_ycbcr(:, :, 3);

a = ones(size(Iorig, 1), size(Iorig, 2)) * mean(A, "all");
b = ones(size(Iorig, 1), size(Iorig, 2)) * mean(B, "all");
c = ones(size(Iorig, 1), size(Iorig, 2)) * mean(C, "all");

AA = I_ycbcr;
BB = I_ycbcr;
CC = I_ycbcr;

AA(:, :, 1) = a;
BB(:, :, 2) = b;
CC(:, :, 3) = c;

AA = ycbcr2rgb(AA); % ! -> comment for rgb demonstration
BB = ycbcr2rgb(BB); % ! -> comment for rgb demonstration
CC = ycbcr2rgb(CC); % ! -> comment for rgb demonstration

AAA = AA;
BBB = BB;
CCC = CC;

AAA(:, :, 1) = 0; % = [0, 0.5, 1] not that interesting result for 0
BBB(:, :, 2) = 1; % = [0, 0.5, 1] interesting result for 1
CCC(:, :, 3) = 0.5; % = [0, 0.5, 1] interesting result for 0.5

AAA = ycbcr2rgb(AAA); % ! -> comment for rgb demonstration
BBB = ycbcr2rgb(BBB); % ! -> comment for rgb demonstration
CCC = ycbcr2rgb(CCC); % ! -> comment for rgb demonstration

AAA4 = AAA; % line needed for ycbcr/rgb demonstration
BBB4 = BBB; % line needed for ycbcr/rgb demonstration
CCC4 = CCC; % line needed for ycbcr/rgb demonstration

AAA4 = rgb2ycbcr(AAA); % ! -> comment for rgb demonstration
BBB4 = rgb2ycbcr(BBB); % ! -> comment for rgb demonstration
CCC4 = rgb2ycbcr(CCC); % ! -> comment for rgb demonstration

neu = 0.07;
AAAA = (1 - neu) * AAA4 + neu * (I_ycbcr - AAA4);
BBBB = (1 - neu) * BBB4 + neu * (I_ycbcr - BBB4);
CCCC = (1 - neu) * CCC4 + neu * (I_ycbcr - CCC4);

AAAA = ycbcr2rgb(AAAA); % ! -> comment for rgb demonstration
BBBB = ycbcr2rgb(BBBB); % ! -> comment for rgb demonstration
CCCC = ycbcr2rgb(CCCC); % ! -> comment for rgb demonstration

X = I_ycbcr;
Y = I_ycbcr;
Z = I_ycbcr;

X(:, : , [1 2]) = 0;
Y(:, : , [1 3]) = 0;
Z(:, : , [2 3]) = 0;

X = ycbcr2rgb(X); % ! -> comment for rgb demonstration
Y = ycbcr2rgb(Y); % ! -> comment for rgb demonstration
Z = ycbcr2rgb(Z); % ! -> comment for rgb demonstration

%% seeing the results
figure;
imshowpair(Iorig, A, 'montage');
figure;
imshowpair(B, C, 'montage');

figure;
imshowpair(Iorig, AA, 'montage');
figure;
imshowpair(BB, CC, 'montage');

figure;
imshowpair(Iorig, AAA, 'montage');
figure;
imshowpair(BBB, CCC, 'montage');

figure;
imshowpair(Iorig, AAAA, 'montage');
figure;
imshowpair(BBBB, CCCC, 'montage');

% figure;
% imshowpair(Iorig, X, 'montage');
% figure;
% imshowpair(Y, Z, 'montage');