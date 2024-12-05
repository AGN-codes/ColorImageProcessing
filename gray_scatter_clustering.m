%% color image kmm
close all;
clear;
clc;

%% loading the original image
Iorig = im2double(imread('goya.jpeg'));
Iorig = rgb2gray(Iorig); % uncomment if Iorig is colored
%Iorig = imresize(Iorig, 0.3); % for resizing

figure;
imshow(Iorig);
title('Original image');

%% vectorizing the original image
Iorig_vec = reshape(Iorig, [], 1);

%% histogram plot for the image
figure;
imhist(Iorig);
title('Image histogram');
[ih_counts, ih_binLocations] = imhist(Iorig);
ylim([0 max(ih_counts,[],'all')]);
% for the ylim, the following algorithm was not good
% [hc_N,hc_edges] = histcounts(Iorig, 256);
% ylim([0 max(hc_N,[],'all')]);

%% kmeans clustering
% parameters and dataset
no_cl = 3 % number of clusters
cl_data = Iorig_vec; % clustering data

% clustering the data
[idx,C,sumd,D] = kmeans(cl_data, no_cl);

% representing the data as centroids
cl_data_rep = C(idx,:);
cl_data_rep_2d = reshape(cl_data_rep, size(Iorig, 1), size(Iorig, 2));
figure;
imshow(cat(2, Iorig, cl_data_rep_2d));
title(cat(2, 'representing the data as ', num2str(no_cl), ' centroids'));

% scatter plot for centroids
figure;
imhist(cl_data_rep);
title('Image histogram');
title(cat(2, 'Image centroids histogram for ', num2str(no_cl), ' clusters'));
[ih_counts, ih_binLocations] = imhist(cl_data_rep);
ylim([0 1.05*max(ih_counts,[],'all')]);
xlim([-0.05 1.05]);



%% replacing the darkest and lightest centroids with black and white
C_BW = C;
[~,indx1] = max(C_BW);
[~,indx2] = min(C_BW);
C_BW(indx1) = 1;
C_BW(indx2) = 0;

% representing the data as centroids
cl_data_rep_BW = C_BW(idx,:);
cl_data_rep_2d_BW = reshape(cl_data_rep_BW, size(Iorig, 1), size(Iorig, 2));
figure;
imshow(cat(2, Iorig, cl_data_rep_2d, cl_data_rep_2d_BW));
title(cat(2, 'representing the data as ', num2str(no_cl), ' centroids, last one replacing max and min with BW'));

% scatter plot for centroids
figure;
imhist(cl_data_rep_BW);
title('Image histogram');
title(cat(2, 'Image centroids histogram for ', num2str(no_cl), ' clusters with BW'));
[ih_counts, ih_binLocations] = imhist(cl_data_rep_BW);
ylim([0 1.05*max(ih_counts,[],'all')]);
xlim([-0.05 1.05]);

%% stretching out the clustered picture
cl_data_rep_strch = cl_data_rep - min(cl_data_rep);
cl_data_rep_strch = cl_data_rep_strch / max(cl_data_rep_strch);
cl_data_rep_2d_strch = reshape(cl_data_rep_strch, size(Iorig, 1), size(Iorig, 2));

figure;
imshow(    cat(  1, cat(2, Iorig, cl_data_rep_2d), ...
    cat(2, cl_data_rep_2d_BW, cl_data_rep_2d_strch)  )     );
title(cat(2, 'representing the data as ', num2str(no_cl), ' centroids, streched to fill the spectrum'));

% scatter plot for centroids
figure;
imhist(cl_data_rep_2d_strch);
title('Image histogram');
title(cat(2, 'Image centroids histogram for ', num2str(no_cl), ' clusters streched to fill the spectrum'));
[ih_counts, ih_binLocations] = imhist(cl_data_rep_2d_strch);
ylim([0 1.05*max(ih_counts,[],'all')]);
xlim([-0.05 1.05]);