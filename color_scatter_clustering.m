%% color image kmm
close all;
clear;
clc;

%% loading the original image
Iorig = im2double(imread('AGN_D.M.png'));

figure;
imshow(Iorig);
title('Original image');

%% vectorizing the original image
% Iorig = zeros(3,3,3);
% Iorig(:,:,1) = reshape(1:9, [3,3]);
% Iorig(:,:,2) = 2*reshape(1:9, [3,3]);
% Iorig(:,:,3) = 3*reshape(1:9, [3,3]);
Iorig_vec = reshape(Iorig, [], 3);

%% RGB scatter3 for the image
% resizing
sca_no = 5e3; % number of data points in the scatter3 plot
if sca_no > (size(Iorig, 1)*size(Iorig, 2))
    sca_no = size(Iorig, 1)*size(Iorig, 2);
end
scale_fac = sqrt(sca_no ./ (size(Iorig, 1)*size(Iorig, 2)));
Iresize = imresize(Iorig, scale_fac);

figure;
imshow(Iresize);
title(cat(2, 'Resized image with scale factor of ', num2str(scale_fac)));

%vectorizing the resized image
Iresize_vec = reshape(Iresize, [], 3);

% plotting
figure;
scatter3(Iresize_vec(:,1),Iresize_vec(:,2),Iresize_vec(:,3),'.');
title('Resized image RGB scatter plot');
xlabel('Red');
ylabel('Green');
zlabel('Blue');
xlim([0 1]);
ylim([0 1]);
zlim([0 1]);

%% HSV scatter3 the image
% turning rgb to hsv
Iorig_hsv = rgb2hsv(Iorig);

% vectorizing the image
Iorig_hsv_vec = reshape(Iorig_hsv, [], 3);

% resizing
sca_no = 5e4; % number of data points in the scatter3 plot
if sca_no > (size(Iorig, 1)*size(Iorig, 2))
    sca_no = size(Iorig, 1)*size(Iorig, 2);
end
scale_fac = sqrt(sca_no ./ (size(Iorig, 1)*size(Iorig, 2)));
Iresize_hsv = imresize(Iorig_hsv, scale_fac);

figure;
imshow(hsv2rgb(Iresize_hsv));
title(cat(2, 'Resized image with scale factor of ', num2str(scale_fac)));

%vectorizing the resized image
Iresize_hsv_vec = reshape(Iresize_hsv, [], 3);

% plotting
figure;
scatter3(Iresize_hsv_vec(:,1),Iresize_hsv_vec(:,2),Iresize_hsv_vec(:,3),'.');
title('Resized image HSV scatter plot');
xlabel('Hue');
ylabel('Saturation');
zlabel('Value');
xlim([0 1]);
ylim([0 1]);
zlim([0 1]);

%% HSV / HSV3D color space scatter3 the image
% turning rgb to hsv
Iorig_hsv = rgb2hsv(Iorig);

% converting hsv image to hsv color space coordinates (cartesian coordinates)
Iorig_hsv3D = zeros(size(Iorig, 1), size(Iorig,2), 3);
Iorig_hsv3D(:,:,3) = Iorig_hsv(:,:,3);
Iorig_hsv3D(:,:,1) = Iorig_hsv(:,:,2) .* cos(Iorig_hsv(:,:,1) * 2*pi);
Iorig_hsv3D(:,:,2) = Iorig_hsv(:,:,2) .* sin(Iorig_hsv(:,:,1) * 2*pi);

% vectorizing the image
Iorig_hsv3D_vec = reshape(Iorig_hsv3D, [], 3);

% resizing
sca_no = 5e4; % number of data points in the scatter3 plot
if sca_no > (size(Iorig, 1)*size(Iorig, 2))
    sca_no = size(Iorig, 1)*size(Iorig, 2);
end
scale_fac = sqrt(sca_no ./ (size(Iorig, 1)*size(Iorig, 2)));
Iresize_hsv3D_inRGB = imresize(Iorig, scale_fac);
Iresize_hsv3D = imresize(Iorig_hsv3D, scale_fac);

figure;
imshow(Iresize_hsv3D_inRGB);
title(cat(2, 'Resized image with scale factor of ', num2str(scale_fac)));

%vectorizing the resized image
Iresize_hsv3D_vec = reshape(Iresize_hsv3D, [], 3);

% plotting
figure;
scatter3(Iresize_hsv3D_vec(:,1),Iresize_hsv3D_vec(:,2),Iresize_hsv3D_vec(:,3),'.');
title('Resized image HSV3D / HSV color space scatter plot');
xlabel('Saturation Cos(Hue)');
ylabel('Saturation Sin(Hue)');
zlabel('Value');
xlim([-1 1]);
ylim([-1 1]);
zlim([0 1]);

%% YCbCr scatter3 the image
% turning rgb to hsv
Iorig_ycbcr = rgb2ycbcr(Iorig);

% vectorizing the image
Iorig_ycbcr_vec = reshape(Iorig_ycbcr, [], 3);

% resizing
sca_no = 50000; % number of data points in the scatter3 plot
if sca_no > (size(Iorig, 1)*size(Iorig, 2))
    sca_no = size(Iorig, 1)*size(Iorig, 2);
end
scale_fac = sqrt(sca_no ./ (size(Iorig, 1)*size(Iorig, 2)));
Iresize_ycbcr = imresize(Iorig_ycbcr, scale_fac);

figure;
imshow(ycbcr2rgb(Iresize_ycbcr));
title(cat(2, 'Resized image with scale factor of ', num2str(scale_fac)));

%vectorizing the resized image
Iresize_ycbcr_vec = reshape(Iresize_ycbcr, [], 3);

% plotting
figure;
scatter3(Iresize_ycbcr_vec(:,1),Iresize_ycbcr_vec(:,2),Iresize_ycbcr_vec(:,3),'.');
title('Resized image YCbCr scatter plot');
xlabel('Y');
ylabel('Cb');
zlabel('Cr');
xlim([0 1]);
ylim([0 1]);
zlim([0 1]);

%% kmeans clustering with RGB
% parameters and dataset
no_cl = 9; % number of clusters
cl_data = Iorig_vec; % clustering data

% clustering the data
[idx,C,sumd,D] = kmeans(cl_data, no_cl);

% representing the data as centroids
cl_data_rep = C(idx,:);
cl_data_rep_2d = reshape(cl_data_rep, size(Iorig, 1), size(Iorig, 2), 3);
figure;
imshow(cat(2, Iorig, cl_data_rep_2d));
title(cat(2, 'representing the data as ', num2str(no_cl), ' RGB centroids'));

% scatter plot for centroids
figure;
scatter3(C(:,1),C(:,2),C(:,3),'.');
title(cat(2, 'RGB centroids scatter plot for ', num2str(no_cl), ' clusters'));
xlabel('Red');
ylabel('Green');
zlabel('Blue');
xlim([0 1]);
ylim([0 1]);
zlim([0 1]);

%% kmeans clustering with HSV
% parameters and dataset
no_cl = 9; % number of clusters
cl_data = Iorig_hsv_vec; % clustering data

% clustering the data
[idx,C,sumd,D] = kmeans(cl_data, no_cl);

% representing the data as centroids
cl_data_rep = C(idx,:);
cl_data_rep_2d = reshape(cl_data_rep, size(Iorig, 1), size(Iorig, 2), 3);
figure;
imshow(cat(2, Iorig, hsv2rgb(cl_data_rep_2d)));
title(cat(2, 'representing the data as ', num2str(no_cl), ' HSV centroids'));

% scatter plot for centroids
figure;
scatter3(C(:,1),C(:,2),C(:,3),'.');
title(cat(2, 'HSV centroids scatter plot for ', num2str(no_cl), ' clusters'));
xlabel('Hue');
ylabel('Saturation');
zlabel('Value');
xlim([0 1]);
ylim([0 1]);
zlim([0 1]);

%% kmeans clustering with HSV color space / HSV3D
% parameters and dataset
no_cl = 12; % number of clusters
cl_data = Iorig_hsv3D_vec; % clustering data

% clustering the data
[idx,C,sumd,D] = kmeans(cl_data, no_cl);

% centroids from HSV color space / HSV3D to rgb
% from hsv_mix.m
C_rgb = C;
for i = 1:no_cl
    this_C_hsv3D = C(i,:);
    this_C_in_hsv = zeros(1,3); % place holder
    this_C_in_hsv(2) = sqrt(this_C_hsv3D(1)^2 + this_C_hsv3D(2)^2); % saturation
    this_C_in_hsv(3) = this_C_hsv3D(3); % value
    if this_C_hsv3D(1) == 0 && this_C_hsv3D(2) >= 0 % hue
        this_C_in_hsv(1) = (pi/2) / (2*pi);
    elseif this_C_hsv3D(1) == 0 && this_C_hsv3D(2) < 0
        this_C_in_hsv(1) = (3 * pi/2) / (2*pi);
    elseif this_C_hsv3D(1) >= 0 && this_C_hsv3D(2) == 0
        this_C_in_hsv(1) = (0) / (2*pi);
    elseif this_C_hsv3D(1) < 0 && this_C_hsv3D(2) == 0
        this_C_in_hsv(1) = (pi) / (2*pi);
    elseif this_C_hsv3D(1) > 0 && this_C_hsv3D(2) > 0
        this_C_in_hsv(1) = atan(this_C_hsv3D(2)/this_C_hsv3D(1));
        this_C_in_hsv(1) = this_C_in_hsv(1) / (2*pi);
    elseif this_C_hsv3D(1) > 0 && this_C_hsv3D(2) < 0
        this_C_in_hsv(1) = atan(this_C_hsv3D(2)/this_C_hsv3D(1));
        this_C_in_hsv(1) = (2*pi+this_C_in_hsv(1)) / (2*pi);
    elseif this_C_hsv3D(1) < 0 && this_C_hsv3D(2) > 0
        this_C_in_hsv(1) = atan(this_C_hsv3D(2)/this_C_hsv3D(1));
        this_C_in_hsv(1) = (pi + this_C_in_hsv(1)) / (2*pi);
    elseif this_C_hsv3D(1) < 0 && this_C_hsv3D(2) < 0
        this_C_in_hsv(1) = atan(this_C_hsv3D(2)/this_C_hsv3D(1));
        this_C_in_hsv(1) = (pi + this_C_in_hsv(1)) / (2*pi);
    end

    C_rgb(i,:) = hsv2rgb(this_C_in_hsv); 
end

% representing the data as centroids
cl_data_rep = C_rgb(idx,:);
cl_data_rep_2d = reshape(cl_data_rep, size(Iorig, 1), size(Iorig, 2), 3);
figure;
imshow(cat(2, Iorig, cl_data_rep_2d));
title(cat(2, 'representing the data as ', num2str(no_cl), ' HSV3D centroids'));

% scatter plot for centroids
figure;
scatter3(C(:,1),C(:,2),C(:,3),'.');
title(cat(2, 'HSV3D centroids scatter plot for ', num2str(no_cl), ' clusters'));
xlabel('Saturation Cos(Hue)');
ylabel('Saturation Sin(Hue)');
zlabel('Value');
xlim([-1 1]);
ylim([-1 1]);
zlim([0 1]);

% scatter plot for centroids
figure;
scatter3(C_rgb(:,1),C_rgb(:,2),C_rgb(:,3),'.');
title(cat(2, 'HSV3D centroids scatter plot for ', num2str(no_cl), ' clusters in RGB'));
xlabel('Red');
ylabel('Green');
zlabel('Blue');
xlim([0 1]);
ylim([0 1]);
zlim([0 1]);

%% kmeans clustering with YCbCr
% parameters and dataset
no_cl = 9; % number of clusters
cl_data = Iorig_ycbcr_vec; % clustering data

% clustering the data
[idx,C,sumd,D] = kmeans(cl_data, no_cl);

% representing the data as centroids
cl_data_rep = C(idx,:);
cl_data_rep_2d = reshape(cl_data_rep, size(Iorig, 1), size(Iorig, 2), 3);
figure;
imshow(cat(2, Iorig, ycbcr2rgb(cl_data_rep_2d)));
title(cat(2, 'representing the data as ', num2str(no_cl), ' YCbCr centroids'));

% scatter plot for centroids
figure;
scatter3(C(:,1),C(:,2),C(:,3),'.');
title(cat(2, 'YCbCr centroids scatter plot for ', num2str(no_cl), ' clusters'));
xlabel('Y');
ylabel('Cb');
zlabel('Cr');
xlim([0 1]);
ylim([0 1]);
zlim([0 1]);

%% kmeans clustering with RGB/HSV/YCbCr
%%% sca_no should be equal in scatter plotting sections for RGB/HSV/YCbCr
% parameters and dataset
no_cl = 20; % number of clusters
cl_data = cat(2, Iorig_vec, Iorig_hsv_vec, Iorig_ycbcr_vec); % clustering data

% clustering the data
[idx,C,sumd,D] = kmeans(cl_data, no_cl);

% representing the data as centroids
cl_data_rep = C(idx,:);
cl_data_rep_2d = reshape(cl_data_rep, size(Iorig, 1), size(Iorig, 2), []);
cl_data_rep_2d_rgb = cl_data_rep_2d(:,:,1:3);
cl_data_rep_2d_hsv = hsv2rgb(cl_data_rep_2d(:,:,4:6));
cl_data_rep_2d_ycbcr = ycbcr2rgb(cl_data_rep_2d(:,:,7:9));

% figure;
% imshow(  cat( 1, cat(2, Iorig, cl_data_rep_2d_rgb), ...
%     cat(2, cl_data_rep_2d_hsv, cl_data_rep_2d_ycbcr))  );
% title(cat(2, 'representing the data as ', num2str(no_cl), ' RGB/HSV/YCbCr centroids'));

% figure;
% imshow(cat(2, Iorig, (cl_data_rep_2d_rgb+cl_data_rep_2d_hsv+cl_data_rep_2d_ycbcr)/3));
% title(cat(2, 'representing the data as average', num2str(no_cl), ' RGB/HSV/YCbCr centroids'));

figure;
imshow(cat(2, Iorig, cl_data_rep_2d_rgb));
title(cat(2, 'representing the data as RGB', num2str(no_cl), ' RGB/HSV/YCbCr centroids'));
