%% dct2 components
close all;
clear;
clc;

%% 
num_elem = 5; % specifying number of elements

mat_ind = zeros(num_elem);
mat_ind_vec = reshape(mat_ind, 1, []);

resize_pix = factorial(6); % size of the final image
resize_scale = floor(resize_pix/num_elem);

figure;
tiledlayout(num_elem, num_elem);

for i = 1:num_elem^2
    mat_ind_vec(i) = 1;

    nexttile;
    imshow(imresize(rescale(idct2(reshape(mat_ind_vec, num_elem, num_elem)')), resize_scale, "nearest"));

    mat_ind_vec(i) = 0;
end

%% index and reshape formulas
% num_elem = 3;
% mat_ind = zeros(num_elem);
% mat_ind_vec = reshape(mat_ind, 1, []);
% mat_ind_vec(num_elem+1) = 1;
% mat_ind = reshape(mat_ind_vec, num_elem, num_elem);