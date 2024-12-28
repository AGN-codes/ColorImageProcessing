%% dct2 components
close all;
clear;
clc;

%% 
num_elem1 = 50; % specifying number of elements
num_elem2 = 13; % % specifying number of elements


mat_ind = zeros(num_elem1, num_elem2);
mat_ind_vec = reshape(mat_ind, 1, []);


% size of the final image
resize_pix = factorial(6);
resize_num1 = floor(resize_pix/num_elem1) * num_elem1;
resize_num2 = floor(resize_pix/num_elem2) * num_elem2;
% resize_num2 = resize_num1/num_elem1 * num_elem2; % comment for square final image


figure;
tiledlayout(num_elem2, num_elem1);

for i = 1:num_elem1*num_elem2
    mat_ind_vec(i) = 1;

    nexttile;
    imshow(imresize(rescale(idct2(reshape(mat_ind_vec, num_elem1, num_elem2)')), [resize_num1, resize_num2], "nearest"));

    mat_ind_vec(i) = 0;
end
