%% signle-hue projection in hsv3D
close all;
clear;
clc;

%% loading the original image
Iorig = im2double(imread('goya.jpeg'));

figure;
imshow(Iorig);
title('Original image');

%% initials
% hue shift starts from this <----------------------- insert hue here
% between 0 and 1
hue_of_proj = HSV3D_mix(Iorig);

% the size of the steps
hue_num_steps = 100;

%% showing the video
FrameRate = 60;
figure;
ax = axes;

videoFrameSequence = 1:hue_num_steps;
videoFrameSequence = cat(2, videoFrameSequence, videoFrameSequence);
for i = videoFrameSequence
    imshow(cat(2,Iorig,monochromeHSV3D(Iorig, hue_of_proj)), 'Parent', ax);

    %pause(1/FrameRate);
    pause(0);

    hue_of_proj = hue_of_proj + (hue_num_steps^-1);
end
close;

%% hue spectrum for tests
hue_spec_inhsv = zeros(1,256,3);
hue_spec_inhsv(:,:,1) = (0:255)/255;
hue_spec_inhsv(:,:,[2,3]) = 1;
hue_spec_inhsv = repmat(hue_spec_inhsv, 256, 1, 1);
hue_spec_inhsv = cat(2, hue_spec_inhsv, hue_spec_inhsv);
hue_spec_inrgb = hsv2rgb(hue_spec_inhsv);
Iorig = hue_spec_inrgb;

%% monochromeHSV3D
function mean_hue = HSV3D_mix(Iorig)
    % converting to hsv
    Iorig_hsv = rgb2hsv(Iorig);
    
    % converting hsv image to hsv color space coordinates (cartesian coordinates)
    Iorig_hsv3D = zeros(size(Iorig, 1), size(Iorig,2), 3);
    Iorig_hsv3D(:,:,3) = Iorig_hsv(:,:,3);
    Iorig_hsv3D(:,:,1) = Iorig_hsv(:,:,2) .* cos(Iorig_hsv(:,:,1) * 2*pi);
    Iorig_hsv3D(:,:,2) = Iorig_hsv(:,:,2) .* sin(Iorig_hsv(:,:,1) * 2*pi);
    
    % mix of image in hsv color space system (mix = mean = average)
    mean_Iorig_mun_v = squeeze(mean(Iorig_hsv3D, [1 2]));
    
    % turning the mix, a single color, from hsv color space coordinates to hsv
    mean_Iorig_hsv = zeros(1,3); % place holder
    mean_Iorig_hsv(2) = sqrt(mean_Iorig_mun_v(1)^2 + mean_Iorig_mun_v(2)^2); % saturation
    mean_Iorig_hsv(3) = mean_Iorig_mun_v(3); % value
    if mean_Iorig_mun_v(1) == 0 && mean_Iorig_mun_v(2) >= 0 % hue
        mean_Iorig_hsv(1) = (pi/2) / (2*pi);
    elseif mean_Iorig_mun_v(1) == 0 && mean_Iorig_mun_v(2) < 0
        mean_Iorig_hsv(1) = (3*pi/2) / (2*pi);
    elseif mean_Iorig_mun_v(1) >= 0 && mean_Iorig_mun_v(2) == 0
        mean_Iorig_hsv(1) = (0) / (2*pi);
    elseif mean_Iorig_mun_v(1) < 0 && mean_Iorig_mun_v(2) == 0
        mean_Iorig_hsv(1) = (pi) / (2*pi);
    elseif mean_Iorig_mun_v(1) > 0 && mean_Iorig_mun_v(2) > 0
        mean_Iorig_hsv(1) = atan(mean_Iorig_mun_v(2)/mean_Iorig_mun_v(1));
        mean_Iorig_hsv(1) = mean_Iorig_hsv(1) / (2*pi);
    elseif mean_Iorig_mun_v(1) > 0 && mean_Iorig_mun_v(2) < 0
        mean_Iorig_hsv(1) = atan(mean_Iorig_mun_v(2)/mean_Iorig_mun_v(1));
        mean_Iorig_hsv(1) = (2*pi+mean_Iorig_hsv(1)) / (2*pi);
    elseif mean_Iorig_mun_v(1) < 0 && mean_Iorig_mun_v(2) > 0
        mean_Iorig_hsv(1) = atan(mean_Iorig_mun_v(2)/mean_Iorig_mun_v(1));
        mean_Iorig_hsv(1) = (pi + mean_Iorig_hsv(1)) / (2*pi);
    elseif mean_Iorig_mun_v(1) < 0 && mean_Iorig_mun_v(2) < 0
        mean_Iorig_hsv(1) = atan(mean_Iorig_mun_v(2)/mean_Iorig_mun_v(1));
        mean_Iorig_hsv(1) = (pi + mean_Iorig_hsv(1)) / (2*pi);
    end

    mean_hue = mean_Iorig_hsv(1);
end

function Iproj_rgb = monochromeHSV3D(Iorig, hue_of_proj)
    hue_of_proj = mod(hue_of_proj, 1);
    
    % converting to hsv
    Iorig_hsv = rgb2hsv(Iorig);
    
    % converting hsv image to hsv color space coordinates (cartesian coordinates)
    Iorig_hsv3D = zeros(size(Iorig, 1), size(Iorig,2), 3);
    Iorig_hsv3D(:,:,3) = Iorig_hsv(:,:,3);
    Iorig_hsv3D(:,:,1) = Iorig_hsv(:,:,2) .* cos(Iorig_hsv(:,:,1) * 2*pi);
    Iorig_hsv3D(:,:,2) = Iorig_hsv(:,:,2) .* sin(Iorig_hsv(:,:,1) * 2*pi);
    
    %
    Iproj_hsv3D = Iorig_hsv3D;
    Iproj_hsv3D(:,:,1) = sqrt(Iorig_hsv3D(:,:,1).^2 + Iorig_hsv3D(:,:,2).^2) .* ...
        cos(Iorig_hsv(:,:,1) * 2*pi - hue_of_proj * 2*pi) .* ...
        cos(hue_of_proj * 2*pi);
    Iproj_hsv3D(:,:,2) = sqrt(Iorig_hsv3D(:,:,1).^2 + Iorig_hsv3D(:,:,2).^2) .* ...
        cos(Iorig_hsv(:,:,1) * 2*pi - hue_of_proj * 2*pi) .* ...
        sin(hue_of_proj * 2*pi);
    
    %
    Iproj_hsv3D = reshape(Iproj_hsv3D, [], 3);
    
    %
    Iproj_hsv = zeros(size(Iorig, 1), size(Iorig,2), 3);
    Iproj_hsv = reshape(Iproj_hsv, [], 3);
    
    %
    i = 1;
    Iproj_hsv_p = zeros(1,3); % place holder
    for i = 1:size(Iproj_hsv3D,1)
        this_hsv3D = Iproj_hsv3D(i, :);
        Iproj_hsv_p(2) = sqrt(this_hsv3D(1)^2 + this_hsv3D(2)^2); % saturation
        Iproj_hsv_p(3) = this_hsv3D(3); % value
        if this_hsv3D(1) == 0 && this_hsv3D(2) >= 0 % hue
            Iproj_hsv_p(1) = (pi/2) / (2*pi);
        elseif this_hsv3D(1) == 0 && this_hsv3D(2) < 0
            Iproj_hsv_p(1) = (3*pi/2) / (2*pi);
        elseif this_hsv3D(1) >= 0 && this_hsv3D(2) == 0
            Iproj_hsv_p(1) = (0) / (2*pi);
        elseif this_hsv3D(1) < 0 && this_hsv3D(2) == 0
            Iproj_hsv_p(1) = (pi) / (2*pi);
        elseif this_hsv3D(1) > 0 && this_hsv3D(2) > 0
            Iproj_hsv_p(1) = atan(this_hsv3D(2)/this_hsv3D(1));
            Iproj_hsv_p(1) = Iproj_hsv_p(1) / (2*pi);
        elseif this_hsv3D(1) > 0 && this_hsv3D(2) < 0
            Iproj_hsv_p(1) = atan(this_hsv3D(2)/this_hsv3D(1));
            Iproj_hsv_p(1) = (2*pi+Iproj_hsv_p(1)) / (2*pi);
        elseif this_hsv3D(1) < 0 && this_hsv3D(2) > 0
            Iproj_hsv_p(1) = atan(this_hsv3D(2)/this_hsv3D(1));
            Iproj_hsv_p(1) = (pi + Iproj_hsv_p(1)) / (2*pi);
        elseif this_hsv3D(1) < 0 && this_hsv3D(2) < 0
            Iproj_hsv_p(1) = atan(this_hsv3D(2)/this_hsv3D(1));
            Iproj_hsv_p(1) = (pi + Iproj_hsv_p(1)) / (2*pi);
        end
    
        Iproj_hsv (i, :) = Iproj_hsv_p;
        i = i + 1;
    end
    
    % the hue of the monochromatic image is not exactly equal to hue_of_proj
    % so epsilon is introduced as a threshold for the conversion error
    epsilon = 1/256;
    Iproj_hsv_n = abs(Iproj_hsv(:,1) - hue_of_proj) > epsilon;
    Iproj_hsv(Iproj_hsv_n, 1) = hue_of_proj;
    Iproj_hsv(Iproj_hsv_n, 2) = 0;
    
    %
    Iproj_hsv = reshape(Iproj_hsv, size(Iorig, 1), size(Iorig, 2), 3);
    Iproj_rgb = hsv2rgb(Iproj_hsv);

end