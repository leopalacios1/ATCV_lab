clear; clc; close all;

%% 

Im = imread( "color pictures/flowers.tif");
% Im = imread( "color pictures/papillary_carciroma.png");
% Im = imread( "color pictures/peppers_color.tif");
% Im = imread( "color pictures/tucan1.jpg");

X = double(reshape(Im, [size(Im,1)*size(Im,2), size(Im,3)]));

k = 7;

[indx, centers] = k_means(X, k, @eucl_distance);

Y = zeros(size(X));
for i = 1:k
    Y(indx == i,:) = repmat(centers(i,:), [sum(indx == i),1] );
end

Y = uint8(reshape(Y,[size(Im,1),size(Im,2),size(Im,3)]));

figure;
subplot(1,2,1)
imshow(Im)
title("Original image")
subplot(1,2,2)
imshow(Y)
title(strcat(" Clusterized k = ", num2str(k)))