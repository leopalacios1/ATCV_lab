clear; clc; close all; 

Im = imread( "color pictures/flowers.tif");
% Im = imread( "color pictures/papillary_carciroma.png");
% Im = imread( "color pictures/peppers_color.tif");
% Im = imread( "color pictures/tucan1.jpg");

%% 1. Euclidean distance
X = double(reshape(Im, [size(Im,1)*size(Im,2), size(Im,3)]));

k = 5;

[indx, centers] = k_means(X, k, @eucl_distance, @RGB_mean);

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

%% HSL space

X = double(reshape(Im, [size(Im,1)*size(Im,2), size(Im,3)]));
XX = rgb2HSL(X);

k = 5;

[indx, centers] = k_means(XX, k, @HSL_dist, @HSL_mean);
centers = HSL2rgb(centers);

Y = zeros(size(X));
for i = 1:k
    Y(indx == i,:) = repmat( centers(i,:) , [sum(indx == i),1] );
end

Y = uint8(  reshape(Y,  [size(Im,1),size(Im,2),size(Im,3)])  );

figure;
subplot(1,2,1)
imshow(Im) 
title("Original image")
subplot(1,2,2)
imshow(Y)
title(strcat(" Clusterized k = ", num2str(k)))


%% HSL example

X = double(reshape(Im, [size(Im,1)*size(Im,2), size(Im,3)]));
XX = rgb2HSL(X);
Y = HSL2rgb(XX);
Y = uint8(  reshape(Y,  [size(Im,1),size(Im,2),size(Im,3)])  );

figure;
subplot(1,2,1)
imshow(Im) 
title("Original image")
subplot(1,2,2)
imshow(Y)
title(strcat(" Clusterized k = ", num2str(k)))


%% HSV distance

XX = rgb2hsv(Im);
X = double(reshape(XX, [size(Im,1)*size(Im,2), size(Im,3)]));

k = 5;

[indx, centers] = k_means(X, k, @HSV_dist, @HSV_mean);
centers = hsv2rgb([centers])*255;

Y = zeros(size(X));
for i = 1:k
    Y(indx == i,:) = repmat( centers(i,:) , [sum(indx == i),1] );
end

Y = uint8(  reshape(Y,  [size(Im,1),size(Im,2),size(Im,3)])  );

figure;
subplot(1,2,1)
imshow(Im) 
title("Original image")
subplot(1,2,2)
imshow(Y)
title(strcat(" Clusterized k = ", num2str(k)))


%% CIE LAB distance

XX = rgb2lab(Im);
X = double(reshape(XX, [size(Im,1)*size(Im,2), size(Im,3)]));

k = 5;

[indx, centers] = k_means(X, k, @lab_dist, @RGB_mean);
centers = lab2rgb([centers]);

Y = zeros(size(X));
for i = 1:k
    Y(indx == i,:) = repmat( centers(i,:) , [sum(indx == i),1] );
end

Y =   reshape(Y,  [size(Im,1),size(Im,2),size(Im,3)])  ;

figure;
subplot(1,2,1)
imshow(Im) 
title("Original image")
subplot(1,2,2)
imshow(Y)
title(strcat(" Clusterized k = ", num2str(k)))

