clear; clc; close all; 

% Im = imread( "color pictures/flowers.tif");
  Im = imread( "color pictures/papillary_carciroma.png");
% Im = imread( "color pictures/peppers_color.tif");
% Im = imread( "color pictures/tucan1.jpg");

%% 0. Visualize Image
h = size(Im,1);
w = size(Im,2);
v = h*w;

% COLORS
C = double([reshape(Im(:,:,1),[v 1]), reshape(Im(:,:,2),[v 1]), reshape(Im(:,:,3),[v 1])])./255;

figure;
scatter3(reshape(Im(:,:,1),[v 1]),reshape(Im(:,:,2),[v 1]),reshape(Im(:,:,3),[v 1]),36,C,'filled')
%
% Histograms
figure; 
histogram(Im(:,:,1)); title('Red values histogram')
figure; 
histogram(Im(:,:,2)); title('Blue values histogram')
figure; 
histogram(Im(:,:,3)); title('Green values histogram')
%% 1. RGB Space
X = double(reshape(Im, [size(Im,1)*size(Im,2), size(Im,3)]));

k    = 3;
% EUCLIDEAN 3
dist = @eucl_distance; %@cityblock_dist @eucl_distance ; 

[indx, centers] = k_means(X, k, dist, @RGB_mean);

Y = zeros(size(X));
for i = 1:k
    Y(indx == i,:) = repmat(centers(i,:), [sum(indx == i),1] );
end

Y_euc3 = uint8(reshape(Y,[size(Im,1),size(Im,2),size(Im,3)]));

% CITY BLOCK 3
dist = @cityblock_dist; %@cityblock_dist @eucl_distance ; 

[indx, centers] = k_means(X, k, dist, @RGB_mean);

Y = zeros(size(X));
for i = 1:k
    Y(indx == i,:) = repmat(centers(i,:), [sum(indx == i),1] );
end

Y_cb3 = uint8(reshape(Y,[size(Im,1),size(Im,2),size(Im,3)]));


k    = 4;
% EUCLIDEAN 4
dist = @eucl_distance; %@cityblock_dist @eucl_distance ; 

[indx, centers] = k_means(X, k, dist, @RGB_mean);

Y = zeros(size(X));
for i = 1:k
    Y(indx == i,:) = repmat(centers(i,:), [sum(indx == i),1] );
end

Y_euc4 = uint8(reshape(Y,[size(Im,1),size(Im,2),size(Im,3)]));

% CITY BLOCK 4
dist = @cityblock_dist; %@cityblock_dist @eucl_distance ; 

[indx, centers] = k_means(X, k, dist, @RGB_mean);

Y = zeros(size(X));
for i = 1:k
    Y(indx == i,:) = repmat(centers(i,:), [sum(indx == i),1] );
end

Y_cb4 = uint8(reshape(Y,[size(Im,1),size(Im,2),size(Im,3)]));


k    = 5;
% EUCLIDEAN 5
dist = @eucl_distance; %@cityblock_dist @eucl_distance ; 

[indx, centers] = k_means(X, k, dist, @RGB_mean);

Y = zeros(size(X));
for i = 1:k
    Y(indx == i,:) = repmat(centers(i,:), [sum(indx == i),1] );
end

Y_euc5 = uint8(reshape(Y,[size(Im,1),size(Im,2),size(Im,3)]));

% CITY BLOCK 5
dist = @cityblock_dist; %@cityblock_dist @eucl_distance ; 

[indx, centers] = k_means(X, k, dist, @RGB_mean);

Y = zeros(size(X));
for i = 1:k
    Y(indx == i,:) = repmat(centers(i,:), [sum(indx == i),1] );
end

Y_cb5 = uint8(reshape(Y,[size(Im,1),size(Im,2),size(Im,3)]));


figure;
subplot(2,3,1)
imshow(Y_euc3)
title(strcat(" k = 3, with euclidean distance."))
subplot(2,3,2)
imshow(Y_euc4)
title(strcat(" k = 4, with euclidean distance."))
subplot(2,3,3)
imshow(Y_euc5)
title(strcat(" k = 5, with euclidean distance."))

subplot(2,3,4)
imshow(Y_cb3)
title(strcat(" k = 3, with city block distance."))
subplot(2,3,5)
imshow(Y_cb4)
title(strcat(" k = 4, with city block distance."))
subplot(2,3,6)
imshow(Y_cb5)
title(strcat(" k = 5, with city block distance."))
%% Visualize Image
close all;
Y = Y_euc5; % Select best parameters

h = size(Y,1);
w = size(Y,2);
v = h*w;

% Color Clustering from original
C = double([reshape(Y(:,:,1),[v 1]), reshape(Y(:,:,2),[v 1]), reshape(Y(:,:,3),[v 1])])./255;

figure;
imshow(Y)
title('Image after Cluster')

figure;
scatter3(reshape(Im(:,:,1),[v 1]),reshape(Im(:,:,2),[v 1]),reshape(Im(:,:,3),[v 1]),36,C,'filled');
title('Original RGB values and their assigned cluster')

% Colors obtained
Y_reshaped = double(reshape(Y, [size(Y,1)*size(Y,2), size(Y,3)]));
Colors = unique(Y_reshaped,'rows');
nColors = k;
% Sort Colors
distToWhite = sum((repmat(255,size(Colors,1),size(Colors,2))-Colors).^2,2);
[B,I]       = sort(distToWhite);
Colors      = Colors(I,:);
% Creat column square array displaying colors
sizeSq = 100;
Palette = zeros(sizeSq*nColors,sizeSq*nColors,3);
for i=1:nColors
    Palette(  ((i-1)*sizeSq+1)  :  (i)*sizeSq,:,1) = Colors(i,1);
    Palette(  ((i-1)*sizeSq+1)  :  (i)*sizeSq,:,2) = Colors(i,2);
    Palette(  ((i-1)*sizeSq+1)  :  (i)*sizeSq,:,3) = Colors(i,3);
end
Palette = uint8(Palette);

figure; 
imshow(Palette);
title('Color Palette Obtained')

%% Morphology 
close all;
Y = Y_euc5;
% Select color from the color palette
colorIdx = 5;
Color = Colors(colorIdx,:);
Y_bin = ~(Y(:,:,1) == Color(1) & Y(:,:,2) == Color(2) & Y(:,:,3) == Color(3));
figure(1);
imshow(Y_bin)

% NOT SURE HOW TO BEST IMPLEMENT THIS PART

SE = strel('disk',2);
%SE = strel('arbitrary',ones(2));

Y_e    = imerode(Y_bin,SE);
Y_d    = imdilate(Y_bin,SE);
Y_dd   = imdilate(Y_d,SE);
Y_dde  = imerode(Y_dd,SE);
Y_ddee = imerode(Y_dde,SE);
figure(2);
imshow(imclose(Y_bin,SE));
%%
figure();
imshow(imnclose(imclose(Y_bin,SE),SE));

%% HSL space

X = double(reshape(Im, [size(Im,1)*size(Im,2), size(Im,3)]));
XX = rgb2HSL(X);

k    = 5;
dist = @HSL_dist;

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
clc

XX = rgb2hsv(double(Im)./255);
X = double(reshape(XX, [size(Im,1)*size(Im,2), size(Im,3)]));
k = 5;

[indx, centers] = k_means(X, k, @HSV_dis2, @HSV_mean2);
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


%% HSV & HSL analysis

% HSL 3
X = double(reshape(Im, [size(Im,1)*size(Im,2), size(Im,3)]));
XX = rgb2HSL(X);

k    = 3;
dist = @HSL_dist;

[indx, centers] = k_means(XX, k, @HSL_dist, @HSL_mean);
centers = HSL2rgb(centers);

Y = zeros(size(X));
for i = 1:k
    Y(indx == i,:) = repmat( centers(i,:) , [sum(indx == i),1] );
end

Y_hsl3 = uint8(  reshape(Y,  [size(Im,1),size(Im,2),size(Im,3)])  );

% HSL 4
X = double(reshape(Im, [size(Im,1)*size(Im,2), size(Im,3)]));
XX = rgb2HSL(X);

k    = 4;
dist = @HSL_dist;

[indx, centers] = k_means(XX, k, @HSL_dist, @HSL_mean);
centers = HSL2rgb(centers);

Y = zeros(size(X));
for i = 1:k
    Y(indx == i,:) = repmat( centers(i,:) , [sum(indx == i),1] );
end

Y_hsl4 = uint8(  reshape(Y,  [size(Im,1),size(Im,2),size(Im,3)])  );

% HSL 5
k    = 5;
dist = @HSL_dist;

[indx, centers] = k_means(XX, k, @HSL_dist, @HSL_mean);
centers = HSL2rgb(centers);

Y = zeros(size(X));
for i = 1:k
    Y(indx == i,:) = repmat( centers(i,:) , [sum(indx == i),1] );
end

Y_hsl5 = uint8(  reshape(Y,  [size(Im,1),size(Im,2),size(Im,3)])  );

% HSV 3
XX = rgb2hsv(double(Im)./255);
X = double(reshape(XX, [size(Im,1)*size(Im,2), size(Im,3)]));

k = 3;

[indx, centers] = k_means(X, k, @HSV_dis2, @HSV_mean2);
centers = hsv2rgb([centers])*255;

Y = zeros(size(X));
for i = 1:k
    Y(indx == i,:) = repmat( centers(i,:) , [sum(indx == i),1] );
end

Y_hsv3 = uint8(  reshape(Y,  [size(Im,1),size(Im,2),size(Im,3)])  );

% HSV 4

k = 4;

[indx, centers] = k_means(X, k, @HSV_dis2, @HSV_mean2);
centers = hsv2rgb([centers])*255;

Y = zeros(size(X));
for i = 1:k
    Y(indx == i,:) = repmat( centers(i,:) , [sum(indx == i),1] );
end

Y_hsv4 = uint8(  reshape(Y,  [size(Im,1),size(Im,2),size(Im,3)])  );

% HSV 5

k = 5;

[indx, centers] = k_means(X, k, @HSV_dis2, @HSV_mean2);
centers = hsv2rgb([centers])*255;

Y = zeros(size(X));
for i = 1:k
    Y(indx == i,:) = repmat( centers(i,:) , [sum(indx == i),1] );
end

Y_hsv5 = uint8(  reshape(Y,  [size(Im,1),size(Im,2),size(Im,3)])  );



figure;
subplot(2,3,1)
imshow(Y_hsl3)
title(strcat(" k = 3, in HSL color space."))
subplot(2,3,2)
imshow(Y_hsl4)
title(strcat(" k = 4, in HSL color space."))
subplot(2,3,3)
imshow(Y_hsl5)
title(strcat(" k = 5, in HSL color space."))

subplot(2,3,4)
imshow(Y_hsv3)
title(strcat(" k = 3, in HSV color space."))
subplot(2,3,5)
imshow(Y_hsv4)
title(strcat(" k = 4, in HSV color space."))
subplot(2,3,6)
imshow(Y_hsv5)
title(strcat(" k = 5, in HSV color space."))

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

