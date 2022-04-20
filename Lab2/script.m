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
%% 1. Euclidean distance
X = double(reshape(Im, [size(Im,1)*size(Im,2), size(Im,3)]));

k    = 5;
dist = @eucl_distance; %@mahal; 

[indx, centers] = k_means(X, k, dist, @RGB_mean);

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
%% 1.1 Visualize Image

h = size(Y,1);
w = size(Y,2);
v = h*w;

% Color Clustering from original
C = double([reshape(Y(:,:,1),[v 1]), reshape(Y(:,:,2),[v 1]), reshape(Y(:,:,3),[v 1])])./255;

figure;
scatter3(reshape(Im(:,:,1),[v 1]),reshape(Im(:,:,2),[v 1]),reshape(Im(:,:,3),[v 1]),36,C,'filled')

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
figure; imshow(Palette); title('Color Segmentation Results')

%% Morphology 

% Select color from the color palette
colorIdx = 5;
Color = Colors(colorIdx,:);
Y_bin = Y(:,:,1) == Color(1) & Y(:,:,2) == Color(2) & Y(:,:,3) == Color(3);
imshow(Y_bin)
%% alternative

carcicoma_blue = [112, 56, 129];
dist_vec = dist( carcicoma_blue, centers);
[m, carcicoma_idx] = min(dist_vec);
% use idx of the segmentation

BW_image = zeros(h,w);
BW_image(indx == carcicoma_idx) = 1;
figure;
imshow(BW_image)

figure;
se = strel('disk',3);
BW_im2 = imopen(BW_image, se);
se = strel('disk',10);
BW_im2 = imdilate(BW_im2, se);
BW_im2 = imdilate(BW_im2, se);
BW_im2 = imfill(BW_im2);
imshow(BW_im2)

figure;
Im_2 = reshape(Im, [h*w, 3]) ;
bool_indices = reshape(BW_im2, [h*w,1]) > 0.5;
Im_2( bool_indices , :  ) = repmat(carcicoma_blue, [sum(bool_indices),1]);
Im_2 = reshape( Im_2, [h,w,3] );
imshow(Im_2)

%% HSL space

X = double(reshape(Im, [size(Im,1)*size(Im,2), size(Im,3)]));
XX = rgb2HSL(X);

k    = 5;
dist = @HSL_dist;

[indx, centers] = k_means(XX, k, @HSL_dist, @HSL_mean);
RGB_centers = HSL2rgb(centers);

Y = zeros(size(X));
for i = 1:k
    Y(indx == i,:) = repmat( RGB_centers(i,:) , [sum(indx == i),1] );
end

Y = uint8(  reshape(Y,  [size(Im,1),size(Im,2),size(Im,3)])  );

figure;
subplot(1,2,1)
imshow(Im) 
title("Original image")
subplot(1,2,2)
imshow(Y)
title(strcat(" Clusterized k = ", num2str(k)))

%% HSL color space

h = size(Y,1);
w = size(Y,2);
v = h*w;
close all

% Color Clustering from original
C = double([reshape(Im(:,:,1),[v 1]), reshape(Im(:,:,2),[v 1]), reshape(Im(:,:,3),[v 1])])./255;

figure;
scatter3( XX(:,2).*cos(XX(:,1))  ,XX(:,2).*sin(XX(:,1)) ,XX(:,3),36,C,'filled')

xlim([-1,1])
ylim([-1,1])
zlim([0,1])

% Color Clustering from original
CC = double([reshape(Y(:,:,1),[v 1]), reshape(Y(:,:,2),[v 1]), reshape(Y(:,:,3),[v 1])])./255;

figure;
scatter3( XX(:,2).*cos(XX(:,1))  ,XX(:,2).*sin(XX(:,1)) ,XX(:,3),36,CC,'filled')
xlim([-1,1])
ylim([-1,1])
zlim([0,1])


%% Carcicoma morphology

carcicoma_blue = [112, 56, 129];
dist_vec = dist( rgb2HSL(carcicoma_blue) , centers);
[m, carcicoma_idx] = min(dist_vec);
% use idx of the segmentation

BW_image = zeros(h,w);
BW_image(indx == carcicoma_idx) = 1;
figure;
imshow(BW_image)

figure;
se = strel('disk',6);
BW_im2 = imopen(BW_image, se);
se = strel('disk',7);
BW_im2 = imdilate(BW_im2, se);
BW_im2 = imdilate(BW_im2, se);
BW_im2 = imfill(BW_im2);
imshow(BW_im2)

figure;
Im_2 = reshape(Im, [h*w, 3]) ;
bool_indices = reshape(BW_im2, [h*w,1]) > 0.5;
Im_2( bool_indices , :  ) = repmat(carcicoma_blue, [sum(bool_indices),1]);
Im_2 = reshape( Im_2, [h,w,3] );
imshow(Im_2)



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

