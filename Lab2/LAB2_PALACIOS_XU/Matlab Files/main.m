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
%Y = Y_euc5; % Select best parameters

h = size(Y,1);
w = size(Y,2);
v = h*w;

% Color Clustering from original
C = double([reshape(Y(:,:,1),[v 1]), reshape(Y(:,:,2),[v 1]), reshape(Y(:,:,3),[v 1])])./255; %REMOVE OR ADD ./255

figure;
imshow(Y)
title('Image after Cluster','FontSize',18)

figure;
scatter3(reshape(Im(:,:,1),[v 1]),reshape(Im(:,:,2),[v 1]),reshape(Im(:,:,3),[v 1]),36,C,'filled');
title('Original RGB values and their assigned cluster','FontSize',18)

% Colors obtained
Y_reshaped = double(reshape(Y, [size(Y,1)*size(Y,2), size(Y,3)]));
Colors = unique(Y_reshaped,'rows');
nColors = k;
%%
% Sort Colors
distToWhite = sum((repmat(255,size(Colors,1),size(Colors,2))-Colors).^2,2);
[B,I]       = sort(distToWhite);
TEMP        = Colors;
Colors      = Colors(I,:)
% Creat column square array displaying colors
sizeSq = 100;
Palette = zeros(sizeSq*nColors,sizeSq*nColors,3);
for i=1:nColors
    Palette(  ((i-1)*sizeSq+1)  :  (i)*sizeSq,:,1) = Colors(i,1);
    Palette(  ((i-1)*sizeSq+1)  :  (i)*sizeSq,:,2) = Colors(i,2);
    Palette(  ((i-1)*sizeSq+1)  :  (i)*sizeSq,:,3) = Colors(i,3);
end
Palette = uint8(Palette);
close all
figure; 
imshow(Palette);
title('Color Palette Obtained','FontSize',18)

%% Morphology 
close all;
Y = Y_euc5;
% Select color from the color palette
colorIdx = 5;
Color = Colors(colorIdx,:);
Y_bin = ~(Y(:,:,1) == Color(1) & Y(:,:,2) == Color(2) & Y(:,:,3) == Color(3));
figure(1);
imshow(Y_bin)
%% alternative


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

carcicoma_blue = [112, 56, 129];
dist_vec = dist( carcicoma_blue, centers);
[m, carcicoma_idx] = min(dist_vec);
% use idx of the segmentation
%%
BW_image = zeros(h,w);
BW_image(indx == carcicoma_idx) = 1;
figure;
imshow(BW_image)

figure;
se = strel('disk',3);
BW_im2 = imopen(BW_image, se);
se = strel('disk',4);
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


%% Carcicoma morphology HSL

carcicoma_blue = [112, 56, 129];
dist_vec = dist( rgb2HSL(carcicoma_blue) , centers);
[m, carcicoma_idx] = min(dist_vec);
% use idx of the segmentation
%%
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

%% DISPLAY

clc, close all
figure;
subplot(1,2,1)
imshow(Im) 
title("Original image",'FontSize',18)
subplot(1,2,2)
imshow(Y)
title(strcat(" Clusterized k = ", num2str(k)),'FontSize',18)

% Color Clustering from original
C = double([reshape(Y(:,:,1),[v 1]), reshape(Y(:,:,2),[v 1]), reshape(Y(:,:,3),[v 1])]);

figure;
scatter3(reshape(Im(:,:,1),[v 1]),reshape(Im(:,:,2),[v 1]),reshape(Im(:,:,3),[v 1]),36,C,'filled');
title('Original RGB values and their assigned cluster','FontSize',18)

%% Carcicoma morphology HSV
clc, close all;
carcicoma_blue = [112, 56, 129]./255;
dist_vec = dist( rgb2hsv(carcicoma_blue) , centers);
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

%% LAB Analysis

XX = rgb2lab(Im);
X = double(reshape(XX, [size(Im,1)*size(Im,2), size(Im,3)]));

% K =3
k = 3;

[indx, centers] = k_means(X, k, @lab_dist, @RGB_mean);
centers = lab2rgb([centers]);

Y = zeros(size(X));
for i = 1:k
    Y(indx == i,:) = repmat( centers(i,:) , [sum(indx == i),1] );
end

Y_lab3 =   reshape(Y,  [size(Im,1),size(Im,2),size(Im,3)])  ;

% K = 4
k = 4;

[indx, centers] = k_means(X, k, @lab_dist, @RGB_mean);
centers = lab2rgb([centers]);

Y = zeros(size(X));
for i = 1:k
    Y(indx == i,:) = repmat( centers(i,:) , [sum(indx == i),1] );
end

Y_lab4 =   reshape(Y,  [size(Im,1),size(Im,2),size(Im,3)])  ;

% K =5
k = 5;

[indx, centers] = k_means(X, k, @lab_dist, @RGB_mean);
centers = lab2rgb([centers]);

Y = zeros(size(X));
for i = 1:k
    Y(indx == i,:) = repmat( centers(i,:) , [sum(indx == i),1] );
end

Y_lab5 =   reshape(Y,  [size(Im,1),size(Im,2),size(Im,3)])  ;




figure;
subplot(1,3,1)
imshow(Y_lab3)
title(strcat(" k = 3, in CIE LAB color space."),'FontSize',18)
subplot(1,3,2)
imshow(Y_lab4)
title(strcat(" k = 4, in CIE LAB color space."),'FontSize',18)
subplot(1,3,3)
imshow(Y_lab5)
title(strcat(" k = 5, in CIE LAB color space."),'FontSize',18)


%% Carcicoma morphology LAB
clc, close all;
carcicoma_blue = [112, 56, 129]./255;
dist_vec = dist(carcicoma_blue , centers);
[m, carcicoma_idx] = min(dist_vec);
% use idx of the segmentation
BW_image = zeros(h,w);
BW_image(indx == carcicoma_idx) = 1;
figure;
imshow(BW_image)

figure;
se = strel('disk',5);
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
