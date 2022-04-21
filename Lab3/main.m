clear; clc; close all; 

Im = imread( "images/C arm -45.png");

%%
% circle_color = impixel(Im);
w = size(Im, 1);
h = size(Im, 2);
th = 0.2;

circle_color = [165.6667   45.0000   38.2222];

distance_matrix = vecnorm(   double(reshape(Im, [w*h, 3]))  - repmat(circle_color, [h*w,1]   ), 2, 2);
BW_im = reshape(distance_matrix/max(distance_matrix, [],'all'),  [w,h]);

imshow( BW_im < th)


stats = regionprops('table',BW_im < th ,'Centroid','Area','Circularity')


[M,i] = max(stats.Area);


imshow(Im)
hold on
x = stats.Centroid(i,1);

y = stats.Centroid(i,2);
plot(x,y, '+')