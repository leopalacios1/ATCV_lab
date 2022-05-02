clear; clc; close all; 

%%
Im = imread( "images/C arm -60.png");

%%
% circle_color = impixel(Im);
w = size(Im, 1);
h = size(Im, 2);
th = 0.15;
area_th = 500;

circle_color = [191.5593   63.1649   59.9497]; % [165.6667   45.0000   38.2222];

% lab space distance
circle_color_lab = rgb2lab(circle_color);
lab_im = rgb2lab(Im);
lab_im = reshape(lab_im, [w*h, 3]);

distance_matrix = lab_dist(lab_im, circle_color_lab);
max_dist = max(distance_matrix, [],'all');
min_dist = min(distance_matrix, [],'all');
BW_im = reshape( (distance_matrix-min_dist)/(max_dist-min_dist),  [w,h]);


bin_image =  BW_im < th;
SE = strel('disk',2);
% bin_image = imclose(bin_image, SE);
% imshow( bin_image)

stats = regionprops('table',bin_image ,'Centroid','Area','Circularity');

stats = stats( stats.Area > area_th,:)
[M,i] = max(stats.Circularity);
% the circle will have the centroid in the maximum circularity
centroid = stats(i,:).Centroid;

% show results
imshow( bin_image)
hold on
plot(centroid(1), centroid(2), 'r+')
hold off
%% get the circle color with lab space
names = ["images/C arm 90.png", "images/C arm 60.png", "images/C arm 45.png",
        "images/C arm 30.png", "images/C arm 0.png", "images/C arm -30.png",
        "images/C arm -45.png", "images/C arm -60.png", "images/C arm -90.png"];

point_colors = [];
for i = 1:9
    Im = imread( names(i) );
    point_colors = [point_colors; impixel(Im)];
end
save('point_colors.mat','point_colors')
lab_points = rgb2lab(point_colors);
circle_color_lab = mean(lab_points, 1);
circle_color = lab2rgb(circle_color_lab);

%%

imshow(Im)
hold on
x = stats.Centroid(i,1);

y = stats.Centroid(i,2);
plot(x,y, '+')