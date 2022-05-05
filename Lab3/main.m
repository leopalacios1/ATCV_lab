clear; clc; close all; 

%%
Im = imread( "images/C arm -90.png");

% ["images/C arm 90.png", "images/C arm 60.png", "images/C arm 45.png",
%  "images/C arm 30.png", "images/C arm 0.png", "images/C arm -30.png",
%  "images/C arm -45.png", "images/C arm -60.png", "images/C arm -90.png"];
%
% circle_color = impixel(Im);
w = size(Im, 1);
h = size(Im, 2);
th = 0.15;  
% THRESHOLDS:
% 0,90: 0.25
% 60:   0.2
% m90:  0.15
area_th = 500;

circle_color = [191.5593   63.1649   59.9497]; % [165.6667   45.0000   38.2222];
%circle_color = [165.6667   45.0000   38.2222];
circle_color = [122.4832   31.8623   21.5139]; % 90,m60,m90
%circle_color = [164.9948   45.4730   38.6707]; % 0
%circle_color = [220.3644   48.2486   54.4132]; %-30
%circle_color = [233.8445  122.8051  132.2756] %-45 [BAD]

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
names = ["images/C arm 90.png", "images/C arm 60.png", "images/C arm 45.png", ...
        "images/C arm 30.png", "images/C arm 0.png", "images/C arm -30.png", ...
        "images/C arm -45.png", "images/C arm -60.png", "images/C arm -90.png"];

point_colors = [];
for i = 1%:9
    names(i)
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


%% GET CAMERA POSE (USING estimateWorldCameraPose(): P3P algorithm)

% estimateWorldCameraPose Estimate camera pose from 3-D to 2-D point correspondences
%     [worldOrientation, worldLocation] = estimateWorldCameraPose(imagePoints, worldPoints, cameraParams)
%     returns the orientation and location of a calibrated camera in the world coordinate
%     system in which worldPoints are defined. 
%   
%     The function solves the Perspective-n-Point (PnP) problem using the P3P algorithm. 
%     The function eliminates spurious correspondences using the M-estimator SAmple Consensus 
%     (MSAC) algorithm
%   
%     Inputs            Description
%     ------            -----------
%     imagePoints       M-by-2 array of [x,y] coordinates of undistorted image points,
%                       with M >= 4.  
%    
%     worldPoints       M-by-3 array of [x,y,z] coordinates of world points.
%    
%     cameraParams      a cameraParameters or cameraIntrinsics object

% PROBLEM PARAMETERS
radiusRot     = 0.745; % m (Estimated)
angles        = [90,60,45,30,0,-30,-45,-60,-90];  % Degrees (Exact)

% CIRCLE CENTROIDS FOR EACH ANGLE (FROM PREVIOUS RESULTS)
C90  = [1024.3 , 532.78]; % Surpisingly Decent
C60  = [1065.4 , 390.99]; % Decent
C45  = [1111   , 328.99]; % Not That Great
C30  = [1175.7 , 282.47]; % Decent
C00  = [1366.1 ,  235.6]; % Good
Cm30 = [1560.4 , 289.39]; % Good
Cm45 = [1113.6 , 331.15]; % Good
Cm60 = [1699.1 , 442.75]; % Good
Cm90 = [1721.2 , 633.76]; % Good
idxBEST = [5,6,7,8,9];
imagePointsAll  = [C90;C60;C45;C30;C00;Cm30;Cm45;Cm60;Cm90];
imagePointsBest = imagePointsAll(idxBEST);

% ASSUME CIRCLES LIES IN A PLANE PERPENDICULAR TO THE FLOOR IN A REFERENCE
% FRAME WERE THE CENTER OF ROTATION IS THE ORIGIN. HENCE, THE CENTROIDS
% WILL Z=0 AND X AND Y CAN BE KNOWN FROM THE ANGLE AND THE RADIUS OF
% ROTATION. (SEE PDF IMAGE)

worldPointsAll  = [sind(-angles'),cosd(angles'),zeros(length(angles),1)];
worldPointsBest = worldPointsAll(idxBEST);

% Camera Params
Xres = 1920;
Yres = 1080;
IntrinsicMatrix = [1060,0,Xres/2   ; 0,1060,Yres/2   ;  0,0,1]'; % NOTE: TRANSPOSED
cameraParams = cameraParameters('IntrinsicMatrix',IntrinsicMatrix);
[worldOrientation, worldLocation] = estimateWorldCameraPose(imagePointsAll, worldPointsAll, cameraParams)
%[worldOrientation, worldLocation] = estimateWorldCameraPose(imagePointsBest, worldPointsBest, cameraParams)

%% OBTAINED CAMERA POSE TO transformation from world coordinates to camera coordinates.
% Create rigid3d object:
cameraPose = rigid3d(worldOrientation,worldLocation);
% see tform and rigid3d documentation
tform = cameraPoseToExtrinsics(cameraPose);

% THE TRANSFORMATION IS THEN: [x y z] = [X Y Z]*R + t

%% DISPLAY OBTAINED ROTATION CENTER ON IMAGE

origin   = [0 0 0];
originIm = origin*tform.Rotation + tform.Translation;
pixels   = IntrinsicMatrix'*originIm';
pixels   = pixels./pixels(3)

Im = imread( "images/C arm 0.png");
imshow(Im)
hold on;
plot(pixels(1),pixels(2), 'r+','MarkerSize',20)











