%% DOUBTS AND OBSERVATIONSforma

% No key Points detection, using whole patch...
%% SET UP _____________________________________________________________________________
clear; clc; 
displaying   = 1;
currentPath = pwd;
f           = filesep;
ImagesPath  = strcat(currentPath,f,"traffic models/Meta/",f);

% Load Traffic Images
nImages = 43;
Ims = cell(1,nImages);           % library for images
for i = 1:nImages
    Im = imread(strcat(ImagesPath,int2str(i-1), ".png"));
    Ims{i} = Im;
end

%% RESIZE IMAGES TO 100x100 PIXELS ____________________________________________________

for i = 1:nImages
    Ims{i} = imresize(Ims{i},[100 100]);
end

%% SMOOTH - GAUSSIAN FILTER ___________________________________________________________
% Parameters
sigma = 2; %(Filter Size is then 9x9)

% Gaussian Filter Images

Ims_G = cell(1,nImages); 
for i = 1:nImages
    Ims_G{i} = imgaussfilt(rgb2gray(Ims{i}),sigma);
end
%% HERE FOR TESTING AND AVERAGING ACROSS RANDOM PATTERNS FOR EVALUATION
clc
nTrials = 1000;
displaying = 0;
DISTminmax_norm_mean = zeros(nImages,2);
for i =1:nTrials
%% BRIEF PATTERN SELECTION/CREATION ____________________________________________________

% Parameters
nBits = 1024; %128,256,512,1024
im_center = 50;
S = 100;        % Entire Image is our patch

%% (A) i.i.d. Gaussian(0, 1/25*S^2) (From paper: sigma = 1/5*S)
POINTS = zeros(nBits,4); 
nCoord = 4*nBits;   % For each bite in the descriptor, 2 points of 2 coordinates are used
count = 1;
pointMin = 1;       % Select points inside image
pointMax = 100;      % Smallest image dimension, could have also resized images
while count<=nCoord
    point = im_center + round(randn(1)*1/5*S);
    if point >= pointMin && point <= pointMax
        POINTS(count) = point;
        count = count +1;
    end
end
% Change from coordinates (u,v) to linear index. New Point MTX is then nBits by 2
POINTS_lin = zeros(nBits,2);
for pair = 1:nBits
    POINTS_lin(pair,1) = sub2ind([S S],POINTS(pair,1),POINTS(pair,2));
    POINTS_lin(pair,2) = sub2ind([S S],POINTS(pair,3),POINTS(pair,4));
end



%% (B) Random polar pairs
% nBits = 1024; %128,256,512,1024
% Psize = 100;
% POINTS = zeros(nBits,4); 
% for i=1:nBits
%     theta       = rand()*360;
%     r1          = rand()*100-50;
%     r2          = rand()*100-50;
%     POINTS(i,1) = ceil(50+cosd(theta)*r1); % x1
%     POINTS(i,2) = ceil(50+sind(theta)*r1); % y1
%     POINTS(i,3) = ceil(50+cosd(theta)*r2); % x2
%     POINTS(i,4) = ceil(50+sind(theta)*r2); % y2
% end
% if displaying
%     figure; hold on;
%     for i = 1:nBits
%         plot( [POINTS(i,1);POINTS(i,2)], [POINTS(i,3);POINTS(i,4)], '-k')
%     end
% end
% % Change from coordinates (u,v) to linear index. New Point MTX is then nBits by 2
% POINTS_lin = zeros(nBits,2);
% for pair = 1:nBits
%     POINTS_lin(pair,1) = sub2ind([Psize Psize],POINTS(pair,1),POINTS(pair,2));
%     POINTS_lin(pair,2) = sub2ind([Psize Psize],POINTS(pair,3),POINTS(pair,4));
% end

%% IMAGE DESCRIPTORS ____________________________________________________________

Ims_D = zeros(nImages,nBits); % Descriptor Array: Row for Image, Columns for bits
for i = 1:nImages
    Ims_D(i,:) = Ims_G{i}(POINTS_lin(:,1)) > Ims_G{i}(POINTS_lin(:,2));
end


%% EVALUATION ___________________________________________________________________

DistanceMTX = zeros(nImages,nImages);

for i= 1:nImages
    % Image Descriptor = Ims_D(i,:)
    DescriptorMTX    = repmat(Ims_D(i,:),nImages,1);
    ComparisonMTX    = xor(Ims_D,DescriptorMTX);
    DistanceVec      = sum(ComparisonMTX,2)';
    DistanceMTX(i,:) = DistanceVec;
end
DistanceMTX_norm     = round(DistanceMTX./nBits*100);
if displaying
    % Display compactly (Taken from google)
    fprintf([repmat(sprintf('%% %dd',max(floor(log10(abs(DistanceMTX(:)))))+2+any(DistanceMTX(:)<0)),1,size(DistanceMTX,2)) '\n'],DistanceMTX');
    disp('____________________________________________________________________________________________________________________________________________________________________________')
    fprintf([repmat(sprintf('%% %dd',max(floor(log10(abs(DistanceMTX_norm(:)))))+2+any(DistanceMTX_norm(:)<0)),1,size(DistanceMTX_norm,2)) '\n'],DistanceMTX_norm');
end
%% RESULTS _______________________________________________________________

DISTminmax = zeros(nImages,2); % Rows: model, Columns: Min Dist, Max Dist

for i = 1:nImages
    tempResult    = DistanceMTX(i,:);
    % Convert 0 into NaN, hard fix to avoid selecting same image
    tempResult(tempResult==0)=NaN;
    [closesDist, closestIm]    = min(tempResult);
    [furthestDist, furthestIm] = max(tempResult);
    DISTminmax(i,:) = [closesDist, furthestDist];
    if displaying
        disp(['Image ', num2str(i-1), ':'])
    
        disp(['Closest classification is image ', num2str(closestIm-1), ...
            ', with a hamming distance of ', num2str(closesDist)])
        disp(['Furthest classification is image ', num2str(furthestIm-1), ...
            ', with a hamming distance of ', num2str(furthestDist)])
        disp(' ')
    end
end
DISTminmax_norm = round(DISTminmax./ nBits*100,2);


DISTminmax_norm_mean = DISTminmax_norm_mean + DISTminmax_norm;

%% END of TESTING
end
DISTminmax_norm_mean = DISTminmax_norm_mean./nTrials;
nBits
DISTminmax_norm_mean
