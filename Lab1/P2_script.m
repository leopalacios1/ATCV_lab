%% DOUBTS AND OBSERVATIONSforma

% No key Points detection, using whole patch...
%% SET UP _____________________________________________________________________________
clear; clc; 
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

%% BRIEF PATTERN SELECTION/CREATION ____________________________________________________

% Parameters
nBits = 1024;
im_center = 50;
S = 100;        % Entire Image is our patch


% (II) i.i.d. Gaussian(0, 1/25*S^2) (From paper: sigma = 1/5*S)
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



%% (IV) coarse polar grid
X = [];
Y = [];
for i=1:round(nBits/2)
    theta = rand()*360;
    r1    = rand()*100-50;
    r2    = rand()*100-50;
    X     = [X, 50+cosd(theta)*r1];
end

POINTS = [X,Y]; % NOTE: X and Y will be used as the pixel coordinates u,v of a given image. (X>Y)

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
% Display compactly (Taken from google)
fprintf([repmat(sprintf('%% %dd',max(floor(log10(abs(DistanceMTX(:)))))+2+any(DistanceMTX(:)<0)),1,size(DistanceMTX,2)) '\n'],DistanceMTX');
%% PRINT RESULTS _______________________________________________________________


for i = 1:nImages
    tempResult    = DistanceMTX(i,:);
    % Convert 0 into NaN, hard fix to avoid selecting same image
    tempResult(tempResult==0)=NaN;
    [closesDist, closestIm]    = min(tempResult);
    [furthestDist, furthestIm] = max(tempResult);
    disp(['Image ', num2str(i-1), ':'])

    disp(['Closest classification is image ', num2str(closestIm-1), ...
        ', with a hamming distance of ', num2str(closesDist)])
    disp(['Furthest classification is image ', num2str(furthestIm-1), ...
        ', with a hamming distance of ', num2str(furthestDist)])
    disp(' ')

end
