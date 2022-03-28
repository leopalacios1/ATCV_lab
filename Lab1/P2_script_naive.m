%% DOUBTS AND OBSERVATIONSforma
% Using only one BRIEF point pattern for each path/image
% No key Points detection, using whole patch... naive approach
%% SET UP _____________________________________________________________________________
clear; clc; 
currentPath = pwd;
f           = filesep;
ImagesPath  = strcat(currentPath,f,"traffic models/Meta/",f);

% Load Traffic Images
nImages = 43;
Ims = cell(1,nImages);           % library for imaes
for i = 1:nImages
    Im = imread(strcat(ImagesPath,int2str(i-1), ".png"));
    Ims{i} = Im;
end

%% SMOOTH _____________________________________________________________________________

% Parameters
sigma = 2; %(Filter Size is then 9x9)

% Gaussian Filter Images

Ims_G = cell(1,nImages); 
for i = 1:nImages
    Ims_G{i} = imgaussfilt(rgb2gray(Ims{i}),sigma);
end

%% BRIEF PATTERN SELECTION/CREATION ____________________________________________________

% Parameters
nBits = 256;
im_center = 50;
S = 100;

% (II) i.i.d. Gaussian(0, 1/25*S^2) (From paper: sigma = 1/5*S)
X = zeros(nBits,1);
count = 1;
while count<=nBits
    point = 50 + round(randn(1)*1/5*S);
    if point >= 1 && point <= 100
        X(count) = point;
        count = count +1;
    end
end

Y = zeros(nBits,1);
count = 1;
while count<=nBits
    point = 50 + round(randn(1)*1/5*S);
    if point >= 1 && point <= 100
        Y(count) = point;
        count = count +1;
    end
end
POINTS = [X,Y]; % NOTE: X and Y will be used as the pixel coordinates u,v of a given image. (X>Y)

% (IV) coarse polar grid

%% IMAGE DESCRIPTORS ____________________________________________________________

Ims_D = zeros(nImages,nBits); % Descriptor Array: Row for Image, Columns for bits
for i = 1:nImages
    Ims_D(i,:) = Ims_G{i}(X) > Ims_G{i}(Y);
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
% Displacy compactly (Taken from google)
fprintf([repmat(sprintf('%% %dd',max(floor(log10(abs(DistanceMTX(:)))))+2+any(DistanceMTX(:)<0)),1,size(DistanceMTX,2)) '\n'],DistanceMTX');
%% RESULTS MATRIX _______________________________________________________________

