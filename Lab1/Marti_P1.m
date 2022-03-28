

I1=rgb2gray(imread('img1.ppm'));I2=rgb2gray(imread('img2.ppm'));I3=rgb2gray(imread('img3.ppm'));I4=rgb2gray(imread('img4.ppm'));I5=rgb2gray(imread('img5.ppm'));I6=rgb2gray(imread('img6.ppm'));
N=1200;
features = detectFASTFeatures(I1);
keypoints_features=features.selectStrongest(N)

imshow(I1); 
hold on;
plot(keypoints_features);

%Homography
load('H1to2p');
load('H1to3p');
load('H1to4p');
load('H1to5p');
load('H1to6p');
Hom=zeros(3,3,6);
Hom(:,:,2)=H1to2p;
Hom(:,:,3)=H1to3p; 
Hom(:,:,4)=H1to4p; 
Hom(:,:,5)=H1to5p; 
Hom(:,:,6)=H1to6p;

% Transformation
Points=zeros(N,2,6);
Points(:,:,1)=keypoints_features.Location;
for j=2:6
    for i=1:N
    norm=Hom(3,1,j)*keypoints_features.Location(i,1)+Hom(3,2,j)*keypoints_features.Location(i,2)+1;
    Points(i,1,j)=floor((Hom(1,1,j)*keypoints_features.Location(i,1)+Hom(1,2,j)*keypoints_features.Location(i,2)+Hom(1,3,j))/norm);
    Points(i,2,j)=floor((Hom(2,1,j)*keypoints_features.Location(i,1)+Hom(2,2,j)*keypoints_features.Location(i,2)+Hom(2,3,j))/norm);
    end
end

%% Patern and wall features

[X,Y] = Pattern();

F1=Features(I1,Points(:,:,1),X,Y);
F2=Features(I2,Points(:,:,2),X,Y);
F3=Features(I3,Points(:,:,3),X,Y);
F4=Features(I4,Points(:,:,4),X,Y);
F5=Features(I5,Points(:,:,5),X,Y);
F6=Features(I6,Points(:,:,6),X,Y);



%% Hamming distance

Classifier=zeros(N,1);
C=[];
W=[];

Wall=F2;
for j=1:N  
    V=F1(j,:);
    [minDistance,index] = HammingDistance(V,Wall);
    if index==j
        Classifier(j)=1; 
        C=[C;minDistance];
    else
         W=[W;minDistance];
    end

end

Recognition_Rate=sum(Classifier(:))/N*100
%% Histogram of the error
figure;
h1 = histogram(C,33,'Normalization','pdf')
hold on
b1=h1.Values;
h2 = histogram(W,33,'Normalization','pdf')
b2=h2.Values;

%% ROC CURVE
figure
plot(cumsum(b2),cumsum(b1));grid
%%
function  [result] = BRIEFDescriptor(img,rndX,rndY)
    
    result = zeros(128,1);
    img = imgaussfilt(img,2);
    for i = 1:128
        
        x1 = rndX(i,1);
        y1 = rndY(i,1);
        x2 = rndX(i,2);
        y2 = rndY(i,2);
        if img(x1,y1)< img(x2,y2)
            result(i,1) = 1;
        else
            result(i,1) = 0;
        end


   end
end

function [M] = Features(I,points,rndX,rndY)

    Img= padarray(I,[8 8],0,'both');
    points(points<0)=1;
    points=points+8;
    M=zeros(size(points,1),128);
    X=points(:,1);
    Y=points(:,2);
    X(X>size(I,2))=size(I,2)-1;
    Y(Y>size(I,1))=size(I,1)-1;
   
    for i=1:size(points,1) 
      
    SubImg=Img(Y(i)-8:Y(i)+8,X(i)-8:X(i)+8);
    M(i,:) = BRIEFDescriptor(SubImg,rndX,rndY)';

    end

end

function [minDistance,index] = HammingDistance(V,Wall)
    N=size(Wall,1);
    DistanceVector=zeros(N,1);
    for i=1:N
        DistanceVector(i)=sum(xor(V,Wall(i,:)));
    end   
    [minDistance,index]=min(DistanceVector);
       
end

function [X,Y] = Pattern()
    m=16;
    n=16;

    rng(1)
    %Random part
    X = round(randi([1,n],128,2));
    Y = round(randi([1,m],128,2));
    
    X(X==0)=1;Y(Y==0)=1; % remove img(0,0), MATLAB starts at 1 the vectors
   
       % Plot positions
    figure;hold on;grid on
    for i = 1:128
        plot([X(i,1) Y(i,1)],[X(i,2) Y(i,2)],'Color','black')
        hold on
    end

end


