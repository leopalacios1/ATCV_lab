clear; clc; close all;

S = 32;                     % patch size

Ims = {};           % library for imaes
for i = 1:6
    Im = imread(strcat("wall/img", int2str(i), ".ppm"));
    Ims{i} = rgb2gray(Im);
end

A_mat = {};         % library of homographic matrices
for i = 2:6
    FileName = strcat('wall/H1to', int2str(i), 'p');
    fid = fopen(FileName, 'r');
    C = textscan(fid, '%s');
    C = C{1};
    A = zeros(3,3); 
    for j = 1:9
        A(floor((j-1)/3)+1, mod(j-1,3)+1) = str2num(C{j});
    end
    A_mat{i} = A;
end

%% select points

N_points = 1000;            % Number of corners selected

points1 = detectHarrisFeatures(Ims{1});


valid_points1 = and(and(points1.Location(:,1) > S/2+1, points1.Location(:,2) > S/2+1),  ...
            and(points1.Location(:,1) < size(Ims{1},2)-S/2, points1.Location(:,2) < size(Ims{1},1)-S/2)); 
% valid points are those that have some distance to the edge

points1 = points1(valid_points1);               % subset of valid points
points1 = points1.selectStrongest(N_points);    % choose the strongest 

points_matrix = zeros(6,2,N_points);            % dimensions are (which_im, (x,y), Points)
points_matrix(1,:,:) = points1.Location';       % save the ones from the first image
valid_points_mat = zeros(6,N_points);           % valid points for (which_im, which_point)
valid_points_mat(1,:) = ones(1,N_points);       % we have already made sure that the first image where valid

p_use = ones(3,N_points);
p_use(1:2,:) = points_matrix(1,:,:);            % p_use are the points (x,y,1)' in a matrix
for i = 2:6
    p_aux = A_mat{i}*p_use;                     % hom transf for the image i
    points_matrix(i,1,:) = p_aux(1,:)./p_aux(3,:);
    points_matrix(i,2,:) = p_aux(2,:)./p_aux(3,:);
    valid_points_mat(i,:) = and(and(points_matrix(i,1,:) > S/2+1, points_matrix(i,2,:) > S/2+1),  ...
            and(points_matrix(i,1,:) < size(Ims{i},2)-S/2, points_matrix(i,2,:) < size(Ims{i},1)-S/2)); 
    % check which are valid in the image i
end


% r_ind = randi(N_points);
% 
figure()
imshow(Ims{1})
hold on
plot( squeeze(points_matrix(1,1,:)), squeeze(points_matrix(1,2,:)), '+g' )
% plot(points_matrix(1, 1, r_ind), points_matrix(1, 2, r_ind), 'ro')
hold off


figure()
imshow(Ims{3})
hold on
plot( squeeze(points_matrix(3,1,:)), squeeze(points_matrix(3,2,:)), '+g' )
% plot(points_matrix(3, 1, r_ind), points_matrix(3, 2, r_ind), 'ro')
hold off


%% start exercise
n = 128;                        % number of descriptors pairs

X_descriptor_points = randi(S, [2,n])-S/2; % each column [p1(1), p1(2)]'
Y_descriptor_points = randi(S, [2,n])-S/2; % each column [p2(1), p2(2)]'

% from [-S/2, +S/2]
figure; hold on;
for i = 1:n
    plot( [X_descriptor_points(1,i);Y_descriptor_points(1,i)], [X_descriptor_points(2,i);Y_descriptor_points(2,i)], '-k')
end
%%

% First filter the image
for i = 1:6
    Im = imread(strcat("wall/img", int2str(i), ".ppm"));
    Ims{i} = rgb2gray(Im);
    Ims{i} = imgaussfilt(Ims{i}, 2);
end

% Obtain the desciptors of each valid point in each image
Im1_descriptor_points = zeros(6, N_points, n);    % (Im, point, description)
for i = 1:6
    for j = 1:N_points
        if(valid_points_mat(i,j))
            c_point = round(squeeze(points_matrix(i,:,j)));
            idx1 = sub2ind(size(Ims{i}), c_point(2)+X_descriptor_points(2,:), c_point(1)+X_descriptor_points(1,:));
            idx2 = sub2ind(size(Ims{i}), c_point(2)+Y_descriptor_points(2,:), c_point(1)+Y_descriptor_points(1,:));
            Im1_descriptor_points(i,j,:) = Ims{i}(idx1) > Ims{i}(idx2);
        else
            Im1_descriptor_points(i,j,:) = -1;
        end
    end
end
% descriptor comparison
predictions_correct = zeros(6,N_points);             % save in a matrix if it hited right
predictions_distance = zeros(6,N_points);           % save in a matrix the distances
for j = 1:N_points % for all points in the first image
    xor_dist = sum(xor(Im1_descriptor_points(1,j,:), Im1_descriptor_points( ...
        logical(valid_points_mat(:,j)),:,:)), 3);       % xor=> matrix of size (valid_points, 1, n) thend add third dimension
    [m,m_ind] = min(xor_dist ,[],2 );                   % find the minimum
    predictions_correct(logical(valid_points_mat(:,j)),j) = (m_ind == j); % did it hit?
    idx1 = sub2ind(size(xor_dist), [1:size(xor_dist, 1)], m_ind');        % indices of the descriptor
    predictions_distance(logical(valid_points_mat(:,j)),j) = xor_dist(idx1); % compute distance
end

%% results

which_image = 3;

True_pos_dist = predictions_distance(which_image,and(logical(predictions_correct(which_image,:)), logical(valid_points_mat(which_image,:)))   ) ;
False_pos_dist  = predictions_distance(which_image,and(~ logical(predictions_correct(which_image,:)), logical(valid_points_mat(which_image,:)))   ) ;


figure;
h1 = histogram(True_pos_dist,40,'Normalization','pdf')
hold on
h2 = histogram(False_pos_dist,40,'Normalization','pdf')
legend('True positive', 'False positive')
hold off

correct_hits = sum( and(logical(predictions_correct(which_image,:)), logical(valid_points_mat(which_image,:)))  )
bad_hits = sum(  and(~logical(predictions_correct(which_image,:)), logical(valid_points_mat(which_image,:)))  )

TPR = hist(True_pos_dist , 40, 'Normalization', 'pdf');
FPR = hist(False_pos_dist , 40, 'Normalization', 'pdf');

hit_rate = correct_hits/(correct_hits+bad_hits)

figure()
title('Roc curve')
plot(cumsum(FPR), cumsum(TPR))






%% new descriptor proposed

theta = linspace(0,2*pi, 100);
circle_matrix = zeros(33,33);

r = S/2;
circle_points = [r*cos(theta);r*sin(theta)];            % radius S/2 outer circle
circle_points = unique( round(circle_points)' , 'rows' )';% round and eliminate the repeated points

circle_matrix( sub2ind([33,33] , circle_points(1,:)+S/2+1, circle_points(2,:)+S/2+1 ) ) = 1;
n_left = n-round(n/3);
X_descriptor_points_new = circle_points(:,randi(size(circle_points,2),round(n/3),1));
Y_descriptor_points_new = circle_points(:,randi(size(circle_points,2),round(n/3),1)); % one third to the outer circle

% midle circle
r = S/3;
circle_points = [r*cos(theta);r*sin(theta)];            % radius S/3 outer circle
circle_points = unique( round(circle_points)' , 'rows' )';% round and eliminate the repeated points

circle_matrix( sub2ind([33,33] , circle_points(1,:)+S/2+1 , circle_points(2,:)+S/2+1) ) = 1;
n_left = n_left - round(n/2);
X_descriptor_points_new = [X_descriptor_points_new, circle_points(:,randi(size(circle_points,2),round(n/2),1))];
Y_descriptor_points_new = [Y_descriptor_points_new, circle_points(:,randi(size(circle_points,2),round(n/2),1))]; % one half to the middle circle

% inner circle
r = S/6;
circle_points = [r*cos(theta);r*sin(theta)];            % radius S/3 outer circle
circle_points = unique( round(circle_points)' , 'rows' )';% round and eliminate the repeated points

circle_matrix( sub2ind([33,33] , circle_points(1,:)+S/2+1 , circle_points(2,:)+S/2+1 ) ) = 1;
X_descriptor_points_new = [ X_descriptor_points_new, circle_points(:,randi(size(circle_points,2),n_left,1) ) ];
Y_descriptor_points_new = [ Y_descriptor_points_new, circle_points(:,randi(size(circle_points,2),n_left,1) ) ]; % one half to the middle circle


figure
imshow(circle_matrix)
hold on
for i = 1:n
    plot( [X_descriptor_points_new(1,i);Y_descriptor_points_new(1,i)]+S/2+1, [X_descriptor_points_new(2,i);Y_descriptor_points_new(2,i)]+S/2+1, '-r')
end
%%

% First filter the image
for i = 1:6
    Im = imread(strcat("wall/img", int2str(i), ".ppm"));
    Ims{i} = rgb2gray(Im);
    Ims{i} = imgaussfilt(Ims{i}, 2);
end

% Obtain the desciptors of each valid point in each image
Im1_descriptor_points = zeros(6, N_points, n);    % (Im, point, description)
for i = 1:6
    for j = 1:N_points
        if(valid_points_mat(i,j))
            c_point = round(squeeze(points_matrix(i,:,j)));
            idx1 = sub2ind(size(Ims{i}), c_point(2)+X_descriptor_points_new(2,:), c_point(1)+X_descriptor_points_new(1,:));
            idx2 = sub2ind(size(Ims{i}), c_point(2)+Y_descriptor_points_new(2,:), c_point(1)+Y_descriptor_points_new(1,:));
            Im1_descriptor_points(i,j,:) = Ims{i}(idx1) > Ims{i}(idx2);
        else
            Im1_descriptor_points(i,j,:) = -1;
        end
    end
end
% descriptor comparison
predictions_correct = zeros(6,N_points);             % save in a matrix if it hited right
predictions_distance = zeros(6,N_points);           % save in a matrix the distances
for j = 1:N_points % for all points in the first image
    xor_dist = sum(xor(Im1_descriptor_points(1,j,:), Im1_descriptor_points( ...
        logical(valid_points_mat(:,j)),:,:)), 3);       % xor=> matrix of size (valid_points, 1, n) thend add third dimension
    [m,m_ind] = min(xor_dist ,[],2 );                   % find the minimum
    predictions_correct(logical(valid_points_mat(:,j)),j) = (m_ind == j); % did it hit?
    idx1 = sub2ind(size(xor_dist), [1:size(xor_dist, 1)], m_ind');        % indices of the descriptor
    predictions_distance(logical(valid_points_mat(:,j)),j) = xor_dist(idx1); % compute distance
end

%% results

which_image = 6;

True_pos_dist = predictions_distance(which_image,and(logical(predictions_correct(which_image,:)), logical(valid_points_mat(which_image,:)))   ) ;
False_pos_dist  = predictions_distance(which_image,and(~ logical(predictions_correct(which_image,:)), logical(valid_points_mat(which_image,:)))   ) ;


% figure;
% h1 = histogram(True_pos_dist,40,'Normalization','pdf')
% hold on
% h2 = histogram(False_pos_dist,40,'Normalization','pdf')
% legend('True positive', 'False positive')
% hold off

correct_hits = sum( and(logical(predictions_correct(which_image,:)), logical(valid_points_mat(which_image,:)))  )
bad_hits = sum(  and(~logical(predictions_correct(which_image,:)), logical(valid_points_mat(which_image,:)))  )

TPR = hist(True_pos_dist , 33, 'Normalization', 'pdf');
FPR = hist(False_pos_dist , 33, 'Normalization', 'pdf');

hit_rate = correct_hits/(correct_hits+bad_hits)

figure()
title('Roc curve')
plot(cumsum(FPR), cumsum(TPR))


%% new descriptor

circle_matrix = zeros(33,33);
theta = linspace(0,2*pi,100);

% outer circle
r = S/2;
circle_points = [r*cos(theta);r*sin(theta)];            % radius S/2 outer circle
circle_points = unique( round(circle_points)' , 'rows' )';% round and eliminate the repeated points
circle_matrix( sub2ind([33,33] , X_descriptor_points_new(1,:)+S/2+1 , X_descriptor_points_new(1,:)+S/2+1) ) = 1;

n_left = n-round(n/3);
X_descriptor_points_new = circle_points(:,randi(size(circle_points,2),round(n/3)));% one third to the outer circle
y_descriptor_points_new = circle_points(:,randi(size(circle_points,2),round(n/3)));

% middle circle
r = S/3;
circle_points = [r*cos(theta);r*sin(theta)];            % radius S/2 outer circle
circle_points = unique( round(circle_points)' , 'rows' )';% round and eliminate the repeated points
circle_matrix( sub2ind([33,33] , X_descriptor_points_new(1,:)+S/2+1 , X_descriptor_points_new(1,:)+S/2+1) ) = 1;

n_left = n_left-round(n/2);
X_descriptor_points_new = circle_points(:,randi(size(circle_points,2),round(n/2)));% one third to the outer circle
y_descriptor_points_new = circle_points(:,randi(size(circle_points,2),round(n/2)));
% inner circle
r = S/6;
circle_points = [r*cos(theta);r*sin(theta)];            % radius S/2 outer circle
circle_points = unique( round(circle_points)' , 'rows' )';% round and eliminate the repeated points
circle_matrix( sub2ind([33,33] , X_descriptor_points_new(1,:)+S/2+1 , X_descriptor_points_new(1,:)+S/2+1) ) = 1;

n_left = n-round(n/3);
X_descriptor_points_new = circle_points(:,randi(size(circle_points,2),n_left));% one third to the outer circle
y_descriptor_points_new = circle_points(:,randi(size(circle_points,2),n_left));

figure






