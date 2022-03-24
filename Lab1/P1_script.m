clear; clc; 

S = 10;                     % patch size

Ims = {};           % library for imaes
for i = 1:6
    Im = imread(strcat("wall/img", int2str(i), ".ppm"));
    Ims{i} = imgaussfilt(rgb2gray(Im), S/2);
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
% valid points are those that have some distance from 

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
% figure()
% imshow(Ims{1})
% hold on
% plot( squeeze(points_matrix(1,1,:)), squeeze(points_matrix(1,2,:)), '+g' )
% plot(points_matrix(1, 1, r_ind), points_matrix(1, 2, r_ind), 'ro')
% hold off
% figure()
% imshow(Ims{3})
% hold on
% plot( squeeze(points_matrix(3,1,:)), squeeze(points_matrix(3,2,:)), '+g' )
% plot(points_matrix(3, 1, r_ind), points_matrix(3, 2, r_ind), 'ro')
% hold off


%% start exercise
n = 128;                        % number of descriptors pairs
descriptor_points = randi(S, [4,n])-S/2; % each column [p1(1), p1(2), p2(1), p2(2)]'
% from [-S/2, +S/2]

%%
Im1_descriptor_points = zeros(6, N_points, n);    % (Im, point, description)
for i = 1:6
    for j = 1:N_points
        if(valid_points_mat(i,j))
            c_point = round(squeeze(points_matrix(i,:,j)));
            idx1 = sub2ind(size(Ims{i}), c_point(2)+descriptor_points(2,:), c_point(1)+descriptor_points(1,:));
            idx2 = sub2ind(size(Ims{i}), c_point(2)+descriptor_points(4,:), c_point(1)+descriptor_points(3,:));
            Im1_descriptor_points(i,j,:) = Ims{i}(idx1) > Ims{i}(idx2);
        end
    end
end
%% descriptor comparison
predictions_correct = ones(6,N_points);             % save in a matrix if it hited right
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

y = predictions_distance(which_image,and(logical(predictions_correct(which_image,:)), logical(valid_points_mat(which_image,:)))   ) ;
yy = predictions_distance(which_image,and(~ logical(predictions_correct(which_image,:)), logical(valid_points_mat(which_image,:)))   ) ;


figure()
hist(y, 20);
title('correct')

figure()
hist(yy,20);
title('non correct')

correct_hits = sum( and(logical(predictions_correct(which_image,:)), logical(valid_points_mat(which_image,:)))  )
bad_hits = sum(  and(~logical(predictions_correct(which_image,:)), logical(valid_points_mat(which_image,:)))  )