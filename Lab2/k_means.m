function [indx, centers] = k_means(X, k, dist_funct)
% we assume X is [n, dims] size
% set tolerance and maximum number of iterations;
tol = 1e-4; max_it = 20;
% initialize the indexes [ size(X,1) ]
indx = zeros(size(X,1), 1);    

% choose random k centers from the input [k, size(X,2)]
centers = X( randi(size(X,1), [k,1]) ,:);
prev_centers = zeros(k, size(X,2));

% allocate the distance matrix [n,k]
dist_matrix = zeros( size(X,1), k)

curr_it = 1;

% start the loop
while(curr_it <= max_it && max(dist_funct(centers, prev_centers)) > tol )
    
    prev_centers = centers;

    for i = 1:k         % define the distances
        dist_matrix(:,i) = dist_funct( centers(i,:), X );
    end

    [m,indx] = min(dist_matrix, [], 2); % the index goes from [1,k]

    for i = 1:k         % update the centers
        centers = mean( X( indx==i,:), 2 );   % choose the Euclidean mean function
    end

    curr_it = curr_it + 1;
end
end
