function [dist] = lab_dist(X1, X2)
% assume X1 and X2 are [n,dim] or one of them is [1,dim]
dist = vecnorm( X1-X2, 2, 2 ); % 2-norm, second dimension
end