function [this_mean] = RGB_mean(X)
% we assume that X is of shape [n,3];
this_mean = mean( X, 1 );   % choose the Euclidean mean function
end