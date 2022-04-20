function [dist] = cityblock_dist(X1,X2)
% assume X1 and X2 are [n,dim] or one of them is [1,dim]
dist = sum((abs(X1 - X2)),2);

end