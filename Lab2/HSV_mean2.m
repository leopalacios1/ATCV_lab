function [this_mean] = HSL_mean(X)
% we assume that X is of shape [n,3];
mean_cos = mean( cos(X(:,1)*2*pi) );
mean_sin = mean( sin(X(:,1)*2*pi) );

% the Hue mean is given by the sinus and the cosinus (this could slow
% computations).
this_mean = [atan2(mean_sin, mean_cos)./(2*pi), mean(X(:,2)), mean(X(:,3))];
Idx       = this_mean(:,1)<0;
this_mean(Idx,1) = 1 + this_mean(Idx,1);

end