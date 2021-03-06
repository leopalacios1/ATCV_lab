function [this_mean] = HSV_mean(X)
% we assume that X is of shape [n,3];
mean_cos = mean(  cos(X(:,1)*2*pi)  );
mean_sin = mean(  sin(X(:,1)*2*pi)  );
mean_angle = atan2(mean_sin, mean_cos)/(2*pi);
if(mean_angle < 0)
    mean_angle = mean_angle + 1;
end
% the Hue mean is given by the sinus and the cosinus (this could slow
% computations).
this_mean = [ mean_angle, mean(X(:,2)), mean(X(:,3))  ]; 

end