function [dist] = HSV_dist(X1, X2)


dh = min(    abs(360*X1(1)-360*X2(:,1)),360-abs(360*X1(1)-360.*X2(:,1)))./ 180;
ds = abs(X1(2)-X2(:,2));
dv = abs(X1(3)-X2(:,3));

dist = sqrt(dh.^2+ds.^2+dv.^2);

end