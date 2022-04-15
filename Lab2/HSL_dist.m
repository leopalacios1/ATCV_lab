function [dist] = HSL_dist(X1, X2)
% assuming X1 and X2 are both [n,dim] or one is [1,dim].
% HSL dist has H in [-pi, pi], S in [0,1], L in [0,1].

% sqrt(S1*S2)*sin(  (H1-H2)/2  )^2
aux1 = sqrt( X1(:,2).*X2(:,2) ).*sin( (X1(:,1) - X2(:,1))/2 ).^2;
% sqrt(  (1-S1)(1-S2)  )*(L1-L2)^2
aux2 = sqrt( (1-X1(:,2)).*(1-X2(:,2)) ).*(X1(:,3) - X2(:,3)).^2;
% (S1-S2)^2
aux3 = (X1(:,2) - X2(:,2)).^2;
dist = sqrt( aux1 + aux2 + aux3 );

end