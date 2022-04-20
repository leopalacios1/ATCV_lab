function [RGB] = HSL2rgb(X);
% we assume that X is [n,3] vector of HSL pixels
% with H in [-pi,pi], S in [0,1] and L in [0,1]
rgb = zeros(size(X));

% C = ( 1 - abs(2*L-1) ).*S
C = ( 1 - abs(2*X(:,3)-1) ).*X(:,2);
% x = C*( 1-abs( mod(3*H/pi, 2)  - 1) );
H = mod(X(:,1), 2*pi);
x = C.*( 1 - abs( mod(3*H/pi, 2) - 1 )  );
% m = L-C/2
m = X(:,3) - C/2;

% now for cases 
bool_vector = H < pi/3;
rgb(bool_vector,:) = [C(bool_vector), x(bool_vector), zeros(sum(bool_vector),1)];

bool_vector = and(H >= pi/3 , H < 2*pi/3);
rgb(bool_vector,:) = [x(bool_vector), C(bool_vector), zeros(sum(bool_vector),1)];

bool_vector = and( H >= 2*pi/3 , H < pi) ;
rgb(bool_vector,:) = [zeros(sum(bool_vector),1), C(bool_vector), x(bool_vector)];

bool_vector = and( H < 4*pi/3 , H >= pi) ;
rgb(bool_vector,:) = [zeros(sum(bool_vector),1), x(bool_vector), C(bool_vector)];

bool_vector = and( H < 5*pi/3 , H >= 4*pi/3) ;
rgb(bool_vector,:) = [x(bool_vector), zeros(sum(bool_vector),1), C(bool_vector)];

bool_vector = H >= 5*pi/3;
rgb(bool_vector,:) = [C(bool_vector), zeros(sum(bool_vector),1), x(bool_vector)];

% RGB = ( (r+m)*255, (g+m)*255, (b+m)*255 )
RGB = uint8( [ (rgb(:,1)+m)*255, (rgb(:,2)+m)*255, (rgb(:,3)+m)*255 ]);

end