function [HSL] = rgb2HSL(Im)
% assuming that the image is [n, rgb];
    HSL = zeros(size(Im));
    % H = atan2( (B-G)/sqrt(2)  ,  (2R-B-G)/sqrt(6) )
    HSL(:,1) = atan2( (-Im(:,3) + Im(:,2))/sqrt(2), (2*Im(:,1) - Im(:,3) - Im(:,2))/sqrt(6)   );
    M = max(Im, [], 2)/255.;
    m = min(Im, [], 2)/255.;
    % L = M/(1+M-n)
    HSL(:,3) = M./(1+M-m);
    HSL(:,2) = 2*(M-m)./(  1+abs(M-0.5) + abs(m-0.5)  );
end