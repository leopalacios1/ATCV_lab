function H = lpfilter(tipus, M, N, D0, n)
% LPFILTER computes frequency domain lowpass filters
% H = LPFILTER(TYPE, M,N,D0,N) creates the transfer
% function of a lowpass filter, H, of the specified TYPE and
% size (M-by-N). To view the filter as an image or mesh
% plot, it should be centered using H=fftshift(H).
%
% Valid values for TYPE, DO and n are:
%
% 'ideal' ideal lowpass filter with cutoff frequency D0. n
% need not be supplied. D0 must be positive.
%
% 'btw' Butterworth lowpass filter of order n, and cutoff
% D0. The default value for n is 1.0. D0 must be
% positive.
%
% 'gauss' gaussian lowpass filter with cutoff (standard
% deviation) D0. n need not be supplied. D0 must
% be positive.

[U,V] = dftuv(M,N);

D = sqrt(U.^2 + V.^2);

switch tipus
    case 'ideal'
        H = double(D <= D0);
    case 'btw'
        if( nargin == 4):
            n = 1;
        end
        H = 1./( 1 + (D/D0).^(2^n));
    case 'gaussian'
        N = exp(-D^2./(2*D0^2));