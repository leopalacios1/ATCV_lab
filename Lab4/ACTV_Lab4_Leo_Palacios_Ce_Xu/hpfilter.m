function H = hpfilter(tipus, M, N, D0, n)
% HPFILTER computes frequency domain highpass filters
% H = HPFILTER(TYPE, M,N,D0,N) creates the transfer
% function of a highpass filter, H, of the specified TYPE and
% size (M-by-N). To view the filter as an image or mesh
% plot, it should be centered using H=fftshift(H).
%
% Valid values for TYPE, DO and n are:
%
% 'ideal' ideal highpass filter with cutoff frequency D0. n
% need not be supplied. D0 must be positive.
%
% 'btw' Butterworth highpass filter of order n, and cutoff
% D0. The default value for n is 1.0. D0 must be
% positive.
%
% 'gauss' gaussian highpass filter with cutoff (standard
% deviation) D0. n need not be supplied. D0 must
% be positive.
if( nargin == 4)
    n = 1;
end
h = lpfilter(tipus, M, N, D0, n);
H = 1-h;