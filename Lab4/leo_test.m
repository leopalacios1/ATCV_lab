
%% BASIC STEPS FOR FILTERING USING DFT 
% 1.Obtain the padding parameters:
PQ = paddedsize(size(f));
% 2. Obtain the Fourier transform.
F = fft2(f, PQ(1), PQ(2));
% (Spectrum:)
S = abs(F);
% 3. Generate a filter function, H, of size PQ(1)xPQ(2). The filter must be in the format shown in fig.4 (right).
% If it is centered instead, as in fig.4 (left) let H = fftshift(H) before usign the filter.
% 4. Multiply the trasform by the filter. 
G = H.* F;
% 5. Obtain the real part of the inverse FFT of G: 
g = real(ifft2(G));
% 6. Crop the top, left rectangle to the original size: 
g = g(1:size(f, 1), 1:size(f, 2));

%% Computing and visualizing the 2-D DFT

f = imread('truck.jpg');
f = rgb2gray(f);

PQ = paddedsize(size(f));

%F = fft2(f);
F= fft2(f, PQ(2), PQ(1));
S= abs(F);

Fc = fftshift(F); 
imshow(Fc, [ ])

whos


%% Exercise r 
N = 9; % Size  in one dimension of the input photo
M = 5; % Size  in the other dimension of the output photo
[U,V] = dftuv(N,M);
D = U.^2 + V.^2;
fftshift(D)


%% 4-4 

% function H = lpfilter(tipus, M, N, D0, n)
% LPFILTER computes frequency domain lowpass filters % H = LPFILTER(TYPE, M,N,D0,N) creates the transfer
% function of a lowpass filter, H, of the specified TYPE and % size (M-by-N). To view the filter as an image or mesh
% plot, it should be centered using H=fftshift(H).
%
% Valid values for TYPE, DO and n are:
% ...
[U,V] = dftuv(N,M);
D = U.^2 + V.^2;

switch tipus
    case ideal
        H = double(D<=D0)

