function g = dftfilt(f,H)
%DFTFILT performs frequency domain filtering 
%   G = DFTFILT(F,H) filters F in the frequency domain using the
%   filter transfer function H. The output, G, is the filtered
%   image, which has the same size as F. DFTFILT automatically 
%   pads F to be the same size as H. 

%   DFTFILT assumes that F is real and that H is a real, uncentered,
%   circularly-symmetric filter function.

% Obtain the FFT of the padded input
F = fft2(f, size(H,1), size(H,2));

% Perform filtering 
g = real(ifft2(H.*F));

% Crop to original size
g = g(1:size(f,1), 1:size(f,2));
