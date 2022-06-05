function [clean_im, source_fft, clean_fft] = filter_pattern_noise(source_im, D0, th)

if(nargin == 1)
    D0 = 20;
    th = 0.7;
else if(nargin == 2)
    th = 0.7;
end
PQ = paddedsize(size(source_im));
source_fft = fft2(source_im, PQ(1), PQ(2));

H = hpfilter("ideal", PQ(1),PQ(2),D0);
noise_fft = H.*source_fft;

Fc_norm = abs(noise_fft);
S_patt = log(1+Fc_norm)/log( max(Fc_norm, [], 'all') );
m = max(S_patt, [], 'all');

clean_fft = source_fft;
clean_fft(S_patt>th*m) = 0;


clean_im = ifft2(clean_fft);
clean_im = clean_im(1:size(source_im, 1), 1:size(source_im, 2));
end