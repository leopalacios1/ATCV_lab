close all;
clear all;
clc

im = im2double( imread("moon.jpg") );
imshow(im); title( "original image" ); 
PQ = paddedsize( size(im) );

F = fft2(im, PQ(1), PQ(2));
%%
S = abs(F);
imshow(S, []), title("Spectrum")
%%
S2 = log(1+S)/log(max(S, [], 'all'));
imshow(S2)

%%
