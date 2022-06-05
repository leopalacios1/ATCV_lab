
f = imread("elephant.jpg");

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

%% Exercise 4.1 modify the TF

im1 = imread('truck.jpg');
im1 = rgb2gray(im1);

PQ = paddedsize(size(im1));

%F = fft2(f);
F= fft2(im1, PQ(1), PQ(2));
Fc = fftshift(F); 
FC_norm = abs(Fc);
S = log(1+FC_norm)/log(max(FC_norm, [], 'all'));
% transform

Fc2 = Fc;
% transformation_1
Fc2( int32(PQ(1)/2)-10:int32(PQ(1)/2)+10, int32(PQ(2)/2)-10:int32(PQ(2)/2)+10  ) ...
= Fc( int32(PQ(1)/2)-10:int32(PQ(1)/2)+10, int32(PQ(2)/2)-10:int32(PQ(2)/2)+10  ) * 1.4;
% / FC_norm( int32(PQ(1)/2)-10:int32(PQ(1)/2)+10, int32(PQ(2)/2)-10:int32(PQ(2)/2)+10  ) ;

% % transformation_2
% Fc2( int32(PQ(1)*2/3)-10:int32(PQ(1)*2/3)+10, int32(PQ(2)*2/3)-10:int32(PQ(2)*2/3)+10  )...
% = Fc( int32(PQ(1)*2/3)-10:int32(PQ(1)*2/3)+10, int32(PQ(2)*2/3)-10:int32(PQ(2)*2/3)+10  )...
% / FC_norm( int32(PQ(1)*2/3)-10:int32(PQ(1)*2/3)+10, int32(PQ(2)*2/3)-10:int32(PQ(2)*2/3)+10  ) * max(FC_norm,[], 'all')/500 ;
% 
% Fc2( int32(PQ(1)*1/3)-10:int32(PQ(1)*1/3)+10, int32(PQ(2)*1/3)-10:int32(PQ(2)*1/3)+10  )...
% = Fc( int32(PQ(1)*1/3)-10:int32(PQ(1)*1/3)+10, int32(PQ(2)*1/3)-10:int32(PQ(2)*1/3)+10  )...
% / FC_norm( int32(PQ(1)*1/3)-10:int32(PQ(1)*1/3)+10, int32(PQ(2)*1/3)-10:int32(PQ(2)*1/3)+10  ) * max(FC_norm,[], 'all')/500 ;
% 
% % transformation_3
% Fc2( int32(PQ(1)*2/3)-10:int32(PQ(1)*2/3)+10, int32(PQ(2)*1/3)-10:int32(PQ(2)*1/3)+10  )...
% = Fc( int32(PQ(1)*2/3)-10:int32(PQ(1)*2/3)+10, int32(PQ(2)*1/3)-10:int32(PQ(2)*1/3)+10  )...
% / FC_norm( int32(PQ(1)*2/3)-10:int32(PQ(1)*2/3)+10, int32(PQ(2)*1/3)-10:int32(PQ(2)*1/3)+10  ) * max(FC_norm,[], 'all')/500 ;
% 
% Fc2( int32(PQ(1)*1/3)-10:int32(PQ(1)*1/3)+10, int32(PQ(2)*2/3)-10:int32(PQ(2)*2/3)+10  )...
% = Fc( int32(PQ(1)*1/3)-10:int32(PQ(1)*1/3)+10, int32(PQ(2)*2/3)-10:int32(PQ(2)*2/3)+10  )...
% / FC_norm( int32(PQ(1)*1/3)-10:int32(PQ(1)*1/3)+10, int32(PQ(2)*2/3)-10:int32(PQ(2)*2/3)+10  ) * max(FC_norm,[], 'all')/500 ;


F2 = ifftshift(Fc2);
Fc2_norm = abs(Fc2);
S2 = log(1+Fc2_norm)/log(max(Fc2_norm, [], 'all'));
im2 = real(ifft2(F2)/255);
im2 = im2(1:size(im1, 1), 1:size(im1,2));

figure
subplot(2,2,1)
imshow(im1, [ ])
title("original image")
subplot(2,2,2)
imshow(S)
title(" Fc ")

subplot(2,2,3)
imshow(S2)
title( "transformed Fc" )

subplot(2,2,4)
imshow(im2)
title(" new image ")

%% Ex 4.2 Check the effects in the inverse transform whether fftshift is executed over the spectrum or not

im1 = imread('truck.jpg');
im1 = rgb2gray(im1);

PQ = paddedsize(size(im1));

%F = fft2(f);
F= fft2(im1, PQ(1), PQ(2));
Fc = fftshift(F); 
FC_norm = abs(Fc);
S = log(1+FC_norm)/log(max(FC_norm, [], 'all'));


Fc2 = Fc;

F2 = ifftshift(Fc2);

im2 = real(ifft2(Fc2)/255);
im2 = im2(1:size(im1, 1), 1:size(im1,2));
im3 = real(ifft2(F2)/255);
im3 = im3(1:size(im1, 1), 1:size(im1,2));
figure
subplot(1,3,1)
imshow(im1)
title(" original image " )
subplot(1,3,2)
imshow(im3)
title(" without shift " )
subplot(1,3,3)
imshow(im2)
title(" shifted " )


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


%% Exercise 4.3
N = 9; % Size  in one dimension of the input photo
M = 5; % Size  in the other dimension of the output photo
[U,V] = dftuv(N,M);
D = U.^2 + V.^2
D2 = fftshift(D)

% the coordinates of the center are the half of each dim

%% 4-4 


im = imread('truck.jpg');
im = rgb2gray(im);

M = size(im,1);
N = size(im,2);

figure
subplot(1,3,1)
H = lpfilter("ideal", M,N,100);
imshow(H)
title("ideal")
subplot(1,3,2)
H = lpfilter("btw", M,N,100);
imshow(H)
title("btw")

subplot(1,3,3)
H = lpfilter("gauss", M,N,100);
imshow(H)
title("gauss")
%%
figure
h=lpfilter('ideal',128,128,18, 4);
surf(fftshift(h))
colormap(gray)
grid off
axis off
shading interp

figure
h=lpfilter('btw',128,128,18, 4);
surf(fftshift(h))
colormap(gray)
grid off
axis off
shading interp

figure
h=lpfilter('gauss',128,128,18, 4);
surf(fftshift(h))
colormap(gray)
grid off
axis off
shading interp
%% 4.5

im = imread("moon.jpg");
PQ = paddedsize(size(im));

%F = fft2(f);
F= fft2(im, PQ(1), PQ(2));
D0 = 60;
M = PQ(1); N = PQ(2);

figure
subplot(2,2,1)
H = lpfilter("ideal", M,N,D0);
filtered_F = H.*F;
im2 = real(ifft2(filtered_F)/255);
im2 = im2(1:size(im,1), 1:size(im,2)  );
imshow(im2)
title("ideal")

subplot(2,2,2)
H = lpfilter("btw", M,N,D0, 2);
filtered_F = H.*F;
im2 = real(ifft2(filtered_F)/255);
im2 = im2(1:size(im,1), 1:size(im,2)  );
imshow(im2)
title("btw, 2")

subplot(2,2,3)
H = lpfilter("btw", M,N,D0, 4);
filtered_F = H.*F;
im2 = real(ifft2(filtered_F)/255);
im2 = im2(1:size(im,1), 1:size(im,2)  );
imshow(im2)
title("btw, 4")

subplot(2,2,4)
H = lpfilter("gauss", M,N,D0, 4);
filtered_F = H.*F;
im2 = real(ifft2(filtered_F)/255);
im2 = im2(1:size(im,1), 1:size(im,2)  );
imshow(im2)
title("gauss")

%% Ex 4.6

im = imread('truck.jpg');
im = rgb2gray(im);

M = size(im,1);
N = size(im,2);

figure
subplot(1,3,1)
H = hpfilter("ideal", M,N,100);
imshow(H)
title("ideal")
subplot(1,3,2)
H = hpfilter("btw", M,N,100);
imshow(H)
title("btw")

subplot(1,3,3)
H = hpfilter("gauss", M,N,100);
imshow(H)
title("gauss")

figure
h=hpfilter('ideal',128,128,18, 4);
surf(fftshift(h))
colormap(gray)
grid off
axis off
shading interp

figure
h=hpfilter('btw',128,128,18, 4);
surf(fftshift(h))
colormap(gray)
grid off
axis off
shading interp

figure
h=hpfilter('gauss',128,128,18, 4);
surf(fftshift(h))
colormap(gray)
grid off
axis off
shading interp

%% Ex 4.7


im = imread("moon.jpg");
PQ = paddedsize(size(im));


%F = fft2(f);
F= fft2(im, PQ(1), PQ(2));
D0 = 20;
a = 0.8; b = .4;

M = PQ(1); N = PQ(2);

figure
subplot(2,2,1)
H = hpfilter("ideal", M,N,D0);
H = a + b*H;
filtered_F = H.*F;
im2 = real(ifft2(filtered_F)/255);
im2 = im2(1:size(im,1), 1:size(im,2)  );
imshow(im2)
title("ideal")

subplot(2,2,2)
H = hpfilter("btw", M,N,D0, 2);
H = a + b*H;
filtered_F = H.*F;
im2 = real(ifft2(filtered_F)/255);
im2 = im2(1:size(im,1), 1:size(im,2)  );
imshow(im2)
title("btw, 2")

subplot(2,2,3)
H = hpfilter("btw", M,N,D0, 4);
H = a + b*H;
filtered_F = H.*F;
im2 = real(ifft2(filtered_F)/255);
im2 = im2(1:size(im,1), 1:size(im,2)  );
imshow(im2)
title("btw, 4")

subplot(2,2,4)
H = hpfilter("gauss", M,N,D0, 4);
H = a + b*H;
filtered_F = H.*F;
im2 = real(ifft2(filtered_F)/255);
im2 = im2(1:size(im,1), 1:size(im,2)  );
imshow(im2)
title("gauss")

%% P1

im = imread("moon.jpg");
im = double(im)/255.;
[h,w] = size(im);
x_pattern = 0.5 + 0.5*sin([1:h]*2*pi/30);
pattern1 = ones(h,w).*x_pattern';
pattern = pattern1;
% pattern2 = zeros(h,w);
% for i = 1:h
%     for j = 1:w
%         pattern2(i,j) = 0.5 + 0.5*sin( (3*i+2*j)*2*pi/45 );
%     end
% end
% pattern = 0.5*pattern1 + 0.5*pattern2;
patt_im = 0.7*im + 0.3*pattern;
[clean_im, source_fft, clean_fft] = filter_pattern_noise(patt_im, 20);

PQ = paddedsize(size(patt_im));
F_patt = fft2(pattern, PQ(1), PQ(2));


Fc_patt = fftshift(F_patt); 
Fc_norm = abs(Fc_patt);
S_patt = log(1+abs(Fc_norm))/log(max(abs(Fc_norm), [], 'all'));


Fc_im = fftshift(source_fft); 
S_im = log(1+abs(Fc_im))/log(max(abs(Fc_im), [], 'all'));

Fc_filtered = fftshift(clean_fft); 
S_filtered = log(1+abs(Fc_filtered))/log(max(abs(Fc_filtered), [], 'all'));

figure
subplot(1,3,1)
imshow(S_patt)
title("pattern fft")

subplot(1,3,2)
imshow(S_im)
title("image fft")

subplot(1,3,3)
imshow(S_filtered)
title("filtered fft")



figure
subplot(1,3,1)
imshow(pattern)
title("pattern")

subplot(1,3,2)
imshow(patt_im)
title("image")

subplot(1,3,3)
imshow(clean_im)
title("filtered image")

%% P2 Information
clear;clc;

im1 = imread("truck.jpg");
im2 = imread("elephant.jpg");

im1_grey = rgb2gray(im1);
im2_grey = rgb2gray(im2);
PQ = paddedsize(size(im1_grey));

F_1 = fft2(im1_grey, PQ(1), PQ(2));
F_2 = fft2(im2_grey, PQ(1), PQ(2));

F1_abs = abs(F_1);
S_1 = log(1+abs(F1_abs))/log(max(abs(F1_abs), [], 'all'));
S_1 = fftshift(S_1);
F2_abs = abs(F_2);
S_2 = log(1+abs(F2_abs))/log(max(abs(F2_abs), [], 'all'));
S_2 = fftshift(S_2);

F1_angle = angle(F_1);
F2_angle = angle(F_2);

F1_new = F1_abs.*exp(i*F2_angle);
F2_new = F2_abs.*exp(i*F1_angle);

im1_new = real(ifft2(F1_new)/255);
im1_new = im1_new(1:size(im1, 1), 1:size(im1,2)  );
im2_new = real(ifft2(F2_new)/255);
im2_new = im2_new(1:size(im2, 1), 1:size(im2,2)  );

figure
subplot(3,2,1)
imshow( im1 )
subplot(3,2,2)
imshow( im2 )
subplot(3,2,3)
imshow( S_1 )
title(" log magnitude of the truck ")
subplot(3,2,4)
imshow( S_2 )
title(" log magnitude of the elephant ")
subplot(3,2,5)
imshow( im1_new )
title("Magnitude of truck, phase of elephant")
subplot(3,2,6)
imshow( im2_new )
title("Magnitude of elephant, phase of truck")

%%

%% P2 Information
clear;clc;

im1 = imread("truck.jpg");
im2 = imread("elephant.jpg");

im1_grey = rgb2gray(im1);
im2_grey = rgb2gray(im2);
PQ = paddedsize(size(im1_grey));

F_1 = fft2(im1, PQ(1), PQ(2));
F_2 = fft2(im2, PQ(1), PQ(2));

F1_abs = abs(F_1);
S_1 = log(1+abs(F1_abs))/log(max(abs(F1_abs), [], 'all'));
S_1 = fftshift(S_1);
F2_abs = abs(F_2);
S_2 = log(1+abs(F2_abs))/log(max(abs(F2_abs), [], 'all'));
S_2 = fftshift(S_2);

F1_angle = angle(F_1);
F2_angle = angle(F_2);

F1_new = F1_abs.*exp(i*F2_angle);
F2_new = F2_abs.*exp(i*F1_angle);

im1_new = real(ifft2(F1_new)/255);
im1_new = im1_new(1:size(im1, 1), 1:size(im1,2),:  );
im2_new = real(ifft2(F2_new)/255);
im2_new = im2_new(1:size(im2, 1), 1:size(im2,2),:  );

figure
subplot(1,2,1)
imshow( im1_new )
title("Magnitude of truck, phase of elephant")
subplot(1,2,2)
imshow( im2_new )
title("Magnitude of elephant, phase of truck")






