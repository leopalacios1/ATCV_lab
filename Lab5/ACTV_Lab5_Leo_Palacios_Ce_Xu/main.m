%% Adding the watermark

im = imread("photo2.jpg");
wm = imread("water_mark.jpg");
wm = rgb2gray(wm);
wm = double(wm);

PQ = paddedsize( size(im, 1:2), size(wm) , "PWR2");

yuv_im = rgb2ycbcr(im);
y_im = double(yuv_im(:,:,1));
u_im = yuv_im(:,:,2);
v_im = yuv_im(:,:,3);

y_im_augmented = zeros(PQ);
y_im_augmented(1:size(y_im,1), 1:size(y_im,2)) = y_im;

wm_augmented = zeros(PQ);
start_index = int32( (size(im, [1,2])-size(wm))/2 );
h_ind = start_index(1);
w_ind = start_index(2);
wm_augmented( h_ind:h_ind+size(wm, 1)-1, w_ind:w_ind+size(wm, 2)-1 ) = wm;

hadamard_y_im = fwht(y_im_augmented, PQ(1), 'hadamard');
hadamard_w_wm = fwht(wm_augmented, PQ(1), 'hadamard');

%%

m = 1;  % controlling parameter if 0 => image destroyed, 1 => good wm, 2 => not visible wm
alpha = 1/pi*( atan(mean(  hadamard_y_im , 'all')) + pi/2 )/(10^m)

% alpha = 0.15;

hadamard_y2_im = hadamard_y_im + alpha*hadamard_w_wm;

y_im2 = ifwht(hadamard_y2_im, PQ(1), 'hadamard');
y_im2 = y_im2(1:size(y_im,1), 1:size(y_im,2));
yuv_im2 = cat(3,  y_im2 ,u_im,v_im);
im2 = ycbcr2rgb( yuv_im2 );

figure
imshow(im2)

%% Removing the watermark

% here we are using only the previous watermark matrices and the 
% watermarked image im2. alpha could be known or computed according to m

yuv_im2 = rgb2ycbcr(im2);
y_im2 = yuv_im2(:,:,1);
u_im = yuv_im2(:,:,2);
v_im = yuv_im2(:,:,3);

y_im2_augmented = zeros(PQ);
y_im2_augmented(1:size(y_im2,1), 1:size(y_im2,2)) = y_im2;

hadamard_y_im2 = fwht(y_im2_augmented, PQ(1), 'hadamard');

% alpha = 1/pi*( atan(mean(  hadamard_y_im2 , 'all')) + pi/2 )/(10^m)

hadamard_y_im3 = hadamard_y_im2 - alpha*hadamard_w_wm;

y_im3 = ifwht(hadamard_y_im3, PQ(1), 'hadamard');
y_im3 = y_im3(1:size(y_im2,1), 1:size(y_im2,2));
yuv_im3 = cat(3,y_im3   ,u_im,v_im);
im3 = ycbcr2rgb( yuv_im3 );

% figure
% imshow(im3)



figure

subplot( 1,3,1 )
imshow(im)
title("original")

subplot( 1,3,2 )
imshow(im2)
title("watermarked")

subplot( 1,3,3 )
imshow(im3)
title("recovered")

%% Extracting watermark

% in this section, we are going to assume that we have the original source
% image and the watermarked one, so we have to extract the logo


yuv_im = rgb2ycbcr(im);
y_im = double(yuv_im(:,:,1));
u_im = yuv_im(:,:,2);
v_im = yuv_im(:,:,3);
y_im_augmented = zeros(PQ);
y_im_augmented(1:size(y_im,1), 1:size(y_im,2)) = y_im;
hadamard_y_im = fwht(y_im_augmented, PQ(1), 'hadamard');

yuv_im2 = rgb2ycbcr(im2);
y_im2 = yuv_im2(:,:,1);
y_im2_augmented = zeros(PQ);
y_im2_augmented(1:size(y_im2,1), 1:size(y_im2,2)) = y_im2;
hadamard_y_im2 = fwht(y_im2_augmented, PQ(1), 'hadamard');

% alpha = 1/pi*( atan(mean(  hadamard_y_im , 'all')) + pi/2 )/(10^m)

hadamard_extracted_wm = (hadamard_y_im2-hadamard_y_im)/alpha;

wm_2 = ifwht(hadamard_extracted_wm, PQ(1), 'hadamard');
wm_2 = wm_2(1:size(y_im2,1), 1:size(y_im2,2));

figure
imshow(wm_2/255)



