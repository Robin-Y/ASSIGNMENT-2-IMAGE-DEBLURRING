clc;
close all;
clear all;
for i =1:4
%% original image
org = imread('C:\Users\Shweta\Documents\digital Image processing\assign2\Ground_truth_image.jpg');     % read original image
figure(); imshow(org); title('original image');
org_r = org(:,:,1);       % red component of given image
org_g = org(:,:,2);       % green component of given image
org_b = org(:,:,3);       % blue component of given image
[po,qo] = size(org);      % size of given image

%% blurred image from internet
path1 = 'C:\Users\Shweta\Documents\digital Image processing\assign2\Blurry1_';  % read original image
kernum1 = num2str(i);          % load kernel number
kerex1 = '.jpg';               % load kernel extension
ker1 = strcat(path1,kernum1,kerex1);     % load blur kernel
blrdimg = imread(ker1);       % read blcur kernel
figure(); imshow(blrdimg); title(['image blurred by kernel ',num2str(i),'i']);
blr_r = blrdimg(:,:,1);       % red component of given image

blr_g = blrdimg(:,:,2);       % green component of given image
blr_b = blrdimg(:,:,3);       % blue component of given image
[pb,qb] = size(blrdimg);      % size of given image

%% loading burring kernel
path = 'C:\Users\Shweta\Documents\digital Image processing\assign2\kernel\Kernel';          % load kernel path
kernum = num2str(i);          % load kernel number
kerex = 'G.png';               % load kernel extension
ker = strcat(path,kernum,kerex);     % load blur kernel
kernel1 = imread(ker);       % read blur kernel
kernel = rgb2gray(kernel1);             % convert blur kernel to grayscale
figure()
imshow(kernel)
[p,q] = size(kernel);                  % size of blur kernel
m = p+pb-1;              % size for result of fft calculation 
n = q+qb-1;
fftkernel = fft2(double(kernel),m,n);          % fft of blur kernel
figure()
imshow(abs(fftkernel),[])
title('FFT of kernel')
k = fftshift(fftkernel);
figure()
imshow(abs(k),[])
title('centerd fft of kernel')
figure()
imshow(log(1+abs(k)),[])
title('Log scaledand centered FFT of kernel')
IFF = ifft2(fftkernel);
FINAL_IM = uint8(real(IFF));  
Final_im = FINAL_IM(1:p,1:q);
figure
imshow(FINAL_IM)
title('IFFT of KERNEL')
figure()
imshow(Final_im)
title('IFFT of KERNEL after removing padding')
%fft of kernel red green and blue component
fftkernel_r = fft2(kernel1(:,:,1),m,n);
fftkernel_r = fftshift(fftkernel_r);
fftkernel_g = fft2(kernel1(:,:,2),m,n);
fftkernel_g = fftshift(fftkernel_g);
fftkernel_b = fft2(kernel1(:,:,3),m,n);
fftkernel_b = fftshift(fftkernel_b);



% FFT of blurred image red green and blue components
fftr = fft2(double(blr_r),m,n);
fft1r = fftshift(fftr);
figure()
imshow(mat2gray(log(1+abs(fft1r))))
title('FFT of red componenet of blurred image')
% le1 = ifft2(fft1r);
% figure()
% imshow(mat2gray(abs(le)))

fftg = fft2(double(blr_g),m,n);
fft2g = fftshift(fftg);
figure()
imshow(mat2gray(log(1+abs(fft2g))))
title('FFT of green componenet of blurred image')

fftb = fft2(double(blr_b),m,n);
fft3b = fftshift(fftb);
figure()
imshow(mat2gray(log(1+abs(fft3b))))
title('FFT of blue componenet of blurred image')
pause
close all
%% Deblurring using Inverse filter

fdeb_r = fft1r./fftkernel_r ;
deb_r = ifft2(ifftshift(fdeb_r));
fdeb_g = fft2g./fftkernel_g ;
deb_g = ifft2(ifftshift(fdeb_g));
fdeb_b = fft3b./fftkernel_b ;
deb_b = ifft2(ifftshift(fdeb_b));
inv_result = cat(3,real(deb_r),real(deb_g),real(deb_b));
figure()
imshow(mat2gray((inv_result(1:800,1:800,:))))
title('Result of Inverse filter')
figure()
imshow(mat2gray(log(1+(inv_result(1:800,1:800,:)))))
title('Result of Inverse filter after log transform')

pause
%% Deblurring using truncated filter
[fftkernel_rows, fftkernel_cols] = size(k);

a = (fftkernel_rows/2);
display(a*2)
b =  (fftkernel_cols/2);
display(b*2)
radius = input('Enter raduis of truncation filter')
p = ones(size(k));
for x = 1:fftkernel_rows
    for y = 1:fftkernel_cols
        if ((x-a)^2 + (y-b)^2<=radius^2)
            p(x,y) = 0;
        end
    end
end
figure()
imshow(p)
title('mask in frequency domain for truncatin kernel')
truncated_filter_r =  p+fftkernel_r;
figure()
imshow(truncated_filter_r)
title('truncated filter in frequency domain')
spac_trunc_filter_r = ifft2(ifftshift(truncated_filter_r));
figure()
imshow(mat2gray(abs(spac_trunc_filter_r)))
title('truncated filter in space domain')

truncated_filter_g =  p+fftkernel_g;
truncated_filter_b =  p+fftkernel_b;

trunc_fdeb_r = fft1r./truncated_filter_r ;
figure()
imshow(abs(trunc_fdeb_r))
title('truncated filter image red component frequency domain')
shift_trunc_fdeb_r = ifftshift (trunc_fdeb_r);
deb_r = ifft2(shift_trunc_fdeb_r);
deb_r = deb_r(1:800,1:800,:);
figure()
imshow(mat2gray(real(deb_r)))
title('truncated filter image red component')

trunc_fdeb_g = fft2g./truncated_filter_g ;
figure()
imshow(abs(trunc_fdeb_g))
title('truncated filter image green component frequency domain')
shift_trunc_fdeb_g = ifftshift (trunc_fdeb_g);
deb_g = ifft2(shift_trunc_fdeb_g);
deb_g = deb_g(1:800,1:800,:);
figure()
imshow(mat2gray(real(deb_g)))
title('truncated filter image green component')



trunc_fdeb_b = fft3b./truncated_filter_b ;
figure()
imshow(abs(trunc_fdeb_b))
title('truncated filter image blue component frequency domain')
shift_trunc_fdeb_b = ifftshift (trunc_fdeb_b);
deb_b = ifft2(shift_trunc_fdeb_b);
deb_b =  deb_b(1:800,1:800,:);
figure()
imshow(mat2gray(real(deb_b)))
title('truncated filter image blue component')


final_truncated = cat(3,mat2gray(real(deb_r)),mat2gray(real(deb_g)),mat2gray(real(deb_b)));
figure()
imshow(mat2gray(final_truncated))
title('Deblured Image using Truncated filter')
figure()
imshow(log(1+mat2gray(final_truncated)))
title('Deblured Image using Truncated filter and log transform')
pause
close all

%% Weiner filter
D2 = input('enter a value for ratio of power noise to undegraded image')
D1 = (fftkernel_r.*(conj(fftkernel_r)));
D = fftkernel_r.*(D1+D2);
Weiner_filterimage_fftr = (((fftkernel_r.*(conj(fftkernel_r)))).*fft1r)./D;
Weiner_filterimage_fftr = ifftshift(Weiner_filterimage_fftr);
Weiner_filterimager = ifft2(Weiner_filterimage_fftr);
figure()
imshow(mat2gray(real(Weiner_filterimager(1:800,1:800,:))))
title('Wiener filtered red componenet')

D1 = (fftkernel_g.*(conj(fftkernel_g)));
D = fftkernel_g.*(D1+D2);
Weiner_filterimage_fftg = (((fftkernel_g.*(conj(fftkernel_g)))).*fft2g)./D;
Weiner_filterimage_fftg = ifftshift(Weiner_filterimage_fftg);
Weiner_filterimageg = ifft2(Weiner_filterimage_fftg);
figure()
imshow(mat2gray(abs(Weiner_filterimageg(1:800,1:800,:))))
title('Wiener filtered green componenet')

D1 = (fftkernel_b.*(conj(fftkernel_b)));
D = fftkernel_b.*(D1+D2);
Weiner_filterimage_fftb = (((fftkernel_b.*(conj(fftkernel_b)))).*fft3b)./D;
Weiner_filterimage_fftb = ifftshift(Weiner_filterimage_fftb);
Weiner_filterimageb = ifft2(Weiner_filterimage_fftb);
figure()
imshow(mat2gray(real(Weiner_filterimageb(1:800,1:800,:))))
title('Wiener filtered blue componenet')


Final_wiener_image = cat(3,abs(Weiner_filterimager),abs(Weiner_filterimageg),abs(Weiner_filterimageb));
figure()
imshow(mat2gray(Final_wiener_image(1:800,1:800,:)))
title('Result of Weiner Filter')
pause
close all
%% Constrained least square filter
%laplacian filter
p = [0 1 0;1 -4 1;0 1 0];
fft_p = fft2(p,m,n);
shift_fft_p = fftshift(fft_p);
D1C = k.*(conj(k));
gamma = input('enter a value for gamma')
D2C = shift_fft_p.*conj(shift_fft_p);
DC = (D1C+gamma*D2C);
NU = conj(fftkernel_r).*fft1r;
Ftcs_r = NU./DC;
Ftc_r = ifftshift(Ftcs_r);
tc_r = ifft2(Ftc_r);
figure()
imshow(mat2gray(abs(tc_r(1:800,1:800,:))))
title('Constrained filtered Image red component')

NU = conj(fftkernel_g).*fft2g;
Ftcs_g = NU./DC;
Ftc_g = ifftshift(Ftcs_g);
tc_g = ifft2(Ftc_g);
figure()
imshow(mat2gray(abs(tc_g(1:800,1:800,:))))
title('Constrained filtered Image green component')


NU = conj(fftkernel_b).*fft3b;
Ftcs_b = NU./DC;
Ftc_b = ifftshift(Ftcs_b);
tc_b = ifft2(Ftc_b);
figure()
imshow(mat2gray(abs(tc_b(1:800,1:800,:))))
title('Constrained filtered Image blue component')

Final_cons_image = cat(3,abs(tc_r),abs(tc_g),abs(tc_b));
figure()
imshow(mat2gray(Final_cons_image(1:800,1:800,:)))
title('Constained result ')
figure()
vi = (Final_cons_image(1:800,1:800,:));
vi = log(1+vi);
imshow(mat2gray(vi))
title('Constrained result with log transform appled on filtered image')

pause
close all
%% Just for plotting
figure()
subplot(3,2,1)
imshow(org)
title('Original Image')
subplot(3,2,2)
imshow(blrdimg)
title('Blurred Image')
subplot(3,2,3)
imshow(mat2gray(inv_result(1:800,1:800,:)))
title('Inverse')
[PSNR, SSIM] = my_psnr(double(org),inv_result(1:800,1:800,:));
fprintf(' PSNR for inverse filter is %d',PSNR);
fprintf(' SSIM for inverse filter is %d',SSIM);

subplot(3,2,4)
imshow(mat2gray(final_truncated(1:800,1:800,:)))
title('Truncated')
[PSNR, SSIM] = my_psnr(double(org),final_truncated(1:800,1:800,:));
fprintf(' PSNR for truncated inverse filter is %d \n',PSNR);
fprintf(' SSIM for truncated inverse filter is %d \n',SSIM);
subplot(3,2,5)
imshow(mat2gray(Final_wiener_image(1:800,1:800,:)))
title('Wiener')
[PSNR, SSIM] = my_psnr(double(org),Final_wiener_image(1:800,1:800,:));
fprintf(' PSNR for Wiener filter is %d \n',PSNR);
fprintf(' SSIM for Wiener filter is %d \n',SSIM);
subplot(3,2,6)
imshow(mat2gray(vi))
title('Constrained least square')
[PSNR, SSIM] = my_psnr(double(org),vi(1:800,1:800,:));
fprintf(' PSNR for constrained filter is %d \n',PSNR);
fprintf(' SSIM for constrained filter is %d \n',SSIM);
pause
clear all
close all
end
%% Blind Decon
I = imread('my_image.jpg');
figure()
imshow(I)
title('Blurred Image for applying Blind Decon')
PSF = fspecial('gaussian',20,1);
%PSF = fspecial('motion',5,0);
J = deconvlucy(I,PSF);
figure()
imshow((J))
title('Deblurred Image')

