function [PSNR, SSIM] = my_psnr(Original_image,Estimated_image)
c1 = 10;
c2 = 5;
[r1,co1,dim1] = size(Original_image);
MSE_image = sqrt(sum(sum((Original_image-Estimated_image).^2))/r1*c1);
muxr = sum(sum(Original_image(:,:,1)))/(r1*co1);
muxg = sum(sum(Original_image(:,:,2)))/(r1*co1);
muxb = sum(sum(Original_image(:,:,3)))/(r1*co1);



[r2,co2,dim2] = size(Estimated_image);
muyr = sum(sum(Estimated_image(:,:,1)))/(r2*co2);
NUMr = max(max(Original_image(:,:,1)));
muyg = sum(sum(Estimated_image(:,:,2)))/(r2*co2);
NUMg = max(max(Original_image(:,:,2)));
muyb = sum(sum(Estimated_image(:,:,3)))/(r2*co2);
NUMb = max(max(Original_image(:,:,3)));

PSNRr = 20*log(NUMr/MSE_image(1,1,1));
PSNRg = 20*log(NUMg/MSE_image(1,1,2));
PSNRb = 20*log(NUMb/MSE_image(1,1,3));
PSNR = (PSNRr+PSNRg+PSNRb)/3;

Ar = cov(Original_image(:,:,1),Estimated_image(:,:,1));
sigmaxyr = Ar(1,2);
sigmaxr = Ar(1,1);
sigmayr = Ar(2,2);
SSIMr = ((2*muxr*muyr+c1)*(2*sigmaxyr+c2))/((muxr.^2+muyr.^2+c1)*(sigmaxr^2 + sigmayr^2 +c2));

Ag = cov(Original_image(:,:,2),Estimated_image(:,:,2));
sigmaxyg = Ag(1,2);
sigmaxg = Ag(1,1);
sigmayg = Ag(2,2);
SSIMg = ((2*muxg*muyg+c1)*(2*sigmaxyg+c2))/((muxg.^2+muyg.^2+c1)*(sigmaxg^2 + sigmayg^2 +c2));

Ab = cov(Original_image(:,:,3),Estimated_image(:,:,3));
sigmaxyb = Ab(1,2);
sigmaxb = Ab(1,1);
sigmayb = Ab(2,2);
SSIMb = ((2*muxb*muyb+c1)*(2*sigmaxyb+c2))/((muxb.^2+muyb.^2+c1)*(sigmaxb^2 + sigmayb^2 +c2));

SSIM = (SSIMr+SSIMg+SSIMb)/3;

end


