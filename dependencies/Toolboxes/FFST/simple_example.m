% simple example for FFST
% computes the shearlet transform of some geometric image
%
%--------------------------------------------------------------------------
% Sören Häuser ~ FFST ~ 2014-07-22 ~ last edited: 2014-07-22 (Sören Häuser)

% create image
A = myPicture2();

% shearlet transform
[ST, Psi] = shearletTransformSpect(A);

% inverse shearlet transform
C = inverseShearletTransformSpect(ST,Psi);

% plot results
subplot(2,2,1)
imagesc(A)
axis image off
colormap(gray)
title('original image')

subplot(2,2,2)
imagesc(abs(ST(:,:,18)))
axis image off
colormap(gray)
title('shearlet coefficients')

subplot(2,2,3)
imagesc(Psi(:,:,18))
axis image off
colormap(gray)
title('shearlet')

subplot(2,2,4)
imagesc(C)
axis image off
colormap(gray)
title('reconstructed image')