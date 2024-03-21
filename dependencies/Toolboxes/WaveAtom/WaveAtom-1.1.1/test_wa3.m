% test_wa3.m - this is a demo file which plots a 3d wave atom at the finest scale

N = 128;
pat = 'p';
tp = 'directional';

z = zeros([N,N,N]);
c = fwa3(z,pat,tp);

ns = size(c,1);

is = ns;
d = c;
nw = size(d{is,1},1);
nb = size(d{is,1}{nw/2,nw/2,nw/2},1);
d{is,1}{nw/2,nw/2,nw/2}(nb/2+1,nb/2+1,nb/2+1) = 1;
y = iwa3(d,pat,tp);

figure; imagesc(squeeze(real(y(:,:,N/2+1)))); colormap(1-gray); axis equal; axis tight;
title('spatial cross section: 1st and 2nd coordinates');
figure; imagesc(squeeze(real(y(N/2+1,:,:)))); colormap(1-gray); axis equal; axis tight;
title('spatial cross section: 2nd and 3rd coordinates');
figure; imagesc(squeeze(real(y(:,N/2+1,:)))); colormap(1-gray); axis equal; axis tight;
title('spatial cross section: 1st and 3rd coordinates');

f = fftshift(fftn(y));
figure; imagesc(squeeze(sum(abs(f),3))); colormap(1-gray); axis equal; axis tight;
title('frequency domain projection: 1st and 2nd coordinates');
figure; imagesc(squeeze(sum(abs(f),1))); colormap(1-gray); axis equal; axis tight;
title('frequency domain projection: 2nd and 3rd coordinates');
figure; imagesc(squeeze(sum(abs(f),2))); colormap(1-gray); axis equal; axis tight;
title('frequency domain projection: 1st and 3rd coordinates');
