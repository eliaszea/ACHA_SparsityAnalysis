% test_wa2.m - this is a demo file which plots several 2d wave atoms at different scales

N = 512;
pat = 'p';
tp = 'directional';

z = zeros(N,N);
c = fwa2(z,pat,tp);
ns = size(c,1);

for is=ns-2:ns
  d = c;
  nw = size(d{is,1},1);
  nb = size(d{is,1}{nw,nw},1);
  ia = ceil(nw/2);  ib = ceil(nw/4);
  d{is,1}{ia,ib}(nb/2+1,nb/2+1) = 1;
  y = iwa2(d,pat,tp);
  figure; imagesc(real(y)); axis equal; axis tight; colormap(1-gray); colorbar;
  title(sprintf('j=%d, m=(%d,%d), spatial domain',is,ia,ib));
  figure; imagesc(abs(fftshift(fft2(y)))); axis equal; axis tight; colormap(1-gray); colorbar;
  title(sprintf('j=%d, m=(%d,%d), frequency domain',is,ia,ib));
end

