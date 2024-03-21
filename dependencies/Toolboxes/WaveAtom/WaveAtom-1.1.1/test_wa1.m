% test_wa1.m - this is a demo file which plots several 1d wave atoms at different scales

N = 512;
pat = 'p';
tp = 'ortho';

z = zeros(N,1);
c = fwa1(z,pat,tp);
ns = size(c,1);

for is=ns-2:ns
  d = c;
  nw = size(d{is},1);
  nb = size(d{is}{nw},1);
  ia = ceil(nw/2);
  d{is,1}{ia}(nb/2+1) = 1;
  y = iwa1(d,pat,tp);
  figure; plot(real(y));  title(sprintf('j=%d, m=(%d), spatial domain',is,ia));
  figure; plot(abs(fftshift(fft2(y))));  title(sprintf('j=%d, m=(%d), frequency domain',is,ia));
end
