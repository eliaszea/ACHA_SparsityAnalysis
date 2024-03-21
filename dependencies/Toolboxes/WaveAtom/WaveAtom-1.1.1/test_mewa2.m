% test_mewa2.m - this is a demo file which plots several mirror-extended wave atoms at different scales

N = 512;
pat = 'p';
tp = 'directional';

z = zeros(N,N);
c = mefwa2(z,pat,tp);
ns = size(c,1);

for is=ns-2:ns
  nw = size(c{is,1},1);
  m1 = nw-2;  m2 = nw-2;
  [n1,n2] = size(c{is,1}{m1,m2});
  
  d = c;
  d{is,1}{m1,m2}(ceil(n1/4),1) = 1;
  y = meiwa2(d,pat,tp);
  figure; imagesc(real(y)); axis xy; axis equal; axis tight; colormap(1-gray); colorbar;
  title(sprintf('wave atom near the boundary, j=%d, m=(%d,%d)',is,m1,m2));
  
  d = c;
  d{is,1}{m1,m2}(ceil(n1/4),ceil(n2/2)) = 1;
  y = meiwa2(d,pat,tp);
  figure; imagesc(real(y)); axis xy; axis equal; axis tight; colormap(1-gray); colorbar;
  title(sprintf('wave atom at the center, j=%d, m=(%d,%d)',is,m1,m2));
end

