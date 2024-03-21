if(1)
  N = 512;
  pat = 'p';
  tp = 'ortho';
  
  z = zeros(N,1);
  c = fwa1sym(z,pat,tp);
  [x1,ww,k1,bb] = pwa1sym(N,pat,tp);
  
  for idx=1:N
    d = zeros(size(c));
    d(idx) = 1;
    x = iwa1sym(d,pat,tp);
    plot(x);
    fprintf(1, '%d %d %d %d\n',x1(idx),ww(idx),k1(idx),bb(idx));
    pause;
  end
end

if(0)
  N = 512;
  pat = 'p';
  tp = 'ortho';
  
  z = zeros(N,1);
  c = fwa1(z,pat,tp);
  ns = size(c,1);
  
  [x1,k1] = pwa1(N,pat,tp);
  
  for is=ns-2:ns
    d = c;
    nw = size(d{is},1);
    nb = size(d{is}{nw},1);
    ia = ceil(nw/2);
    d{is,1}{ia}(nb/2+1) = 1;
    y = iwa1(d,pat,tp);
    figure; plot([0:N-1]/N, real(y));  title(sprintf('j=%d, m=(%d), spatial domain',is,ia));
    figure; plot([-N/2:N/2-1], abs(fftshift(fft2(y))));  title(sprintf('j=%d, m=(%d), frequency domain',is,ia));
    
    x1{is,1}{ia}(nb/2+1)
    k1{is,1}{ia}(nb/2+1)
    pause;
  end
end

if(0)
  N = 512;
  pat = 'p';
  tp = 'ortho';
  
  z = zeros(N,N);
  c = fwa2(z,pat,tp);
  ns = size(c,1);
  
  [x1,x2,k1,k2] = pwa2(N,pat,tp);
  
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
    
    [ x1{is,1}{ia,ib}(nb/2+1,nb/2+1) x2{is,1}{ia,ib}(nb/2+1,nb/2+1) ]
    [ k1{is,1}{ia,ib}(nb/2+1,nb/2+1) k2{is,1}{ia,ib}(nb/2+1,nb/2+1) ]
    pause;
  end

  [a1,a2,b1,b2] = pwa2sym(N,pat,tp);
end

if(0)
  N = 128;
  pat = 'p';
  tp = 'ortho';
  
  z = zeros([N,N,N]);
  c = fwa3(z,pat,tp);
  
  ns = size(c,1);
  [x1,x2,x3,k1,k2,k3] = pwa3(N,pat,tp);
  
  is = ns;
  d = c;
  nw = size(d{is,1},1);
  nb = size(d{is,1}{nw/2,nw/2,nw/2},1);
  d{is,1}{nw/2,nw/2,nw/2}(nb/2+1,nb/2+1,nb/2+1) = 1;
  y = iwa3(d,pat,tp);
  
  [ x1{is,1}{nw/2,nw/2,nw/2}(nb/2+1,nb/2+1,nb/2+1)  x2{is,1}{nw/2,nw/2,nw/2}(nb/2+1,nb/2+1,nb/2+1)  x3{is,1}{nw/2,nw/2,nw/2}(nb/2+1,nb/2+1,nb/2+1) ]
  [ k1{is,1}{nw/2,nw/2,nw/2}(nb/2+1,nb/2+1,nb/2+1)  k2{is,1}{nw/2,nw/2,nw/2}(nb/2+1,nb/2+1,nb/2+1)  k3{is,1}{nw/2,nw/2,nw/2}(nb/2+1,nb/2+1,nb/2+1) ]
  
  
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
  
end