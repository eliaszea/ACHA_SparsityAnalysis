% test_wa2_reconstruction.m -
% image reconstruction using wave atoms and wavelets

dwtmode('per');

for runid=1:2
  if(runid==1)
    x = double(imread('data_seismic.jpg'));
    ci = 1536;
  else
    x = double(imread('data_fingerprint.jpg'));
    ci = 512;
  end
  
  if(length(size(x))==3)
    x = x/max(x(:));      x = rgb2gray(x);
  end
  x = (x-min(min(x)))/max(max(x-min(min(x))));
  
  N = size(x,1);
  
  figure; clf; imagesc(x); axis equal; axis tight; colormap gray; colorbar;
  title('input image');
  
  %wave atom frame reconstruct
  pat = 'p';
  tp = 'directional';
  c = real(fwa2sym(x,pat,tp));
  cfs = c(:); cfs = sort(abs(cfs)); cfs = cfs(end:-1:1);
  val = cfs(ci);
  d = c .* (abs(c)>=val);
  y = real(iwa2sym(d,pat,tp));
  figure; clf; imagesc(y); axis equal; axis tight; colormap gray; colorbar;
  title(sprintf('wave atom frame reconstruction with %d coefficients',ci));
  
  %wavelet reconstruct
  [c,s] = wavedec2(x, log2(N), 'db5');
  cfs = c(:); cfs = sort(abs(cfs)); cfs = cfs(end:-1:1);
  val = cfs(ci);
  d = c .* (abs(c)>=val);
  y = waverec2(d, s, 'db5');
  figure; clf; imagesc(y); axis equal; axis tight; colormap gray; colorbar;
  title(sprintf('wavelet reconstruction with %d coefficients',ci));
end

