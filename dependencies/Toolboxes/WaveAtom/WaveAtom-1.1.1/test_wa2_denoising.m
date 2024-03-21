% test_wa2_denoising.m -
% image denoising using wave atoms and wavelets

dwtmode('per');

for runid=1:2
  if(runid==1)
    x = double(imread('data_seismic.jpg'));
    sigma = 0.1;
  else
    x = double(imread('data_fingerprint.jpg'));
    sigma = 0.2;
  end
  
  if(length(size(x))==3)
    x = x/max(x(:));      x = rgb2gray(x);
  end
  x = (x-min(min(x)))/max(max(x-min(min(x))));
  
  N = size(x,1);
  noisy_x = x + sigma * randn(size(x));
  
  figure; clf; imagesc(x); axis equal; axis tight; colormap gray; colorbar;
  title('input image');
  
  figure; clf; imagesc(noisy_x); axis equal; axis tight; colormap gray; colorbar;
  title('noisy image');
  
  threshold = 3 * sigma;
  shifts = [1 2 4];
  
  %wave atom frame denoising
  pat = 'p';
  tp = 'directional';
  val = zeros(size(noisy_x));
  for is=shifts
    for js=shifts
      xx = circshift(noisy_x,[is,js]);
      cfs = real(fwa2sym(xx,pat,tp));
      cfs = cfs .* (abs(cfs)> sqrt(2)/2*threshold);
      tmp = real(iwa2sym(cfs,pat,tp));
      tmp = circshift(tmp,[-is,-js]);
      val = val + tmp;
    end
  end
  val = val / length(shifts)^2;
  figure; clf; imagesc(val); axis equal; axis tight; axis off; colormap gray;
  title(sprintf('wave atom frame, psnr = %.2f', -10*log10(sum(sum((val-x).^2))/(N*N))));
  
  %wavelet denoising
  val = zeros(size(noisy_x));
  for is=shifts
    for js=shifts
      xx = circshift(noisy_x,[is,js]);
      [cfs,s] = wavedec2(xx, log2(N), 'db5');
      cfs = cfs .* (abs(cfs)>threshold);
      tmp = waverec2(cfs,s, 'db5');
      tmp = circshift(tmp,[-is,-js]);
      val = val + tmp;
    end
  end
  val = val / length(shifts)^2;
  figure; clf; imagesc(val); axis equal; axis tight; axis off; colormap gray;
  title(sprintf('wavelet db5, psnr = %.2f', -10*log10(sum(sum((val-x).^2))/(N*N))));
end

