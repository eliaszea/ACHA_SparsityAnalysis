% test_mewa2_reconstruction.m -
% image reconstruction using mirror-extended wave atoms and standard wave atoms

if(1)
  N = 512;
  H = N/2;
  [a,b] = ndgrid(0:N-1);
  x = a+b-N;
  x = (x-min(min(x)))/max(max(x-min(min(x))));
  
  ci = 128;
  %-----------
  pat = 'p';
  tp = 'directional';
  c = mefwa2sym(x,pat,tp);
  cfs = c(:); cfs = sort(abs(cfs)); cfs = cfs(end:-1:1);
  val = cfs(ci);
  gd = abs(c)>=val+eps;  cnt = sum(gd(:));
  d = c .* (abs(c)>=val+eps);
  y = meiwa2sym(d,pat,tp);
  figure;  imagesc(y); axis equal; axis tight; colormap gray; colorbar;
  title('Partial reconstruction with mirror-extended wave atoms');
  fprintf(1,'number of coefficients used for mirror-extended wave atoms = %d\n',cnt);
  %-----------
  pat = 'p';
  tp = 'directional';
  c = fwa2sym(x,pat,tp);
  cfs = c(:); cfs = sort(abs(cfs)); cfs = cfs(end:-1:1);
  val = cfs(ci);
  gd = abs(c)>=val+eps;  cnt = sum(gd(:));
  d = c .* (abs(c)>=val+eps);
  y = iwa2sym(d,pat,tp);
  figure;  imagesc(y); axis equal; axis tight; colormap gray; colorbar;
  title('Partial reconstruction with wave atoms');
  fprintf(1,'number of coefficients used for standard wave atoms = %d\n',cnt);
end

