function x = iwa2(c,pat,tp)
% iwa2 - 2D inverse wave atom transform
% -----------------
% INPUT
% --
% c is a cell array which contains the wave atom coefficients. If
% tp=='ortho', then c{j}{m1,m2}(n1,n2) is the coefficient at scale j,
% frequency index (m1,m2) and spatial index (n1,n2). If
% tp=='directional', then c{j,d}{m1,m2}(n1,n2) with d=1,2 are the
% coefficients at scale j, frequency index (m1,m2) and spatial index
% (n1,n2). If tp=='complex', then c{j,d}{m1,m2)(n1,n2) with d=1,2,3,4
% are the coefficients at scale j, frequency index (m1,m2) and spatial
% index (n1,n2).
% --
% pat specifies the type of frequency partition which satsifies
% parabolic scaling relationship. pat can either be 'p' or 'q'.
% --
% tp is the type of tranform.
% 	'ortho': orthobasis
% 	'directional': real-valued frame with single oscillation direction
% 	'complex': complex-valued frame
% -----------------
% OUTPUT
% --
% x is a real N-by-N matrix. N is a power of 2.
% -----------------
% Written by Lexing Ying and Laurent Demanet, 2007
  
  if( ismember(tp, {'ortho','directional','complex'})==0 | ismember(pat, {'p','q','u'})==0 )    error('wrong');  end
  
  if(strcmp(tp, 'ortho')==1)
    %---------------------------------------------------------
    T = 0;
    for s=1:length(c)
      nw = length(c{s});
      for I=1:nw
        for J=1:nw
          T = T + prod(size(c{s}{I,J}));
        end
      end
    end
    N = sqrt(T);
    H = N/2;
    lst = freq_pat(H,pat);
    A = N;
    f = zeros(A,A);
    %------------------
    for s=1:length(lst)
      nw = length(lst{s});
      for I=0:nw-1
        for J=0:nw-1
          if(~isempty(c{s}{I+1,J+1}))
            B = 2^(s-1);
            D = 2*B;
            Ict = I*B;      Jct = J*B; %starting position in freq
            if(mod(I,2)==0)
              Ifm = Ict-2/3*B;        Ito = Ict+4/3*B;
            else
              Ifm = Ict-1/3*B;        Ito = Ict+5/3*B;
            end
            if(mod(J,2)==0)
              Jfm = Jct-2/3*B;        Jto = Jct+4/3*B;
            else
              Jfm = Jct-1/3*B;        Jto = Jct+5/3*B;
            end
            res = fft2(c{s}{I+1,J+1}) / sqrt(prod(size(c{s}{I+1,J+1}))); %res = zeros(D,D);
            for id=0:1
              if(id==0)
                Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
              else
                Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
              end
              for jd=0:1
                if(jd==0)
                  Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
                else
                  Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
                end
                f(mod(Idx,A)+1,mod(Jdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1) + ( Icf.'*Jcf ) .* res(mod(Idx,D)+1,mod(Jdx,D)+1);
              end          
            end
          end
        end
      end
    end
    %------------------
    x = ifft2(f) * sqrt(prod(size(f)));
    
  elseif(strcmp(tp, 'directional')==1)
    %---------------------------------------------------------
    c1 = c(:,1);
    c2 = c(:,2);
    
    T = 0;
    for s=1:length(c1)
      nw = length(c1{s});
      for I=1:nw
        for J=1:nw
          T = T + prod(size(c1{s}{I,J}));
        end
      end
    end
    N = sqrt(T);
    H = N/2;
    lst = freq_pat(H,pat);
    A = N;
    f = zeros(A,A);
    %------------------
    for s=1:length(lst)
      nw = length(lst{s});
      for I=0:nw-1
        for J=0:nw-1
          if(~isempty(c1{s}{I+1,J+1}))
            B = 2^(s-1);
            D = 2*B;
            Ict = I*B;      Jct = J*B; %starting position in freq
            if(mod(I,2)==0)
              Ifm = Ict-2/3*B;        Ito = Ict+4/3*B;
            else
              Ifm = Ict-1/3*B;        Ito = Ict+5/3*B;
            end
            if(mod(J,2)==0)
              Jfm = Jct-2/3*B;        Jto = Jct+4/3*B;
            else
              Jfm = Jct-1/3*B;        Jto = Jct+5/3*B;
            end
            
            res = fft2(c1{s}{I+1,J+1}) / sqrt(prod(size(c1{s}{I+1,J+1}))); %res = zeros(D,D);
            Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
            Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
            f(mod(Idx,A)+1,mod(Jdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1) + ( Icf.'*Jcf ) .* res(mod(Idx,D)+1,mod(Jdx,D)+1);
            Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
            Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
            f(mod(Idx,A)+1,mod(Jdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1) + ( Icf.'*Jcf ) .* res(mod(Idx,D)+1,mod(Jdx,D)+1);
            
            res = fft2(c2{s}{I+1,J+1}) / sqrt(prod(size(c2{s}{I+1,J+1}))); %res = zeros(D,D);
            Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
            Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
            f(mod(Idx,A)+1,mod(Jdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1) + ( Icf.'*Jcf ) .* res(mod(Idx,D)+1,mod(Jdx,D)+1);
            Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
            Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
            f(mod(Idx,A)+1,mod(Jdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1) + ( Icf.'*Jcf ) .* res(mod(Idx,D)+1,mod(Jdx,D)+1);
          end
        end
      end
    end
    x = ifft2(f) * sqrt(prod(size(f)));
    
  elseif(strcmp(tp, 'complex')==1)
    %---------------------------------------------------------
    c1 = c(:,1);
    c2 = c(:,2);
    c3 = c(:,3);
    c4 = c(:,4);
    
    T = 0;
    for s=1:length(c1)
      nw = length(c1{s});
      for I=1:nw
        for J=1:nw
          T = T + prod(size(c1{s}{I,J}));
        end
      end
    end
    N = sqrt(T);
    H = N/2;
    lst = freq_pat(H,pat);
    A = N;
    f = zeros(A,A);
    %------------------
    for s=1:length(lst)
      nw = length(lst{s});
      for I=0:nw-1
        for J=0:nw-1
          if(~isempty(c1{s}{I+1,J+1}))
            B = 2^(s-1);
            D = 2*B;
            Ict = I*B;      Jct = J*B; %starting position in freq
            if(mod(I,2)==0)
              Ifm = Ict-2/3*B;        Ito = Ict+4/3*B;
            else
              Ifm = Ict-1/3*B;        Ito = Ict+5/3*B;
            end
            if(mod(J,2)==0)
              Jfm = Jct-2/3*B;        Jto = Jct+4/3*B;
            else
              Jfm = Jct-1/3*B;        Jto = Jct+5/3*B;
            end
            
            res = fft2(c1{s}{I+1,J+1}) / sqrt(prod(size(c1{s}{I+1,J+1}))); %res = zeros(D,D);
            Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
            Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
            f(mod(Idx,A)+1,mod(Jdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1) + ( Icf.'*Jcf ) .* res(mod(Idx,D)+1,mod(Jdx,D)+1);
            
            res = fft2(c2{s}{I+1,J+1}) / sqrt(prod(size(c2{s}{I+1,J+1}))); %res = zeros(D,D);
            Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
            Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
            f(mod(Idx,A)+1,mod(Jdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1) + ( Icf.'*Jcf ) .* res(mod(Idx,D)+1,mod(Jdx,D)+1);
            
            res = fft2(c3{s}{I+1,J+1}) / sqrt(prod(size(c3{s}{I+1,J+1}))); %res = zeros(D,D);
            Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
            Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
            f(mod(Idx,A)+1,mod(Jdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1) + ( Icf.'*Jcf ) .* res(mod(Idx,D)+1,mod(Jdx,D)+1);
            
            res = fft2(c4{s}{I+1,J+1}) / sqrt(prod(size(c4{s}{I+1,J+1}))); %res = zeros(D,D);
            Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
            Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
            f(mod(Idx,A)+1,mod(Jdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1) + ( Icf.'*Jcf ) .* res(mod(Idx,D)+1,mod(Jdx,D)+1);
          end
        end
      end
    end
    x = ifft2(f) * sqrt(prod(size(f)));
  end
  