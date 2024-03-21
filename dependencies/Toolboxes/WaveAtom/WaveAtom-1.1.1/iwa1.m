function x = iwa1(c,pat,tp)
% iwa1 - inverse wave atom transform
% -----------------
% INPUT
% --
% c is a cell array which contains the wave atom coefficients. If
% tp=='ortho', then c{j}{m}(n) is the coefficient at scale j, frequency
% index m and spatial index n If tp=='complex', then c{j,1}c{j}{m}(n)
% and c{j,2}c{j}{m}(n) are the coefficients at scale j, frequency index
% m and spatial index n.
% --
% pat specifies the type of frequency partition which satsifies
% parabolic scaling relationship. pat can either be 'p' or 'q'.
% --
% tp is the type of tranform.
% 	'ortho': orthobasis
% 	'complex': complex-valued frame with redunancy 2.
% -----------------
% OUTPUT
% --
% x is a real N-by-1 vector. N is a power of 2.
% -----------------
% Written by Lexing Ying and Laurent Demanet, 2007

  if( ismember(tp, {'ortho','directional','complex'})==0 | ismember(pat, {'p','q','u'})==0 )    error('wrong');  end
  
  if(strcmp(tp, 'ortho')==1)
    %---------------------------------------------------------
    T = 0;
    for s=1:length(c)
      nw = length(c{s});
      for I=1:nw
        T = T + prod(size(c{s}{I}));
      end
    end
    N = T;
    H = N/2;
    lst = freq_pat(H,pat);
    A = N;
    f = zeros(A,1);
    %------------------
    for s=1:length(lst)
      nw = length(lst{s});
      for I=0:nw-1
        if(~isempty(c{s}{I+1}))
          B = 2^(s-1);
          D = 2*B;
          Ict = I*B;
          if(mod(I,2)==0)
            Ifm = Ict-2/3*B;        Ito = Ict+4/3*B;
          else
            Ifm = Ict-1/3*B;        Ito = Ict+5/3*B;
          end
          res = fft(c{s}{I+1}) / sqrt(prod(size(c{s}{I+1})));
          for id=0:1
            if(id==0)
              Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
            else
              Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
            end
            f(mod(Idx,A)+1) = f(mod(Idx,A)+1) + ( Icf.' ) .* res(mod(Idx,D)+1);
          end
        end
      end
    end
    %------------------
    x = ifft(f) * sqrt(prod(size(f)));
  elseif(strcmp(tp, 'directional')==1)
    %---------------------------------------------------------
    error('wrong argument for tp');
  elseif(strcmp(tp, 'complex')==1)
    %---------------------------------------------------------
    c1 = c(:,1);
    c2 = c(:,2);
    
    T = 0;
    for s=1:length(c1)
      nw = length(c1{s});
      for I=1:nw
        T = T + prod(size(c1{s}{I}));
      end
    end
    N = T;
    H = N/2;
    lst = freq_pat(H,pat);
    A = N;
    f = zeros(A,1);
    %------------------
    for s=1:length(lst)
      nw = length(lst{s});
      for I=0:nw-1
        if(~isempty(c1{s}{I+1}))
          B = 2^(s-1);
          D = 2*B;
          
          Ict = I*B;
          if(mod(I,2)==0)
            Ifm = Ict-2/3*B;        Ito = Ict+4/3*B;
          else
            Ifm = Ict-1/3*B;        Ito = Ict+5/3*B;
          end
          
          Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
          res = fft(c1{s}{I+1}) / sqrt(prod(size(c1{s}{I+1})));
          f(mod(Idx,A)+1) = f(mod(Idx,A)+1) + ( Icf.' ) .* res(mod(Idx,D)+1);
          
          Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
          res = fft(c2{s}{I+1}) / sqrt(prod(size(c2{s}{I+1})));
          f(mod(Idx,A)+1) = f(mod(Idx,A)+1) + ( Icf.' ) .* res(mod(Idx,D)+1);
        end
      end
    end
    %------------------
    x = ifft(f) * sqrt(prod(size(f)));
  end
  