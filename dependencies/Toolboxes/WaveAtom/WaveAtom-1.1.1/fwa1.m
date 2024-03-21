function c = fwa1(x,pat,tp)
% fwa1 - forward wave atom transform
% -----------------
% INPUT
% --
% x is a real N-by-1 vector. N is a power of 2.
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
% c is a cell array which contains the wave atom coefficients. If
% tp=='ortho', then c{j}{m}(n) is the coefficient at scale j, frequency
% index m and spatial index n If tp=='complex', then c{j,1}c{j}{m}(n)
% and c{j,2}c{j}{m}(n) are the coefficients at scale j, frequency index
% m and spatial index n.
% -----------------
% Written by Lexing Ying and Laurent Demanet, 2007
  
  if( ismember(tp, {'ortho','directional','complex'})==0 | ismember(pat, {'p','q','u'})==0 )    error('wrong');  end
  x = x(:);
  
  if(strcmp(tp, 'ortho')==1)
    %---------------------------------------------------------
    N = length(x);
    H = N/2;
    lst = freq_pat(H,pat);
    %------------------
    f = fft(x) / sqrt(length(x));
    A = N;
    c = cell(length(lst),1);
    %------------------
    for s=1:length(lst)
      nw = length(lst{s});
      c{s} = cell(nw,1);
      for I=0:nw-1
        if(lst{s}(I+1)==0)
          c{s}{I+1} = [];
        else
          B = 2^(s-1);
          D = 2*B;
          Ict = I*B;
          if(mod(I,2)==0)
            Ifm = Ict-2/3*B;        Ito = Ict+4/3*B;
          else
            Ifm = Ict-1/3*B;        Ito = Ict+5/3*B;
          end
          res = zeros(D,1);
          for id=0:1
            if(id==0)
              Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
            else
              Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
            end
            res(mod(Idx,D)+1) = res(mod(Idx,D)+1) + conj( Icf.' ) .* f(mod(Idx,A)+1);
          end
          c{s}{I+1} = ifft(res) * sqrt(prod(size(res)));
        end
      end
    end
  elseif(strcmp(tp, 'directional')==1)
    %---------------------------------------------------------
    error('wrong argument for tp');
  elseif(strcmp(tp, 'complex')==1)
    %---------------------------------------------------------
    N = length(x);
    H = N/2;
    lst = freq_pat(H,pat);
    %------------------
    f = fft(x) / sqrt(length(x));
    A = N;
    c1 = cell(length(lst),1);
    c2 = cell(length(lst),1);
    %------------------
    for s=1:length(lst)
      nw = length(lst{s});
      c1{s} = cell(nw,1);
      c2{s} = cell(nw,1);
      for I=0:nw-1
        if(lst{s}(I+1)==0)
          c1{s}{I+1} = [];
          c2{s}{I+1} = [];
        else
          B = 2^(s-1);
          D = 2*B;
          
          Ict = I*B;
          if(mod(I,2)==0)
            Ifm = Ict-2/3*B;        Ito = Ict+4/3*B;
          else
            Ifm = Ict-1/3*B;        Ito = Ict+5/3*B;
          end
          
          res = zeros(D,1);
          Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
          res(mod(Idx,D)+1) = res(mod(Idx,D)+1) + conj( Icf.' ) .* f(mod(Idx,A)+1);
          c1{s}{I+1} = ifft(res) * sqrt(prod(size(res)));
          
          res = zeros(D,1);
          Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
          res(mod(Idx,D)+1) = res(mod(Idx,D)+1) + conj( Icf.' ) .* f(mod(Idx,A)+1);
          c2{s}{I+1} = ifft(res) * sqrt(prod(size(res)));
        end
      end
    end
    c = [c1 c2];
  end
  