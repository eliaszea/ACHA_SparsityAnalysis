function x = iwa3(c,pat,tp)
% iwa3 - 3D inverse wave atom transform
% -----------------
% INPUT
% --
% c is a cell array which contains the wave atom coefficients. If
% tp=='ortho', then c{j}{m1,m2,m3}(n1,n2,n3) is the coefficient at scale
% j, frequency index (m1,m2,m3) and spatial index (n1,n2,n3). If
% tp=='directional', then c{j,d}{m1,m2,m3}(n1,n2,n3) with d=1,2,3,4 are
% the coefficients at scale j, frequency index (m1,m2,m3) and spatial
% index (n1,n2,n3). If tp=='complex', then c{j,d}{m1,m2,m3)(n1,n2,n3)
% with d=1,2,3,4,5,6,7,8 are the coefficients at scale j, frequency
% index (m1,m2,m3) and spatial index (n1,n2,n3).
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
% x is a real N-by-N-by-N array. N is a power of 2.
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
          for K=1:nw
            T = T + prod(size(c{s}{I,J,K}));
          end
        end
      end
    end
    N = round(T^(1/3));
    H = N/2;
    lst = freq_pat(H,pat);
    A = N;
    f = zeros(A,A,A);
    %------------------
    for s=1:length(lst)
      nw = length(lst{s});
      for I=0:nw-1
        for J=0:nw-1
          for K=0:nw-1
            if(~isempty(c{s}{I+1,J+1,K+1}))
              B = 2^(s-1);
              D = 2*B;
              Ict = I*B;              Jct = J*B;              Kct = K*B;
              %starting position in freq
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
              if(mod(K,2)==0)
                Kfm = Kct-2/3*B;        Kto = Kct+4/3*B;
              else
                Kfm = Kct-1/3*B;        Kto = Kct+5/3*B;
              end
              res = fftn(c{s}{I+1,J+1,K+1}) / sqrt(prod(size(c{s}{I+1,J+1,K+1})));
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
                  for kd=0:1
                    if(kd==0)
                      Kdx = [ceil(Kfm):floor(Kto)];      Kcf = kf_rt(Kdx/B*pi, K);
                    else
                      Kdx = [ceil(-Kto):floor(-Kfm)];      Kcf = kf_lf(Kdx/B*pi, K);
                    end
                    f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) + kron3(Icf,Jcf,Kcf) .* res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1);
                  end          
                end
              end
            end
          end
        end
      end
    end
    %------------------
    x = ifftn(f) * sqrt(prod(size(f)));
    
  elseif(strcmp(tp, 'directional')==1)
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
          for K=1:nw
            T = T + prod(size(c1{s}{I,J,K}));
          end
      end
      end
    end
    N = round(T^(1/3));
    H = N/2;
    lst = freq_pat(H,pat);
    A = N;
    f = zeros(A,A,A);
    %------------------
    for s=1:length(lst)
      nw = length(lst{s});
      for I=0:nw-1
        for J=0:nw-1
          for K=0:nw-1
            if(~isempty(c1{s}{I+1,J+1,K+1}))
              B = 2^(s-1);
              D = 2*B;
              Ict = I*B;              Jct = J*B;              Kct = K*B;
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
              if(mod(K,2)==0)
                Kfm = Kct-2/3*B;        Kto = Kct+4/3*B;
              else
                Kfm = Kct-1/3*B;        Kto = Kct+5/3*B;
              end
              
              res = fftn(c1{s}{I+1,J+1,K+1}) / sqrt(prod(size(c1{s}{I+1,J+1,K+1})));
              Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
              Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
              Kdx = [ceil(Kfm):floor(Kto)];      Kcf = kf_rt(Kdx/B*pi, K);
              f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) + kron3(Icf,Jcf,Kcf) .* res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1);
              Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
              Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
              Kdx = [ceil(-Kto):floor(-Kfm)];      Kcf = kf_lf(Kdx/B*pi, K);
              f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) + kron3(Icf,Jcf,Kcf) .* res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1);
              
              res = fftn(c2{s}{I+1,J+1,K+1}) / sqrt(prod(size(c2{s}{I+1,J+1,K+1})));
              Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
              Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
              Kdx = [ceil(-Kto):floor(-Kfm)];      Kcf = kf_lf(Kdx/B*pi, K);
              f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) + kron3(Icf,Jcf,Kcf) .* res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1);
              Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
              Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
              Kdx = [ceil(Kfm):floor(Kto)];      Kcf = kf_rt(Kdx/B*pi, K);
              f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) + kron3(Icf,Jcf,Kcf) .* res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1);
              
              res = fftn(c3{s}{I+1,J+1,K+1}) / sqrt(prod(size(c3{s}{I+1,J+1,K+1})));
              Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
              Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
              Kdx = [ceil(Kfm):floor(Kto)];      Kcf = kf_rt(Kdx/B*pi, K);
              f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) + kron3(Icf,Jcf,Kcf) .* res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1);
              Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
              Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
              Kdx = [ceil(-Kto):floor(-Kfm)];      Kcf = kf_lf(Kdx/B*pi, K);
              f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) + kron3(Icf,Jcf,Kcf) .* res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1);
              
              res = fftn(c4{s}{I+1,J+1,K+1}) / sqrt(prod(size(c4{s}{I+1,J+1,K+1})));
              Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
              Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
              Kdx = [ceil(-Kto):floor(-Kfm)];      Kcf = kf_lf(Kdx/B*pi, K);
              f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) + kron3(Icf,Jcf,Kcf) .* res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1);
              Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
              Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
              Kdx = [ceil(Kfm):floor(Kto)];      Kcf = kf_rt(Kdx/B*pi, K);
              f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) + kron3(Icf,Jcf,Kcf) .* res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1);
            end
          end
        end
      end
    end
    %------------------
    x = ifftn(f) * sqrt(prod(size(f)));
    
  elseif(strcmp(tp, 'complex')==1)
    %---------------------------------------------------------
    c1 = c(:,1);
    c2 = c(:,2);
    c3 = c(:,3);
    c4 = c(:,4);
    c5 = c(:,5);
    c6 = c(:,6);
    c7 = c(:,7);
    c8 = c(:,8);
    
    T = 0;
    for s=1:length(c1)
      nw = length(c1{s});
      for I=1:nw
        for J=1:nw
          for K=1:nw
            T = T + prod(size(c1{s}{I,J,K}));
          end
        end
      end
    end
    N = round(T^(1/3));
    H = N/2;
    lst = freq_pat(H,pat);
    A = N;
    f = zeros(A,A,A);
    %------------------
    for s=1:length(lst)
      nw = length(lst{s});
      for I=0:nw-1
        for J=0:nw-1
          for K=0:nw-1
            if(~isempty(c1{s}{I+1,J+1,K+1}))
              B = 2^(s-1);
              D = 2*B;
              Ict = I*B;              Jct = J*B;              Kct = K*B;
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
              if(mod(K,2)==0)
                Kfm = Kct-2/3*B;        Kto = Kct+4/3*B;
              else
                Kfm = Kct-1/3*B;        Kto = Kct+5/3*B;
              end
              
              res = fftn(c1{s}{I+1,J+1,K+1}) / sqrt(prod(size(c1{s}{I+1,J+1,K+1})));
              Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
              Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
              Kdx = [ceil(Kfm):floor(Kto)];      Kcf = kf_rt(Kdx/B*pi, K);
              f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) + kron3(Icf,Jcf,Kcf) .* res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1);
              
              res = fftn(c2{s}{I+1,J+1,K+1}) / sqrt(prod(size(c2{s}{I+1,J+1,K+1})));
              Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
              Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
              Kdx = [ceil(-Kto):floor(-Kfm)];      Kcf = kf_lf(Kdx/B*pi, K);
              f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) + kron3(Icf,Jcf,Kcf) .* res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1);
              
              res = fftn(c3{s}{I+1,J+1,K+1}) / sqrt(prod(size(c3{s}{I+1,J+1,K+1})));
              Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
              Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
              Kdx = [ceil(-Kto):floor(-Kfm)];      Kcf = kf_lf(Kdx/B*pi, K);
              f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) + kron3(Icf,Jcf,Kcf) .* res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1);
              
              res = fftn(c4{s}{I+1,J+1,K+1}) / sqrt(prod(size(c4{s}{I+1,J+1,K+1})));
              Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
              Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
              Kdx = [ceil(Kfm):floor(Kto)];      Kcf = kf_rt(Kdx/B*pi, K);
              f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) + kron3(Icf,Jcf,Kcf) .* res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1);

              res = fftn(c5{s}{I+1,J+1,K+1}) / sqrt(prod(size(c5{s}{I+1,J+1,K+1})));
              Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
              Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
              Kdx = [ceil(Kfm):floor(Kto)];      Kcf = kf_rt(Kdx/B*pi, K);
              f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) + kron3(Icf,Jcf,Kcf) .* res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1);

              res = fftn(c6{s}{I+1,J+1,K+1}) / sqrt(prod(size(c6{s}{I+1,J+1,K+1})));
              Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
              Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
              Kdx = [ceil(-Kto):floor(-Kfm)];      Kcf = kf_lf(Kdx/B*pi, K);
              f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) + kron3(Icf,Jcf,Kcf) .* res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1);
              
              res = fftn(c7{s}{I+1,J+1,K+1}) / sqrt(prod(size(c7{s}{I+1,J+1,K+1})));
              Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
              Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
              Kdx = [ceil(-Kto):floor(-Kfm)];      Kcf = kf_lf(Kdx/B*pi, K);
              f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) + kron3(Icf,Jcf,Kcf) .* res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1);
              
              res = fftn(c8{s}{I+1,J+1,K+1}) / sqrt(prod(size(c8{s}{I+1,J+1,K+1})));
              Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
              Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
              Kdx = [ceil(Kfm):floor(Kto)];      Kcf = kf_rt(Kdx/B*pi, K);
              f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1) + kron3(Icf,Jcf,Kcf) .* res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1);
            end
          end
        end
      end
    end
    %------------------
    x = ifftn(f) * sqrt(prod(size(f)));
  end
  
  
function M = kron3(I,J,K)
  tmp = I.'*J;
  tmp = tmp(:)*K;
  M = reshape(tmp,[length(I),length(J),length(K)]);
  
  
  
  
  
  