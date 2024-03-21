function c = fwa3(x,pat,tp)
% fwa3 - 3D forward wave atom transform
% -----------------
% INPUT
% --
% x is a real N-by-N-by-N array. N is a power of 2.
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
% c is a cell array which contains the wave atom coefficients. If
% tp=='ortho', then c{j}{m1,m2,m3}(n1,n2,n3) is the coefficient at scale
% j, frequency index (m1,m2,m3) and spatial index (n1,n2,n3). If
% tp=='directional', then c{j,d}{m1,m2,m3}(n1,n2,n3) with d=1,2,3,4 are
% the coefficients at scale j, frequency index (m1,m2,m3) and spatial
% index (n1,n2,n3). If tp=='complex', then c{j,d}{m1,m2,m3)(n1,n2,n3)
% with d=1,2,3,4,5,6,7,8 are the coefficients at scale j, frequency
% index (m1,m2,m3) and spatial index (n1,n2,n3).
% -----------------
% Written by Lexing Ying and Laurent Demanet, 2007

  if( ismember(tp, {'ortho','directional','complex'})==0 | ismember(pat, {'p','q','u'})==0 )    error('wrong');  end
  
  if(strcmp(tp, 'ortho')==1)
    %---------------------------------------------------------
    N = size(x,1);
    H = N/2;
    lst = freq_pat(H,pat);
    %------------------
    f = fftn(x) / sqrt(prod(size(x)));
    A = N;
    c = cell(length(lst),1);
    %------------------
    for s=1:length(lst)
      nw = length(lst{s});
      c{s} = cell(nw,nw,nw);
      for I=0:nw-1
        for J=0:nw-1
          for K=0:nw-1
            if(lst{s}(I+1)==0 & lst{s}(J+1)==0 & lst{s}(K+1)==0)
              c{s}{I+1,J+1,K+1} = [];
            else
              B = 2^(s-1);
              D = 2*B;
              Ict = I*B;            Jct = J*B;            Kct = K*B;
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
              res = zeros(D,D,D);
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
                    res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) = res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) + conj( kron3(Icf,Jcf,Kcf) ) .* f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1);
                  end
                end
              end
              c{s}{I+1,J+1,K+1} = ifftn(res) * sqrt(prod(size(res)));
            end
          end
        end
      end
    end
    
  elseif(strcmp(tp, 'directional')==1)
    %---------------------------------------------------------
    N = size(x,1);
    H = N/2;
    lst = freq_pat(H,pat);
    %------------------
    f = fftn(x) / sqrt(prod(size(x)));
    A = N;
    c1 = cell(length(lst),1);
    c2 = cell(length(lst),1);
    c3 = cell(length(lst),1);
    c4 = cell(length(lst),1);
    %------------------
    for s=1:length(lst)
      nw = length(lst{s});
      c1{s} = cell(nw,nw,nw);
      c2{s} = cell(nw,nw,nw);
      c3{s} = cell(nw,nw,nw);
      c4{s} = cell(nw,nw,nw);
      for I=0:nw-1
        for J=0:nw-1
          for K=0:nw-1
            if(lst{s}(I+1)==0 & lst{s}(J+1)==0 & lst{s}(K+1)==0)
              c1{s}{I+1,J+1,K+1} = [];
              c2{s}{I+1,J+1,K+1} = [];
              c3{s}{I+1,J+1,K+1} = [];
              c4{s}{I+1,J+1,K+1} = [];
            else
              B = 2^(s-1);
              D = 2*B;
              Ict = I*B;            Jct = J*B;            Kct = K*B;
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
              
              res = zeros(D,D,D);
              Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
              Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
              Kdx = [ceil(Kfm):floor(Kto)];      Kcf = kf_rt(Kdx/B*pi, K);
              res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) = res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) + conj( kron3(Icf,Jcf,Kcf) ) .* f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1);
              Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
              Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
              Kdx = [ceil(-Kto):floor(-Kfm)];      Kcf = kf_lf(Kdx/B*pi, K);
              res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) = res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) + conj( kron3(Icf,Jcf,Kcf) ) .* f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1);
              c1{s}{I+1,J+1,K+1} = ifftn(res) * sqrt(prod(size(res)));
            
              res = zeros(D,D,D);
              Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
              Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
              Kdx = [ceil(-Kto):floor(-Kfm)];      Kcf = kf_lf(Kdx/B*pi, K);
              res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) = res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) + conj( kron3(Icf,Jcf,Kcf) ) .* f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1);
              Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
              Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
              Kdx = [ceil(Kfm):floor(Kto)];      Kcf = kf_rt(Kdx/B*pi, K);
              res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) = res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) + conj( kron3(Icf,Jcf,Kcf) ) .* f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1);
              c2{s}{I+1,J+1,K+1} = ifftn(res) * sqrt(prod(size(res)));

              res = zeros(D,D,D);
              Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
              Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
              Kdx = [ceil(Kfm):floor(Kto)];      Kcf = kf_rt(Kdx/B*pi, K);
              res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) = res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) + conj( kron3(Icf,Jcf,Kcf) ) .* f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1);
              Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
              Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
              Kdx = [ceil(-Kto):floor(-Kfm)];      Kcf = kf_lf(Kdx/B*pi, K);
              res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) = res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) + conj( kron3(Icf,Jcf,Kcf) ) .* f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1);
              c3{s}{I+1,J+1,K+1} = ifftn(res) * sqrt(prod(size(res)));
              
              res = zeros(D,D,D);
              Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
              Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
              Kdx = [ceil(-Kto):floor(-Kfm)];      Kcf = kf_lf(Kdx/B*pi, K);
              res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) = res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) + conj( kron3(Icf,Jcf,Kcf) ) .* f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1);
              Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
              Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
              Kdx = [ceil(Kfm):floor(Kto)];      Kcf = kf_rt(Kdx/B*pi, K);
              res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) = res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) + conj( kron3(Icf,Jcf,Kcf) ) .* f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1);
              c4{s}{I+1,J+1,K+1} = ifftn(res) * sqrt(prod(size(res)));
            end
          end
        end
      end
    end
    
    c = [c1 c2 c3 c4];
    
  elseif(strcmp(tp, 'complex')==1)
    %---------------------------------------------------------
    N = size(x,1);
    H = N/2;
    lst = freq_pat(H,pat);
    %------------------
    f = fftn(x) / sqrt(prod(size(x)));
    A = N;
    c1 = cell(length(lst),1);
    c2 = cell(length(lst),1);
    c3 = cell(length(lst),1);
    c4 = cell(length(lst),1);
    c5 = cell(length(lst),1);
    c6 = cell(length(lst),1);
    c7 = cell(length(lst),1);
    c8 = cell(length(lst),1);
    %------------------
    for s=1:length(lst)
      nw = length(lst{s});
      c1{s} = cell(nw,nw,nw);
      c2{s} = cell(nw,nw,nw);
      c3{s} = cell(nw,nw,nw);
      c4{s} = cell(nw,nw,nw);
      c5{s} = cell(nw,nw,nw);
      c6{s} = cell(nw,nw,nw);
      c7{s} = cell(nw,nw,nw);
      c8{s} = cell(nw,nw,nw);
      for I=0:nw-1
        for J=0:nw-1
          for K=0:nw-1
            if(lst{s}(I+1)==0 & lst{s}(J+1)==0 & lst{s}(K+1)==0)
              c1{s}{I+1,J+1,K+1} = [];
              c2{s}{I+1,J+1,K+1} = [];
              c3{s}{I+1,J+1,K+1} = [];
              c4{s}{I+1,J+1,K+1} = [];
              c5{s}{I+1,J+1,K+1} = [];
              c6{s}{I+1,J+1,K+1} = [];
              c7{s}{I+1,J+1,K+1} = [];
              c8{s}{I+1,J+1,K+1} = [];
            else
              B = 2^(s-1);
              D = 2*B;
              Ict = I*B;            Jct = J*B;            Kct = K*B;
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
              
              res = zeros(D,D,D);
              Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
              Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
              Kdx = [ceil(Kfm):floor(Kto)];      Kcf = kf_rt(Kdx/B*pi, K);
              res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) = res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) + conj( kron3(Icf,Jcf,Kcf) ) .* f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1);
              c1{s}{I+1,J+1,K+1} = ifftn(res) * sqrt(prod(size(res)));
              
              res = zeros(D,D,D);
              Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
              Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
              Kdx = [ceil(-Kto):floor(-Kfm)];      Kcf = kf_lf(Kdx/B*pi, K);
              res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) = res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) + conj( kron3(Icf,Jcf,Kcf) ) .* f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1);
              c2{s}{I+1,J+1,K+1} = ifftn(res) * sqrt(prod(size(res)));
            
              res = zeros(D,D,D);
              Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
              Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
              Kdx = [ceil(-Kto):floor(-Kfm)];      Kcf = kf_lf(Kdx/B*pi, K);
              res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) = res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) + conj( kron3(Icf,Jcf,Kcf) ) .* f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1);
              c3{s}{I+1,J+1,K+1} = ifftn(res) * sqrt(prod(size(res)));
              
              res = zeros(D,D,D);
              Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
              Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
              Kdx = [ceil(Kfm):floor(Kto)];      Kcf = kf_rt(Kdx/B*pi, K);
              res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) = res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) + conj( kron3(Icf,Jcf,Kcf) ) .* f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1);
              c4{s}{I+1,J+1,K+1} = ifftn(res) * sqrt(prod(size(res)));

              res = zeros(D,D,D);
              Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
              Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
              Kdx = [ceil(Kfm):floor(Kto)];      Kcf = kf_rt(Kdx/B*pi, K);
              res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) = res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) + conj( kron3(Icf,Jcf,Kcf) ) .* f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1);
              c5{s}{I+1,J+1,K+1} = ifftn(res) * sqrt(prod(size(res)));
              
              res = zeros(D,D,D);
              Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
              Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
              Kdx = [ceil(-Kto):floor(-Kfm)];      Kcf = kf_lf(Kdx/B*pi, K);
              res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) = res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) + conj( kron3(Icf,Jcf,Kcf) ) .* f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1);
              c6{s}{I+1,J+1,K+1} = ifftn(res) * sqrt(prod(size(res)));
              
              res = zeros(D,D,D);
              Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
              Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
              Kdx = [ceil(-Kto):floor(-Kfm)];      Kcf = kf_lf(Kdx/B*pi, K);
              res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) = res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) + conj( kron3(Icf,Jcf,Kcf) ) .* f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1);
              c7{s}{I+1,J+1,K+1} = ifftn(res) * sqrt(prod(size(res)));
              
              res = zeros(D,D,D);
              Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
              Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
              Kdx = [ceil(Kfm):floor(Kto)];      Kcf = kf_rt(Kdx/B*pi, K);
              res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) = res(mod(Idx,D)+1,mod(Jdx,D)+1,mod(Kdx,D)+1) + conj( kron3(Icf,Jcf,Kcf) ) .* f(mod(Idx,A)+1,mod(Jdx,A)+1,mod(Kdx,A)+1);
              c8{s}{I+1,J+1,K+1} = ifftn(res) * sqrt(prod(size(res)));
            end
          end
        end
      end
    end
    
    c = [c1 c2 c3 c4 c5 c6 c7 c8];
  end
  
  
function M = kron3(I,J,K)
  tmp = I.'*J;
  tmp = tmp(:)*K;
  M = reshape(tmp,[length(I),length(J),length(K)]);
  
  
