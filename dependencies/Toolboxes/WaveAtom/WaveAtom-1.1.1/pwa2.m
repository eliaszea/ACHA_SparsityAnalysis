function [x1,x2,w,k1,k2,z] = pwa2(N,pat,tp)
% pwa2 - get position information
% -----------------
% INPUT
% --
% N -- size
% --
% pat specifies the type of frequency partition which satsifies
% parabolic scaling relationship. pat can either be 'p' or 'q'.
% --
% tp is the type of tranform.
% 	'ortho': orthobasis
% -----------------
% OUTPUT
% --
% x1,x2 are cell structures that contain the centers of the wave atoms
% in spatial domain.
% --
% w is a cell structure that contains the widths of the wave atoms
% in spatial domain.
% --
% k1,k2 are cell structures that contain the centers of the wave atoms
% in frequency domain.
% --
% z is a cell structure that contains the widths of the wave atoms
% in frequency domain.
% --
% -----------------
% Written by Lexing Ying and Laurent Demanet, 2007

  
  if( ismember(tp, {'ortho','directional','complex'})==0 | ismember(pat, {'p','q','u'})==0 )    error('wrong');  end
  
  if(strcmp(tp,'ortho')==1)
    %---------------------------------------------------------
    H = N/2;
    lst = freq_pat(H,pat);
    x1 = cell(length(lst),1);    x2 = cell(length(lst),1);
    w  = cell(length(lst),1);
    k1 = cell(length(lst),1);    k2 = cell(length(lst),1);
    z  = cell(length(lst),1);
    for s=1:length(lst)
      nw = length(lst{s});
      x1{s} = cell(nw,1);      x2{s} = cell(nw,1);
      w{ s} = cell(nw,1);
      k1{s} = cell(nw,1);      k2{s} = cell(nw,1);
      z{ s} = cell(nw,1);
      for I=0:nw-1
        for J=0:nw-1
          if(lst{s}(I+1)==0 & lst{s}(J+1)==0)
            x1{s}{I+1,J+1} = [];            x2{s}{I+1,J+1} = [];
            w{ s}{I+1,J+1} = [];
            k1{s}{I+1,J+1} = [];            k2{s}{I+1,J+1} = [];
            z{ s}{I+1,J+1} = [];
          else
            B = 2^(s-1);
            D = 2*B;
            Ict = I*B;            Jct = J*B;
            Imd = (I+1/2)*B;            Jmd = (J+1/2)*B;
            [t1,t2] = ndgrid([0:D-1]/D);
            x1{s}{I+1,J+1} = t1;            x2{s}{I+1,J+1} = t2;
            w{ s}{I+1,J+1} = 1/D * ones(D,D);
            k1{s}{I+1,J+1} = Imd * ones(D,D);            k2{s}{I+1,J+1} = Jmd * ones(D,D);
            z{ s}{I+1,J+1} = B * ones(D,D);
          end
        end
      end
    end
  else
    %---------------------------------------------------------
    error('wrong argument for tp');
  end
  