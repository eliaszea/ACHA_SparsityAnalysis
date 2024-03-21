function [x1,w,k1,z] = pwa1(N,pat,tp)
% pwa1 - get position information
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
% x1 is a cell structure that contains the centers of the wave atoms
% in spatial domain.
% --
% w is a cell structure that contains the widths of the wave atoms
% in spatial domain.
% --
% k1 is a cell structure that contains the centers of the wave atoms
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
    x1 = cell(length(lst),1);
    w  = cell(length(lst),1);
    k1 = cell(length(lst),1);
    z  = cell(length(lst),1);
    for s=1:length(lst)
      nw = length(lst{s});
      x1{s} = cell(nw,1);
      w{ s} = cell(nw,1);
      k1{s} = cell(nw,1);
      z{ s} = cell(nw,1);
      for I=0:nw-1
        if(lst{s}(I+1)==0)
          x1{s}{I+1} = [];
          w{ s}{I+1} = [];
          k1{s}{I+1} = [];
          z{ s}{I+1} = [];
        else
          B = 2^(s-1);
          D = 2*B;
          Ict = I*B;
          Imd = (I+1/2)*B;
          x1{s}{I+1} = [0:D-1]'/D;
          w{ s}{I+1} = 1/D * ones(1,D);
          k1{s}{I+1} = Imd * ones(D,1);
          z{ s}{I+1} = B * ones(1,D);
        end
      end
    end
  else
    %---------------------------------------------------------
    error('wrong argument for tp');
  end
  