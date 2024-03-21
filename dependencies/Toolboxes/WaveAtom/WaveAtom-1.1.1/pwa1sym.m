function [a1,ww,b1,zz] = pwa1sym(N,pat,tp)
% pwa1sym - get position information
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
% a1 is an array that contains the centers of the wave atoms
% in spatial domain.
% --
% ww is an array that contains the widths of the wave atoms
% in spatial domain.
% --
% b1 is an array that contains the centers of the wave atoms
% in frequency domain.
% --
% zz is an array that contains the widths of the wave atoms
% in frequency domain.
% --
% -----------------
% Written by Lexing Ying and Laurent Demanet, 2007
  
  if( ismember(tp, {'ortho','directional','complex'})==0 | ismember(pat, {'p','q','u'})==0 )    error('wrong');  end
  
  [x1,w,k1,z] = pwa1(N,pat,tp);
  a1 = zeros(N,1);
  ww = zeros(N,1);
  b1 = zeros(N,1);
  zz = zeros(N,1);
  
  for s=1:length(x1)
    D = 2^s;
    nw = length(x1{s});
    for I=0:nw-1
      if(~isempty(x1{s}{I+1}))
        a1( I*D+[1:D] ) = x1{s}{I+1};
        ww( I*D+[1:D] ) = w{ s}{I+1};
        
        b1( I*D+[1:D] ) = k1{s}{I+1};
        zz( I*D+[1:D] ) = z{ s}{I+1};
      end
    end
  end
  