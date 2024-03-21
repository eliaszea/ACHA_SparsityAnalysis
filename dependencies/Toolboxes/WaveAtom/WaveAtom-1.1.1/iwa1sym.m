function x = iwa1sym(res,pat,tp)
% iwa1sym - inverse wave atom transform (symmetric version)
% -----------------
% INPUT
% --
% res is an array containing all the wave atom coefficients. If
% tp=='ortho', then res is of size N-by-1. If tp=='complex', then res is
% of size N-by-2.
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
  
  N = length(res(:,1));
  H = N/2;
  lst = freq_pat(H,pat);
  red = size(res,2);
  
  call = cell(length(lst),red);
  for i=1:red
    y = res(:,i);
    c = cell(length(lst),1);
    for s=1:length(lst)
      B = 2^(s-1);    D = 2*B;
      nw = length(lst{s});
      c{s} = cell(nw,1);
      for I=0:nw-1
        if(lst{s}(I+1)==0)
          c{s}{I+1} = [];
        else
          c{s}{I+1} = y( I*D+[1:D] );
        end
      end
    end
    call(:,i) = c;
  end
  x = iwa1(call,pat,tp);
  
  