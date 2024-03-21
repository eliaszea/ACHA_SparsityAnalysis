function x = iwa2sym(res,pat,tp)
% iwa3sym - 3D inverse wave atom transform (symmetric version)
% -----------------
% INPUT
% --
% res is an array containing all the wave atom coefficients. If
% tp=='ortho', then res is of size N-by-N-by-N. If tp==' If
% tp=='directional', then res is of size N-by-N-by-N-by-4. If
% tp=='complex', then res is of size N-by-N-by-N-by-8.
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
  
  N = size(res,1);
  H = N/2;
  lst = freq_pat(H,pat);
  red = size(res,4);
  
  call = cell(length(lst),red);
  for i=1:red
    y = res(:,:,:,i);
    c = cell(length(lst),1);
    for s=1:length(lst)
      B = 2^(s-1);    D = 2*B;
      nw = length(lst{s});
      c{s} = cell(nw,nw,nw);
      for I=0:nw-1
        for J=0:nw-1
          for K=0:nw-1
            if(lst{s}(I+1)==0 & lst{s}(J+1)==0 & lst{s}(K+1)==0)
              c{s}{I+1,J+1,K+1} = [];
            else
              c{s}{I+1,J+1,K+1} = y( I*D+[1:D], J*D+[1:D], K*D+[1:D] );
            end
          end
        end
      end
    end
    call(:,i) = c;
  end
  x = iwa3(call,pat,tp);
  
  