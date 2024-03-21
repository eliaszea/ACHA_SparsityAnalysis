function res = fwa3sym(x,pat,tp)
% fwa3sym - 3D forward wave atom transform (symmetric version)
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
% res is an array containing all the wave atom coefficients. If
% tp=='ortho', then res is of size N-by-N-by-N. If tp==' If
% tp=='directional', then res is of size N-by-N-by-N-by-4. If
% tp=='complex', then res is of size N-by-N-by-N-by-8.
% -----------------
% Written by Lexing Ying and Laurent Demanet, 2007
  
  if( ismember(tp, {'ortho','directional','complex'})==0 | ismember(pat, {'p','q','u'})==0 )    error('wrong');  end
  
  call = fwa3(x,pat,tp);
  N = size(x,1);
  res = zeros(N,N,N,size(call,2));
  
  for i=1:size(call,2)
    c = call(:,i);
    y = zeros(N,N,N);
    for s=1:length(c)
      D = 2^s;
      nw = length(c{s});
      for I=0:nw-1
        for J=0:nw-1
          for K=0:nw-1
            if(~isempty(c{s}{I+1,J+1,K+1}))
              y( I*D+[1:D], J*D+[1:D], K*D+[1:D] ) = c{s}{I+1,J+1,K+1};
            end
          end
        end
      end
    end
    res(:,:,:,i) = y;
  end
  
  