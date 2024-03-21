function res = mefwa2sym(x,pat,tp)
% mefwa2sym - 2D forward mirror-extended wave atom transform (symmetric
% version)
% -----------------
% INPUT
% --
% x is a real N-by-N matrix. N is a power of 2.
% --
% pat specifies the type of frequency partition which satsifies
% parabolic scaling relationship. pat can either be 'p' or 'q'.
% --
% tp is the type of tranform.
% 	'ortho': frame based on the orthobasis construction of 
% 		the standard wave atom
% 	'directional': real-valued frame with single oscillation direction
% 	'complex': complex-valued frame
% -----------------
% OUTPUT
% --
% res is an array which contains the wave atom coefficients. If
% tp=='ortho', then res is a real array of size 2N-by-2N. If
% tp=='directional', then res is a real array of size 2N-by-N-by-2. If
% tp=='complex', then res is a complex array of size 2N-by-2N.
% -----------------
% Written by Lexing Ying and Laurent Demanet, 2007

  
  if( ismember(tp, {'ortho','directional','complex'})==0 | ismember(pat, {'p','q','u'})==0 )    error('wrong');  end
  
  if(    strcmp(tp, 'ortho')==1)
    %----------------------------
    call = mefwa2(x,pat,tp);
    N = size(x,1)*2;
    res = zeros(N,N,size(call,2));
    
    for i=1:size(call,2)
      c = call(:,i);
      y = zeros(N,N);
      for s=1:length(c)
        D = 2^s;
        nw = length(c{s});
        for I=0:nw-1
          for J=0:nw-1
            if(~isempty(c{s}{I+1,J+1}))
              y( I*D+[1:D], J*D+[1:D] ) = c{s}{I+1,J+1};
            end
          end
        end
      end
      res(:,:,i) = y;
    end
  elseif(strcmp(tp,'complex')==1)
    %----------------------------
    call = mefwa2(x,pat,tp);
    N = size(x,1)*2;
    res = zeros(N,N,size(call,2));
    
    for i=1:size(call,2)
      c = call(:,i);
      y = zeros(N,N);
      for s=1:length(c)
        D = 2^s;
        nw = length(c{s});
        for I=0:nw-1
          for J=0:nw-1
            if(~isempty(c{s}{I+1,J+1}))
              y( I*D+[1:D], J*D+[1:D] ) = c{s}{I+1,J+1};
            end
          end
        end
      end
      res(:,:,i) = y;
    end
    
  elseif(strcmp(tp,'directional')==1)
    %----------------------------
    call = mefwa2(x,pat,tp);
    N = size(x,1)*2;
    res = zeros(N,N/2,size(call,2));
    
    for i=1:size(call,2)
      c = call(:,i);
      y = zeros(N,N/2);
      for s=1:length(c)
        D = 2^s;
        nw = length(c{s});
        for I=0:nw-1
          for J=0:nw-1
            if(~isempty(c{s}{I+1,J+1}))
              y( I*D+[1:D], J*D/2+[1:D/2] ) = c{s}{I+1,J+1};
            end
          end
        end
      end
      res(:,:,i) = y;
    end
  else
    error('wrong');
  end
  