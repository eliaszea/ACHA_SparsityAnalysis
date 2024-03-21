function x = meiwa2sym(res,pat,tp)
% meiwa2sym - 2D inverse mirror-extended wave atom transform (symmetric
% version)
% -----------------
% INPUT
% --
% res is an array which contains the wave atom coefficients. If
% tp=='ortho', then res is a real array of size 2N-by-2N. If
% tp=='directional', then res is a real array of size 2N-by-N-by-2. If
% tp=='complex', then res is a complex array of size 2N-by-2N.
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
% x is a real N-by-N matrix. N is a power of 2.
% -----------------
% Written by Lexing Ying and Laurent Demanet, 2007
  
  if( ismember(tp, {'ortho','directional','complex'})==0 | ismember(pat, {'p','q','u'})==0 )    error('wrong');  end
  
  if(    strcmp(tp,'ortho')==1)
    %-----------------------------------
    A = size(res,1);
    N  = A/2;
    lst = freq_pat(N,pat);
    red = size(res,3);
    
    call = cell(length(lst),red);
    for i=1:red
      y = res(:,:,i);
      c = cell(length(lst),1);
      for s=1:length(lst)
        B = 2^(s-1);    D = 2*B;
        nw = length(lst{s});
        c{s} = cell(nw,nw);
        for I=0:nw-1
          for J=0:nw-1
            if(lst{s}(I+1)==0 & lst{s}(J+1)==0)
              c{s}{I+1,J+1} = [];
            else
              c{s}{I+1,J+1} = y( I*D+[1:D], J*D+[1:D] );
            end
          end
        end
      end
      call(:,i) = c;
    end
    x = meiwa2(call,pat,tp);
  elseif(strcmp(tp,'complex')==1)
    %-----------------------------------
    A = size(res,1);
    N  = A/2;
    lst = freq_pat(N,pat);
    red = size(res,3);
    
    call = cell(length(lst),red);
    for i=1:red
      y = res(:,:,i);
      c = cell(length(lst),1);
      for s=1:length(lst)
        B = 2^(s-1);    D = 2*B;
        nw = length(lst{s});
        c{s} = cell(nw,nw);
        for I=0:nw-1
          for J=0:nw-1
            if(lst{s}(I+1)==0 & lst{s}(J+1)==0)
              c{s}{I+1,J+1} = [];
            else
              c{s}{I+1,J+1} = y( I*D+[1:D], J*D+[1:D] );
            end
          end
        end
      end
      call(:,i) = c;
    end
    x = meiwa2(call,pat,tp);
  elseif(strcmp(tp,'directional')==1)
    %----------------------------
    A = size(res,1);
    N = A/2;
    lst = freq_pat(N,pat);
    red = size(res,3);
    
    call = cell(length(lst),red);
    for i=1:red
      y = res(:,:,i);
      c = cell(length(lst),1);
      for s=1:length(lst)
        B = 2^(s-1);    D = 2*B;
        nw = length(lst{s});
        c{s} = cell(nw,nw);
        for I=0:nw-1
          for J=0:nw-1
            if(lst{s}(I+1)==0 & lst{s}(J+1)==0)
              c{s}{I+1,J+1} = [];
            else
              c{s}{I+1,J+1} = y( I*D+[1:D], J*D/2+[1:D/2] );
            end
          end
        end
      end
      call(:,i) = c;
    end
    x = meiwa2(call,pat,tp);
  else
    %----------------------------
    error('wrong');
  end
  
  