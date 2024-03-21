function res = fwa1sym(x,pat,tp)
% fwa1sym - forward wave atom transform (symmetric version)
% -----------------
% INPUT
% --
% x is a real N-by-1 vector. N is a power of 2.
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
% res is an array containing all the wave atom coefficients. If
% tp=='ortho', then res is of size N-by-1. If tp=='complex', then res is 
% of size N-by-2.
% -----------------
% Written by Lexing Ying and Laurent Demanet, 2007

if( ismember(tp, {'ortho','directional','complex'})==0 | ismember(pat, {'p','q','u'})==0 )    error('wrong');  end
  x = x(:);
  
  call = fwa1(x,pat,tp);
  
  N = numel(x);
  res = zeros(N,size(call,2));
  for i=1:size(call,2)
    c = call(:,i);
    y = zeros(N,1);
    for s=1:length(c)
      D = 2^s;
      nw = length(c{s});
      for I=0:nw-1
        if(~isempty(c{s}{I+1}))
          y( I*D+[1:D] ) = c{s}{I+1};
        end
      end
    end
    res(:,i) = y;
  end
  