function lst = freq_pat(H,pat)
% freq_pat.m - generate frequency partition
%
% input:
%   H:    half length of the frequency span, H needs to be greater than 16
%   pat:  type of frequency partition which satsifies parabolic scaling relationship
%         equal to 'p' or 'q'
%
% output:
%   lst:  data representing the partition of frequency
%
% Written by Lexing Ying and Laurent Demanet, 2006
  
  if(H<16)    error('H needs to be at least 16');  end
  
  if(    pat=='p') %parabolic scaling by lexing
    len = ceil(log2(H)/2)+1;
    lst = cell(len,1);
    lst{1} = [1 1];    lst{2} = [0 1];    lst{3} = [0 1 1 1];
    cnt = 16;
    idx = 3;    rad = 4;
    while(cnt<H)
      old = cnt;
      lst{idx} = [lst{idx} 1 1];
      cnt = cnt + 2*rad;
      idx = idx+1;
      rad = 2*rad;
      trg = min(4*old,H);
      lst{idx} = [zeros(1,cnt/rad), ones(1,(trg-cnt)/rad)];
      cnt = trg;
    end
  elseif(pat=='q') %parabolic scaling by laurent
    len = floor(log2(H)/2)+1;
    lst = cell(len,1);
    lst{1} = [1 1];    lst{2} = [0 1 1 1];
    cnt = 8;
    idx = 2;      rad = 2;
    while(cnt<H)
      old = cnt;
      lst{idx} = [lst{idx} 1 1];
      cnt = cnt + 2*rad;
      idx = idx+1;
      rad = 2*rad;
      trg = min(4*old,H);
      lst{idx} = [zeros(1,cnt/rad), ones(1,(trg-cnt)/rad)];
      cnt = trg;
    end
  elseif(pat=='u')
    len = ceil(log2(H)/2)+1; %uniform partitioning
    lst = cell(len,1);
    B = 2^(len-1);
    lst{end} = ones(1,H/B);
  else
    error('wrong pat');
  end
  
