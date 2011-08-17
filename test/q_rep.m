function [R,H] = q_rep(G,ignoreSign)
  if ( ~exist('ignoreSign', 'var') ) ignoreSign = 0; end;
  [n,d] = size(G);
  R = zeros(n,d);
  s = strrep(num2str(G(1,1), 'v%.10g'), '.', '_');
  if ( ignoreSign ) 
    s = strrep(s, '-', ''); 
  else 
    s = strrep(s, '-', 'n'); 
  end;
  H.(s) = 0; % load value into hashtable
  nextVal = 1;
  for i = 1:n
    for j = 1:d
      s = strrep(num2str(G(i,j), 'v%.10g'), '.', '_');
      if ( ignoreSign ) 
        s = strrep(s, '-', ''); 
      else 
        s = strrep(s, '-', 'n'); 
      end;
      
      if ( ~isfield(H, s) )
        H.(s) = nextVal;
	nextVal = nextVal + 1;
      end;
      R(i,j) = H.(s);
    end;
  end;
return;
