bas_in;

function [D] = double_mat( V )
[n,d] = size(V);
D = zeros(2*n,d);
for i = 1:n
D(2*i-1,:) = V(i,:);
D(2*i,:) = -1 * V(i,:);
endfor;
endfunction;

function [Q] = q_mat( V )
[n,d] = size(V);
Q = zeros(d);
for i = 1:n
Q += V(i,:)'*V(i,:);
endfor;
endfunction;

function [G] = q_gram_wr(V, Q, ignoreSign)
  if ( ~exist('ignoreSign', 'var') ) ignoreSign = 0; end;
  [n,d] = size(V);
  G = zeros(n);
  Qi = cholinv(Q);
  R = chol(Qi);
  W = V*R;
  for i = 1:n
    for j = 1:n
      G(i,j) = W(i,:)*W(j,:)';
      if ( ignoreSign ) G(i,j) = abs( G(i,j) ); end;
    end;
  end;
endfunction;

function [G] = q_gram(V, Q, ignoreSign)
  if ( ~exist('ignoreSign', 'var') ) ignoreSign = 0; end;
  [n,d] = size(V);
  G = zeros(n);
  Qinv = cholinv(Q);
  for i = 1:n
    Wi = V(i,:)*Qinv;
    for j = 1:n
      G(i,j) = Wi*V(j,:)';
      if ( ignoreSign ) G(i,j) = abs( G(i,j) ); end;
    end;
  end;
endfunction;

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
endfunction;

function [R] = q_gramrep(V, ignoreSign)
  if ( ~exist('ignoreSign', 'var') ) ignoreSign = 0; end;
  Q = q_mat(V);
  G = q_gram(V,Q,ignoreSign);
  R = q_rep(G);
endfunction;
