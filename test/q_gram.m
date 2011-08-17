function [G] = q_gram(V, Q, ignoreSign)
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
return;
