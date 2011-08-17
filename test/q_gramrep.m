function [R] = q_gramrep(V, ignoreSign)
  if ( ~exist('ignoreSign', 'var') ) ignoreSign = 0; end;
  Q = q_mat(V);
  G = q_gram(V,Q,ignoreSign);
  R = q_rep(G);
return;
