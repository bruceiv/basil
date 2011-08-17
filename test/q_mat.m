function [Q] = q_mat( V )
[n,d] = size(V);
Q = zeros(d);
for i = 1:n
Q += V(i,:)'*V(i,:);
endfor;
return;
