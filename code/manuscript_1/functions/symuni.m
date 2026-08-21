function out=symuni(A)
% source: https://ieeexplore.ieee.org/document/10319662
% github: https://github.com/YijieLinaMao/BD-RIS-low-complexity
% symmetric unitary projection
% input a square matrix
% outout a symmetric unitary matrix
[U,S,V]=svd(A+A.');
R=rank(S);
N=size(A,1);
U(:,R+1:N)=conj(V(:,R+1:N));
out=U*V';

end

