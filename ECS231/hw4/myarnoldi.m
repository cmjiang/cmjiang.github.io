function [V,H] = myarnoldi(A,v,m,reorth,sh_in,tau)
% Input
%   A: input matrix
%   v: initial vector
%   m: Arnoldi steps
%   reorth: =1 reorthogonalization
%   sh_in: =1 shift-and-invert spectral transformation
%   tau: target eigenvalue used in spectral transformation
% Output
%   V: orthonormal matrix
%   H: upper Hessenberg matrix
n = size(A,1);
V = zeros(n,m+1);
H = zeros(m+1,m);
V(:,1) = v/norm(v);
if sh_in == 1
    [L, R, P] = lu(A-tau*speye(n));
end
for j=1:m
    if sh_in == 1
        w = R \(L \ (P*V(:,j)));
    else
        w = A*V(:,j);
    end
    for i=1:j
        H(i,j) = V(:,i)'*w;
        w = w - H(i,j)*V(:,i);
    end
    if reorth == 1
        for i = 1:j
            temp = (V(:,i))'*w;
            w = w - temp*V(:,i);
            H(i,j)=H(i,j)+temp;
        end
    end
    H(j+1,j) = norm(w);
    if H(j+1,j) < 1e-10
        return
    end
    V(:,j+1) = w/H(j+1,j);
end