function [V, H] = myarnoldi(A,v1,m)
n = size(A,1);
V = zeros(n,m+1);
H = zeros(m+1,m);
V(:,1) = v1;
for j = 1:m
    w = A*V(:,j);
    for i = 1:j
        H(i,j) = V(:,i)'*w;
        w = w - H(i,j)*V(:,i);
    end
    H(j+1,j) = norm(w);
    if (H(j+1,j) < 1e-6 ) % determine if H(j+1,j) == 0
        disp('h(j+1,j) = 0');
        break;
    end
    V(:,j+1) = w / H(j+1,j);
end