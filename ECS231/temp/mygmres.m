function [x, flag, err, iter] = mygmres(A,b,m,tol,maxit,x0)
err = zeros(maxit,1);
k = 1;
r = b - A*x0;
x = x0;
err(k) = norm(r) / norm(b); % relative residual
while (err(k) > tol && k <= maxit)
    beta = norm(r);
    v = r / beta;
    [V,H] = myarnoldi(A,v,m); % call Arnoldi process
    e1 = zeros(size(H,1),1);
    e1(1) = 1;
    y = H \ (beta*e1);
    x = x + V(:,1:m)*y;
    r = b - A*x;
    k = k + 1;
    err(k) = norm(r) / norm(b); % relative residual
end
iter = k;
if (iter < maxit)
    flag = 0; % converge to tol
else
    flag = 1; % not converge but stop at maxit
end
