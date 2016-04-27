function [x, flag, err, iter] = mymr(A,b,tol,maxit,x0)
err = zeros(maxit,1);
k = 1;
r = b - A*x0;
x = x0;
err(k) = norm(r) / norm(b); % relative residual
B = A'*A;
while (err(k) > tol && k <= maxit)
    alpha = r'*A*r/(r'*B*r);
    x = x + alpha*r;
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