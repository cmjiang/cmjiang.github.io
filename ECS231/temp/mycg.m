function [x, flag, err, iter] = mycg(A,b,tol,maxit,x0)
err = zeros(maxit,1);
k = 1;
r = b - A*x0;
x = x0;
p = r;
err(k) = norm(r) / norm(b); % relative residual
while (err(k) > tol && k <= maxit)
    theta = r'*r/(p'*A*p);
    x = x + theta*p;
    rr = r - theta*A*p;
    tau = rr'*rr / (r'*r);
    p = rr + tau*p;
    r = rr;
    k = k + 1;
    err(k) = norm(r) / norm(b); % relative residual
end
iter = k;
if (iter < maxit)
    flag = 0; % converge to tol
else
    flag = 1; % not converge but stop at maxit
end