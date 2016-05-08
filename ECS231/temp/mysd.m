function [x,flag,relres,iter,resvec] = mysd(A,b,tol,maxit)
% mysd attempts to solve the system of linear equations
% A*x = b
% for x.
% Inputs:
%   A: n by n nonsingular matrix;
%   b: n by 1 vector;
%   tol: tolrence of the method;
%   maxit: maximum number of iterations.
% Outputs:
%   x: n by 1 output vector;
%   flag: convergence flag --
%         0-converged to tol within maxit;
%         1-iterated maxit times but did not converge;
%   relres: the relative residual norm(r)/norm(b);
%   iter: the iteration number at which x is computed;
%   resvec: a vector of the residual norms at each iteration including norm(b-A*x0).

n = size(A,1);
x = zeros(n,1);
resvec = zeros(maxit,1);
iter = 1;
r = b - A*x;
resvec(iter) = norm(r);
relres = resvec(iter) / norm(b); % relative residual 
while (relres > tol && iter <= maxit)
    alpha = r'*r / (r'*A*r);
    x = x + alpha*r;
    r = b - A*x;
    iter = iter + 1;
    resvec(iter) = norm(r);
    relres = resvec(iter)/ norm(b); % relative residual
end
if (iter <= maxit)
    flag = 0; % converge to tol
else
    flag = 1; % not converge but stop at maxit
end
resvec = resvec(1:iter);
iter = iter - 1;