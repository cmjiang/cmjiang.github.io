clear all;close all;
%% Test matrices
[A, rows, cols, entries] = mmread('bcsstk15.mtx');
%[A, rows, cols, entries] = mmread('NOS3.mtx');

%% Input parameters: let x_sol = 1, form b = A* x_sol
x_sol = ones(rows,1);
b = A * x_sol;
maxit = 2e2;
tol = 1e-6;
dopcgp = 0; % set to 0 for bcsstk15 and 1 for NOS3

%% Call mysd, mycg and Matlab function pcg without/with a precondioner
[x1,fl1,rr1,it1,rv1] = mysd(A,b,tol,maxit);
[x2,fl2,rr2,it2,rv2] = mycg(A,b,tol,maxit);
[x3,fl3,rr3,it3,rv3] = pcg(A,b,tol,maxit);
if (dopcgp)
    L = ichol(A);
    [x4,fl4,rr4,it4,rv4] = pcg(A,b,tol,maxit,L,L');
end

%% Plot data
figure(1);
semilogy(0:it1,rv1/norm(b),'b-','LineWidth',2);
hold on;
semilogy(0:it2,rv2/norm(b),'r-','LineWidth',2);
semilogy(0:it3,rv3(1:it3+1)/norm(b),'g:','LineWidth',2);
if (dopcgp)
    semilogy(0:it4,rv4/norm(b),'k-','LineWidth',2);
end
hold off;
ylabel('relative residual');
xlabel('number of iterations');
if (~dopcgp)
    title('SD, CG, PCG (no precond.) on bcsstk15');
    legend('mysd','mycg','pcg (no precond.)');
else
    title('SD, CG, PCG (no precond.) and PCG (with precond.) on NOS3');
    legend('mysd','mycg','pcg (no precond.)','pcg (with precond.)');
end