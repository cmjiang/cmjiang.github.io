clear all;close all;
%% Test matrices
%[A, rows, cols, entries] = mmread('bcsstk15.mtx');
[A, rows, cols, entries] = mmread('NOS3.mtx');

%% Input parameters: let x_sol = 1, form b = A* x_sol
x_sol = ones(rows,1);
b = A * x_sol;
maxit = 1e4;
tol = 1e-4;
% set tol to be 1e-8 for mahindas.mtx
% tol = 1e-8; 
x0 = rand(rows,1);
m = 20;

%% Call mymr and plot
[x_mr, flag_mr, err_mr, iter_mr] = mymr(A,b,tol,maxit,x0);
figure(1);
loglog(1:iter_mr,err_mr(1:iter_mr),'b-');
hold on;


%% Call mygmres and plot
[x_gmres, flag_gmres, err_gmres, iter_gmres] = mygmres(A,b,m,tol,maxit,x0);    
loglog(1:iter_gmres,err_gmres(1:iter_gmres),'r-');
hold off;
ylabel('relative residual');
xlabel('number of iterations');
title('MR and GMRES on NOS3');
legend('mymr','mygmres (m=20)');