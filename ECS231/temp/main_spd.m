clear all;close all;
%% Test matrices
%[A, rows, cols, entries] = mmread('bcsstk15.mtx');
[A, rows, cols, entries] = mmread('NOS3.mtx');

%% Input parameters: let x_sol = 1, form b = A* x_sol
x_sol = ones(rows,1);
b = A * x_sol;
maxit = 1e5;
tol = 1e-4;
x0 = rand(rows,1);

%% Call mysd, mycg and Matlab function pcg
[x_sd, flag_sd, err_sd, iter_sd] = mysd(A,b,tol,maxit,x0);
[x_cg, flag_cd, err_cg, iter_cg] = mycg(A,b,tol,maxit,x0);
[x_pcg, flag_pcg, res_pcg, iter_pcg, err_pcg] = pcg(A,b,tol,maxit);

%% Plot data
figure(1);
loglog(1:iter_sd,err_sd(1:iter_sd),'b-');
hold on;
loglog(1:iter_cg,err_cg(1:iter_cg),'r-');
loglog(1:iter_pcg,err_pcg(1:iter_pcg)/norm(b),'g-');
hold off;
ylabel('relative residual');
xlabel('number of iterations');
title('SD, CG and PCG on NOS3');
legend('mysd','mycg','pcg');