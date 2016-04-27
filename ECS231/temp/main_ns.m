clear all;close all;
%% Test matrices
%[A, rows, cols, entries] = mmread('west0479.mtx');
[A, rows, cols, entries] = mmread('mahindas.mtx');

%% Input parameters: let x_sol = 1, form b = A* x_sol
x_sol = ones(rows,1);
b = A * x_sol;
maxit = 1e2;
% tol = 1e-4;
% set tol to be 1e-8 for mahindas.mtx
tol = 1e-8; 
x0 = rand(rows,1);
% m = [20,40,60,100,479]; % west0479
m = [50,100,200,400,1258]; % mahindas

%% Call mymr and plot
[x_mr, flag_mr, err_mr, iter_mr] = mymr(A,b,tol,maxit,x0);
figure(1);
loglog(1:iter_mr,err_mr(1:iter_mr),'b-');

%% Call mygmres and plot
hold on;
[x_gmres, flag_gmres, err_gmres, iter_gmres] = mygmres(A,b,m(1),tol,maxit,x0);    
loglog(1:iter_gmres,err_gmres(1:iter_gmres),'r-');
[x_gmres, flag_gmres, err_gmres, iter_gmres] = mygmres(A,b,m(2),tol,maxit,x0);    
loglog(1:iter_gmres,err_gmres(1:iter_gmres),'g-');
[x_gmres, flag_gmres, err_gmres, iter_gmres] = mygmres(A,b,m(3),tol,maxit,x0);    
loglog(1:iter_gmres,err_gmres(1:iter_gmres),'k-');
[x_gmres, flag_gmres, err_gmres, iter_gmres] = mygmres(A,b,m(4),tol,maxit,x0);    
loglog(1:iter_gmres,err_gmres(1:iter_gmres),'m-');
[x_gmres, flag_gmres, err_gmres, iter_gmres] = mygmres(A,b,m(5),tol,maxit,x0);    
loglog(1:iter_gmres,err_gmres(1:iter_gmres),'c-');
hold off;
ylabel('relative residual');
xlabel('number of iterations');
% title('MR and GMRES on west0479');
% legend('mymr','mygmres (m=20)','mygmres(m=40)',...
%    'mygmres(m=60)','mygmres(m=100)','mygmres(m=1258)');
title('MR and GMRES on mahindas');
legend('mymr','mygmres (m=50)','mygmres(m=100)',...
   'mygmres(m=200)','mygmres(m=400)','mygmres(m=1258)');