clear all;close all;
%% Test matrices
%[A, rows, cols, entries] = mmread('west0479.mtx');
[A, rows, cols, entries] = mmread('mahindas.mtx');

%% Input parameters: let x_sol = 1, form b = A* x_sol
x_sol = ones(rows,1);
b = A * x_sol;
maxit = 1e1;
tol = 1e-8; 
m=60;

%% Call mymr, mygmres(m=60) and Matlab function gmres with incomplete LU
[x1,fl1,rr1,it1,rv1] = mymr(A,b,tol,maxit);
[x2,fl2,rr2,it2,rv2] = mygmres(A,b,m,tol,maxit);
[L,U] = ilu(A,struct('type','ilutp','droptol',tol));
[x3,fl3,rr3,it3,rv3] = gmres(A,b,[],tol,maxit,L,U);

%% Plot
figure(1);
semilogy(0:it1,rv1/norm(b),'b-','LineWidth',2);
hold on;
semilogy(0:it2,rv2/norm(b),'r-','LineWidth',2);
semilogy(0:it3(2),rv3/norm(b),'g-','LineWidth',2);
hold off;
ylabel('relative residual');
xlabel('number of iterations');
%title('MR, GMRES and GMRES (incomplete LU) on west0479');
title('MR, GMRES and GMRES (incomplete LU) on mahindas');
legend('mymr','mygmres (m=60)', 'gmres (incomplete LU)');