clear all; close all;
load west0479;
A = west0479;
n = size(A,1);
lam = eig(full(A));

% plot all the exact eigenvalues
figure(1);
plot(real(lam),imag(lam),'r+');
% title('Exact eigenvalues of west0479 computed by eig');

v = randn(n,1);
m = 60;
reorth = 1;
sh_in = 1;
tau = -10;
[V,H] = myarnoldi(A,v,m,reorth,sh_in,tau);
res1 = norm(A*V(:,1:m)-V*H);
res2 = norm(eye(m+1)-V'*V);
disp([res1,res2]);

% plot the Ritz values
hold on;
ritzv = eig(H(1:m,:));
if sh_in == 1
    lam2 = 1./ritzv + tau;
    plot(real(lam2),imag(lam2),'bo');
else
    plot(real(ritzv),imag(ritzv),'bo');
end
legend('exact eigenvalues','Ritz values');
% title('Eigenvalues and Ritz values with 60-step Arnoldi (no reorth.)');
title('Eigenvalues by shift and inverse (tau=-10)');