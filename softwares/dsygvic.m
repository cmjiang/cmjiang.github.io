function [k, X, lambda] = dsygvic(n, A, B, etol)
%
% Last update -- 4/2015
% Fix-Heiberger reduction algorithm for 
% solving the generalized eigenvalue problem:
%   A*X = lambda*B*X
% where:
%   A is an n by n symmetric matrix;
%   B is an n by n positive semi-definite matrix with 
%   respect to e.
%
% INPUT:
%   A:      n by n symmetric matrix.
%   B:      n by n positive semi-definite matrix.
%   etol:   a prescribed parameter to define zero in 
%           finite precision arithmetics
%
% OUTPUT:
%   k:      the number of tol-stable eigenpairs.
%   X:      eigenvectors, dimension n by k.
%   lambda: eigenvalues, diagonal matrix, dimension k by k.
%
%------------------------------------------------------------------
%
%   Check whether A and B are both symmetric.
%
% if (issymmetric(A)==0 || issymmetric(B) == 0)
%    fprintf('A or B must be symmetric.\n');
%    return;
% end
%
%------------------------------------------------------------------
%                            Phase 1.
%------------------------------------------------------------------
%
%   Diagonalize B and sort eigvalues in the descending order.
% 
[Q_1,B0] = eig(B);
dB0 = diag(B0);
[dsB0,P0] = sort(dB0,'descend');
B0 = diag(dsB0);
%
%   Check B is semi-positive definite.
%
if ( B0(n,n) < -etol*B0(1,1) )
    fprintf('B is not semi-positive definite.\n');
    lambda = [];
    X = [];
    k = 0;
    return;
end
%
%   Find n1 and n2.
%
if( B0(1,1) <= eps(1) )
    n1 = 0;
    n2 = n;
else
    i = 1;
    while ( i < n+1 && B0(i,i) > etol*B0(1,1) )
        i = i + 1;
    end
    n1 = i-1;
    n2 = n-n1;
end
%
%   Early exit: n1 = 0.
%
if ( n1 == 0 )
    if ( det(A) == 0 )
        fprintf('singular pencil, exit.#1\n');
        k = -1;
        X = [];
        lambda = [];
        return;
    else
        fprintf('regular pencil but no finite eigenvalue, exit.#1\n');
        k = 0;
        X = [];
        lambda = [];
        return;
    end    
end
%
%   Update A
%
Q_1 = Q_1(:,P0);
A0 = Q_1'*A*Q_1;
%
%   Reduce B0 to identity + zero.
%
D0 = B0(1:n1,1:n1);
dD0 = diag(D0);
drD0 = 1./sqrt(dD0);
R_1(1:n1,1:n1) = diag(drD0);
R_1(n1+1:n,n1+1:n) = eye(n2);
%
%   Update A0 (and B0)
%
A1 = R_1'*A0*R_1;
%
%   Early exit: n2 = 0. B is a e-well-conditioned matrix.
%
if (n2 == 0)
    fprintf('%d e-finite eigenvalues.#1\n',n);
    [U, lambda] = eig(A1);
    X = Q_1*R_1*U;
    k = n;
    return;
end
%------------------------------------------------------------------
%                            Phase 2.
%------------------------------------------------------------------
%
%   Diagonalize A1_22 and sort eigenvalue in descending order
%
A1_22 = A1(n1+1:n,n1+1:n);
A1_22 = 0.5*(A1_22'+A1_22);
[Q2_22,D02] = eig(A1_22);
[~, P2] = sort(abs(diag(D02)),'descend');
D2 = D02(P2,P2);
%
%   Find n3 and n4
%
if( abs(D2(1,1))<=etol )
    n3 = 0;
    n4 = n2;
else
    i = 1;
    while (i < n2+1 && abs(D2(i,i)) > etol*abs(D2(1,1)))
        i = i + 1;
    end
    n3 = i-1;
    n4 = n2-n3;
end
%
%   Early exit: n3 = 0.
%
if( n3 == 0 )
%
%   If n1 < n2, A-lambda*B is singluar, exit
%
    if ( n1 < n2 )
        fprintf('singular pencil, exit.#2\n');
        k = -1;
        X = [];
        lambda = [];
        return;
    end
%
%   If n1 >= n2
%       reveal the rank of A1_12 by QR with column pivoting
%
    A1_12 = A1(1:n1,n1+1:n);
    [Q2_12,R2_12,P2_12] = qr(A1_12);
%
%   If n1 = n2:
%       if R1_12 rank deficient, then singular pencil, exit;
%       else regular pencil but no e-finite eigenvalues, exit.
%
    if ( n1 == n2 )
        if( abs(R2_12(n2,n2)) <= eps(1))
            fprintf('singular pencil, exit.3\n');
            k = -1;
            X = [];
            lambda = [];
            return;
        else
            fprintf('regular pencil');
            fprintf('but no e-finite eigenvalues,exit.2\n');
            k = 0;
            X = [];
            lambda = [];
            return;
        end
    end
%
%   If n1 > n2:
%       if R2_12 rank deficient, then singular pencil, exit;
%       else n1-n2 e-finite eigenvalues.
%
    if( abs(R2_12(n2,n2)) <= eps(1) )
        fprintf('singular pencil, exit.4\n');
        k = -1;
        X = [];
        lambda = [];
        return;
    else
        fprintf('%d e-finite eigenvalues.#2\n',n1-n2);
        Q_2(1:n1,1:n1) = Q2_12;
        Q_2(n1+1:n,n1+1:n) = P2_12;
        A2 = Q_2'*A1*Q_2;
        A2_22 = A2(n2+1:n1,n2+1:n1);
        A2_13 = A2(1:n2,n1+1:n);
        A2_12 = A2(1:n2,n2+1:n1);
        [U2, lambda] = eig(A2_22);
        U3 = -A2_13\A2_12*U2;
        U = [zeros(n2,n1-n2);U2;U3];
        X = Q_1*R_1*Q_2*U;
        k = n1-n2;
        return;
    end
end
%
%   Update A1 (and B1)
%
Q_2(1:n1,1:n1) = eye(n1);
Q_2(n1+1:n,n1+1:n) = Q2_22(:,P2);
A2 = Q_2'*A1*Q_2;
%
%   Early exit: n4 = 0.
%
if (n4 == 0)
    fprintf('%d e-finite eigenvalues.#3\n',n1);
    A2_11 = A2(1:n1,1:n1);
    A2_12 = A2(1:n1,n1+1:n);
    A2_U1 = A2_11-A2_12/D2*A2_12';
    [U1, lambda] = eig(A2_U1);
    U2 = -D2\A2_12'*U1;
    U = [U1;U2];
    X = Q_1*R_1*Q_2*U;
    k = n1;
    return;
end
%
%------------------------------------------------------------------
%                            Phase 3.
%------------------------------------------------------------------
%
%   Since n3 ~= 0 and n4 ~=0, we can view A2 as a 3 by 3 block:
%
%               n1     n3     n4
%        n1 | A2_11  A2_12  A2_13 |
%   A2 = n3 | A2_12'   D2     0   |
%        n4 | A2_13'   0      0   |
%
%   Early exit: n1 < n4.
%
if ( n1 < n4 )
    fprintf('singular pencil, exit. 5\n');
    k = -1;
    X = [];
    lambda = [];
    return;
end
%
%   If n1 >= n4
%       reveal the rank of A2_13 by QR with column pivoting
%
A2_13 = A2(1:n1,n1+n3+1:n);
[Q3_11,R3_13,P3_13] = qr(A2_13);
%
%   n1 = n4:
%       if R3_13 rank deficient, singular pencil,
%       else regular pencil but no e-finite eigenvalues.
%
if ( n1 == n4 )
    if( abs(R3_13(n4,n4)) <= eps(1))
        fprintf('singular pencil, exit.6\n');
        k = -1;
        X = [];
        lambda = [];
        return;
    else
        fprintf('regular pencil');
        fprintf('but no e-finite eigenvalues,exit.3\n');
        k = 0;
        X = [];
        lambda = [];
        return;
    end
end
%
%   n1 > n4:
%       if R3_13 rank deficient, singular pencil,
%       else n1-n4 e-finite eigenvalues.
%
if( abs(R3_13(n4,n4)) <= eps(1))
    fprintf('singular pencil, exit.#7\n');
    k = -1;
    X = [];
    lambda = [];
    return;
else
    fprintf('%d e-finite eigenvalues.#4\n',n1-n4);
    n5 = n1 - n4;
    Q_3(1:n1,1:n1) = Q3_11;
    Q_3(n1+1:n1+n3,n1+1:n1+n3) = eye(n3);
    Q_3(n1+n3+1:n,n1+n3+1:n) = P3_13;
%
%   Update A2, A3 is a 4 by 4 block matrix:
%
%               n1     n5     n3    n4
%        n1 | A3_11  A3_12  A3_13  A3_14 |
%   A3 = n5 | A3_12' A3_22  A3_23    0   |
%        n3 | A3_13' A3_23' D2_33    0   |
%        n4 | A3_14'   0      0      0   |
%
    A3 = Q_3'*A2*Q_3;
    A3_12 = A3(1:n4,n4+1:n4+n5);
    A3_13 = A3(1:n4,n4+n5+1:n4+n5+n3);
    A3_14 = A3(1:n4,n4+n5+n3+1:n);
    A3_22 = A3(n4+1:n4+n5,n4+1:n4+n5);
    A3_23 = A3(n4+1:n4+n5,n4+n5+1:n4+n5+n3);
    D2_33 = A2(n1+1:n1+n3,n1+1:n1+n3);
    A3_U2 = A3_22-A3_23/D2_33*A3_23';
    A3_U2 = 0.5*(A3_U2+A3_U2');
    [U2, lambda] = eig(A3_U2);
    U3 = -D2_33\(A3_23'*U2);
    U4 = -A3_14\(A3_12*U2+A3_13*U3);
    k = n5;
    U = [zeros(n4,k);U2;U3;U4];
    X = Q_1*R_1*Q_2*Q_3*U;
    return;
end
