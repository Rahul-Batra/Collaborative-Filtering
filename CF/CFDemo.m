clear all
clc

m = 943; n = 1682; f=10;
%matrix=load ('943x1682.txt');
load('matrix.mat');
M = matrix>0;
Y = M.*matrix;
IDX = find(M);

S.type = '()';
S.subs{:} = IDX;

A = @(X) subsref(X,S);
Ah = @(X) subsasgn(zeros(m,n),S,X);
AhA = @(X) X.*M;
%% Non-Negative Matrix Factorization
[P,Q] = nnmf(Y, 20);
XRecon = P*Q;
norm(matrix - XRecon)/norm(matrix)

%% Incremental Rank Power Factorization
XRecon = irpf_operator_cg(A, Ah, AhA, Y(IDX), [m,n], f,f+1);
norm(XRecon-matrix)/norm(matrix)
%% Singular Value Thresholding
tau = 5*sqrt(m*n);
delta = 1.2*length(IDX)/(m*n);
[U,S,V,numiter] = SVT([m n],IDX,matrix(IDX),tau,delta);
XRecon = U*S*V';
norm(matrix-XRecon,'fro')/norm(matrix,'fro')
%% Fixed Point Continuation
mu_final = 0.1;
[U,S,V,numiter] = FPC([m, n],IDX,matrix(IDX),mu_final);
XRecon = U*S*V';
norm(matrix-XRecon,'fro')/norm(matrix,'fro')
