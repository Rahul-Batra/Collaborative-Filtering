clear all
clc
%% Non-Negative Matrix Factorization
m = 943; n = 1682; 
%R0=load ('943x1682.txt');
load('matrix.mat');
M = matrix>0;
Y = M.*matrix;
IDX = find(M);

S.type = '()';
S.subs{:} = IDX;

A = @(X) subsref(X,S);
Ah = @(X) subsasgn(zeros(m,n),S,X);
AhA = @(X) X.*M;

[P,Q] = nnmf(Y, 20);
XRecon = P*Q;
norm(matrix - XRecon)/norm(matrix)
