clear all
clc
%% Incremental Rank Power Factorization
m = 943; n = 1682; f=10;
R0=load ('943x1682.txt');
M = R0>0;
Y = M.*R0;
IDX = find(M);

S.type = '()';
S.subs{:} = IDX;

A = @(X) subsref(X,S);
Ah = @(X) subsasgn(zeros(m,n),S,X);
AhA = @(X) X.*M;

XRecon = irpf_operator_cg(A, Ah, AhA, Y(IDX), [m,n], f,f+1);
norm(XRecon-R0)/norm(R0)