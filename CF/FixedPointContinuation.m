clear all
clc
%% Fixed Point Continuation
m = 943; n = 1682; 
R0=load ('943x1682.txt');
M = R0>0;
Y = M.*R0;
IDX = find(M);

S.type = '()';
S.subs{:} = IDX;

A = @(X) subsref(X,S);
Ah = @(X) subsasgn(zeros(m,n),S,X);
AhA = @(X) X.*M;

mu_final = 2.502370;
[U,S,V,numiter] = FPC([m, n],IDX,R0(IDX),mu_final);
XRecon = U*S*V';
norm(R0-XRecon,'fro')/norm(R0,'fro')
