clear all
clc
%% Singular Value Thresholding
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

tau = 5*sqrt(m*n);
delta = 1.2*length(IDX)/(m*n);
[U,S,V,numiter] = SVT([m n],IDX,R0(IDX),tau,delta);
XRecon = U*S*V';
norm(R0-XRecon,'fro')/norm(R0,'fro')