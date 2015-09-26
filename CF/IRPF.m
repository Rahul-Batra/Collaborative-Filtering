clear all
clc
%% Incremental Rank Power Factorization
m = 943;
n = 1682;
X=load('943x1682.txt');
L = 10;
fract = 3/20;
p = round(m*n*fract);
samplingMask = zeros(m,n);
P = randperm(m*n);
samplingMask(P(1:p)) = 1;
S.type = '()';
S.subs{:} = find(samplingMask);
A = @(X) subsref(X,S);
Ah = @(X) subsasgn(zeros(m,n),S,X);
AhA = @(X) X.*samplingMask;
b =A(X)+randn(p,1)/50;
Xrecon = irpf_operator_cg(A, Ah, AhA, b, [m,n], L,L+1);
norm(X-Xrecon)/norm(X)