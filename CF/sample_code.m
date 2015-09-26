% This sample code shows how to use Incremented Rank PowerFactorization to
% perform rank-constrained matrix recovery for a couple different cases of
% reconstructing a 512x512 matrix from limited measurements
% (this implementation is not necessarily the fastest IRPF implementation for
% all of these cases -- other versions of IRPF available on request).
%
% Justin Haldar (jhaldar@usc.edu)
% May 21, 2012

clear;
close all;
clc;

%% Load Sample Data for Testing
disp('Loading and Displaying Example Data');
m = 943;
n = 1682;
X=load('943x1682.txt');


disp('*******************************************');
%% Matrix Completion Example
disp('Random Matrix Completion Example');
% A standard matrix completion example, with random sampling

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
tic
Xrecon = irpf_operator_cg(A, Ah, AhA, b, [m,n], L,L+1);
toc
NMAE = norm(X-Xrecon)/norm(X)
