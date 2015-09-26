% Function: [X,U,V] =  irpf_operator_cg( A, Ah, b, AhA, sizex, L, MAX_ITERS, epsilon, initU, initV )
%
% Description: Incremented Rank PowerFactorization algorithm to solve (P0), with
%              modifications to use iterative conjugate-gradient matrix
%              inversion (other variations of Incremented Rank
%              PowerFactorization can be faster than this one, depending on the
%              problem structure -- alternative implementations are available
%              on request).
%
% (P0):  Minimize || A(X) - b ||_F
%            s.t. rank(X) <= L,
%   where X is an [m X n] real or complex matrix and b is a [p X 1]
%   real or complex data vector
%
%
% Main Reference:
%
%      *J. P. Haldar and D. Hernando.  "Rank-Constrained Solutions to Linear Matrix Equations using PowerFactorization."  IEEE Signal Processing Letters 16:584-587, 2009.
%
% Supplementary references (for the conjugate-gradient component):
%
%      *J. P. Haldar and Z.-P. Liang. "Spatiotemporal imaging with partially separable functions: A matrix recovery approach."  IEEE International Symposium on Biomedical Imaging, 2010, pp. 716-719.
%      *J. P. Haldar.  "Constrained Imaging: Denoising and Sparse Sampling."  Ph.D. Dissertation, University of Illinois at Urbana-Champaign, May 2011.
%
% Arguments:
%   - A: function handle that performs the forward measurement operator A(X),
%        which maps an [m X n] matrix into a [p X 1] data vector
%   - Ah: function handle that performs the adjoint of A, and maps a [p X 1]
%         vector into an [m X n] matrix
%   - AhA: function handle that performs Ah(A(x)), and maps an [m X n]
%         matrix into an [m X n] matrix
%   - b: [p X 1] data vector
%   - sizeX: A length-2 vector, with sizeX(1) = m and sizeX(2) = n
%   - L: rank of the desired matrix
%   - MAX_ITERS (optional): maximum number of alternation iterations --
%                           more iterations can significantly improve the 
%                           quality of the results (default = 20)
%   - epsilon (optional): stopping threshold for relative error of the fit (so
%                         algorithm will stop when MAX_ITERS is reached or the
%                         relative data error (measured in the 2-norm) goes 
%                         below epsilon) (default = 1e-6)
%   - initU (optional): Initialization for the [m X r] U matrix (default = random)
%   - initV (optional): Initialization for the [r X n] V matrix (default = random)
%
% Returns:
%   - X: [m X n] the estimated matrix
%   - U: [m X r] matrix, with X = U*V
%   - V: [r X n] matrix, with X = U*V
%
% Authors: Justin Haldar (jhaldar AT usc DOT edu) and Diego Hernando
% Date created: October 26, 2008
% Date last modified: May 21, 2012
% Copyright 2012
%
% This code is licensed under the CC Attribution-Noncommercial-Share Alike 
% 3.0 License (http://creativecommons.org/licenses/by-nc-sa/3.0/).  Please 
% cite at least the main IRPF reference (and the related references if 
% appropriate) if you use this code or its derivatives in your own work.
%
function [X,U,V] = irpf_operator_cg( A,Ah, AhA, b, sizeX, L, MAX_ITERS, epsilon, initU, initV )

disp('   IRPF starting:');
if not(exist('epsilon','var'))|| numel(epsilon)==0
    epsilon = 1e-6;
end

if not(exist('MAX_ITERS','var'))||(numel(MAX_ITERS)==0)
    MAX_ITERS = 20;
end

U = zeros(sizeX(1),0);
V = zeros(0,sizeX(2));
for myR = 1:L
    disp(['     Rank ' num2str(myR) '/' num2str(L)]);
    if (not(exist('initU','var'))||(numel(initU)==0))
        initU_up = randn(sizeX(1),1);
    else
        initU_up = initU(:,myR);
    end
    if (not(exist('initV','var'))||(numel(initV)==0))
        initV_up = randn(1,sizeX(2));
    else
        initV_up = initV(myR,:);
    end
    [~,U_up,V_up] = pf_operator_cg( A, Ah, AhA, b-A(reshape(U*V,[],1)), sizeX, 1, MAX_ITERS, epsilon, initU_up, initV_up );
    [X, U, V] = pf_operator_cg( A, Ah, AhA, b, sizeX, myR, MAX_ITERS, epsilon, [U,U_up], [V;V_up] );
end
return;

%% function pf_operator_cg
function [Xfinal,x1,y1] = pf_operator_cg( A,Ah, AhA, b, sizex, r, MAX_ITERS, epsilon, initx, inity )

m = sizex(1);
n = sizex(2);

x1 = initx;
y1 = inity;

xold = x1*10;
yold = y1*10;

nrefine = 1;
curerror = epsilon+1;
while nrefine<=MAX_ITERS && curerror>epsilon && ((norm(x1-xold,'fro')+norm(y1-yold,'fro'))/(norm(x1)+norm(y1))>epsilon)
    xold = x1;
    yold = y1;
    
    x1 = reshape(cgsolver(@(x) reshape(AhA(reshape(x,[m,r])*y1)*y1',[],1), reshape(reshape(Ah(b),[m,n])*y1',[],1), 1e-6, max(r+1,10), x1(:)),[m,r]);
    y1 = reshape(cgsolver(@(y) reshape(x1'*AhA(x1*reshape(y,[r,n])),[],1), reshape(x1'*reshape(Ah(b),[m,n]),[],1), 1e-6, max(r+1,10), y1(:)),[r,n]);
    
    curerror = norm(b -A(reshape(x1*y1,[],1)))/norm(b);
    
    nrefine = nrefine + 1;
end
Xfinal = x1*y1;
return;

%% cgsolver
function out = cgsolver(A,b,epsilon,maxit,x0)
[out,~] = pcg(A,b,epsilon,maxit,[],[],x0);
return;
