m = 943; n = 1682;f=20;
R=zeros(943,1682);
fp = fopen('u.data');
A = fscanf(fp, '%d %d');
for i=1:4:400000
    R(A(i:i,1:1):A(i:i,1:1),A(i+1:i+1,1:1):A(i+1:i+1,1:1)) = A(i+2:i+2,1:1);
end

M = createSamplingScheme ([m n], 'random', 0.2);
Y = M.*R;
IDX = find(M);
S.type = '()';
S.subs{:} = IDX;

Ab = @(X) subsref(X,S);
Ah = @(X) subsasgn(zeros(m,n),S,X);
AhA = @(X) X.*M;
%% Non-Negative Matrix Factorization
[P,Q] = nnmf(R, 60);
XRecon = P*Q;
norm(R - XRecon)/norm(R)

%% Singular Value Thresholding
%tau = 5*sqrt(m*n);
%delta = 1.2*length(IDX)/(m*n);
%[U,S,V,numiter] = SVT([m n],IDX,R(IDX),tau,delta);
%XRecon = U*S*V';
%norm(R-XRecon,'fro')/norm(R,'fro')

%% Incremental Rank Power Factorization
%XRecon = irpf_operator_cg(Ab, Ah, AhA, Y, [n,m], f,f+1);
%norm(XRecon-R)/norm(R)

%% Fixed Point Continuation
%mu_final = 0.1;
% [U,S,V,numiter] = FPC([m, n],IDX,R(IDX),mu_final);
% XRecon = U*S*V';
% norm(R-XRecon,'fro')/norm(R,'fro')