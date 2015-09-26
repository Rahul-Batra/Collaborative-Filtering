clear all
clc

m = 943; n = 1682; f = 100;

P0 = rand(m,f);
Q0 = rand(n,f);
R0 = P0*Q0';
Y=zeros(943,1682);
Y2=zeros(943,1682);
       
fileName  = 'u1.base';
inputfile = fopen(fileName);
delimiter = ' ';  % for tabs: delimiter = sprintf('\t');

% Get the values from the file
values = textscan(inputfile, '%d%d%d%d', 'delimiter', delimiter);

fclose(inputfile); 
h= length(values{1});
for j=1:h
    user = values{1}(j);
    item = values{2}(j);
    rating = values{3}(j);
    Y((item-1)*943+user)=rating;
    
end


fileName  = 'u1.test';
inputfile = fopen(fileName);
delimiter = ' ';  % for tabs: delimiter = sprintf('\t');

% Get the values from the file
values = textscan(inputfile, '%d%d%d%d', 'delimiter', delimiter);

fclose(inputfile); 
h= length(values{1});
for j=1:h
    user = values{1}(j);
    item = values{2}(j);
    rating = values{3}(j);
    Y2((item-1)*943+user)=rating;
    
end




IDX = find(Y);
IDX2= find(Y2);
S.type = '()';
S.subs{:} = IDX;

A = @(X) subsref(X,S);
Ah = @(X) subsasgn(zeros(m,n),S,X);
AhA = @(X) X.*Y;



%% Non-Negative Matrix Factorization
[P,Q] = nnmf(Y, 300);
XRecon = P*Q;
sum=0;
for j=IDX2
    add=abs(Y2(j)-XRecon(j));
end

x=length(IDX2);
ss=0;
for j=1:x
    ss=ss+add(j);
end
ss
res=ss/x


