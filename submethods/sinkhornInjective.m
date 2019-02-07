function [X,v,w] = sinkhornInjective(K,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this code was adapted from the code released by the authors of
% "Entropic Metric Alignment for Correspondence Problems"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.null = [];
maxSinkhornIter = getoptions(params,'maxSinkhornIter',500 );
verbose = getoptions(params,'verbose',0);
sinkhornTol = getoptions(params,'sinkhornTol',10^-6);
n=params.n;
m=params.k;
colSum=ones(1,n);
rowSum=ones(m+1,1);
rowSum(1,1)=n-m;
%check wheter there is a zero row or column
removeCol=(colSum(1)==0);
if removeCol
    K=K(:,2:end);
    colSum=colSum(2:end);
end
removeRow=(rowSum(1)==0);
if removeRow
    K=K(2:end,:);
    rowSum=rowSum(2:end);
end

[a,b]=size(K);
v = getoptions(params,'v',ones(a,1));
w = getoptions(params,'w',ones(b,1));


Kw = K*w;
for i=1:maxSinkhornIter % check this...
    
    if sum(Kw==0)>0
       error('division by zero') 
    end
    v = (1./Kw).*(rowSum);
    Kv = K'*v;
    if sum(Kv==0)>0
       error('division by zero') 
    end
    w = 1./Kv.*(colSum)';
    Kw = K*w;
    % Stop when margins are close to 1
    if norm((v.*Kw-rowSum),1)<sinkhornTol && norm((w.*Kv-colSum'),1)<sinkhornTol
        break;
    end
end
X = diag(v)*K*diag(w);
if verbose
    fprintf('Sinkhorn projection ended after %d iterations...\n',i);
end
%bring back delted row and/or column if necessary
if removeCol
    X=[zeros(size(X,1),1) X];
end
if removeRow
    X=[zeros(1,size(X,2)); X];
end
end
