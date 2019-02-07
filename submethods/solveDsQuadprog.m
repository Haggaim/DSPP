function [X,obj_opt] = solveDsQuadprog(params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve objectives of the form x'Wx + linTerm'*x + <X,log(X)> ):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this code was adapted from the code released by the authors of
% "Entropic Metric Alignment for Correspondence Problems"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global alphaNormalizedVec

% params
n= getoptions(params,'n',6);
k=getoptions(params,'k',n);
maxIter = getoptions(params,'maxIter',140);
eta = getoptions(params,'eta',0.01);
alpha = getoptions(params,'alpha',0.01);
tol = getoptions(params,'tol',10^-8);
if params.injective
    Xinit=ones(k+1,n);
    Xinit=sinkhornInjective(Xinit,params);
    Xinit=getoptions(params,'Xinit',Xinit);
else
    Xinit=getoptions(params,'Xinit',ones(n)/n);             %optional: user defined initialization
end
% display options
%define linear term
linTerm = getoptions(params,'linTerm',colstack(zeros(k,n)));
if params.injective
    helperMat=zeros(k+1,n);
    helperMat(2:k+1,:)=reshape(linTerm,[k n]);
    linTerm=helperMat(:);
end
Wtranslation=params.translation;
% start KL iterations
running = 1;
ii = 0;
oldX = Xinit;
objs = inf(1,maxIter);
while running
    ii = ii+1;
    expArg = -(calcWProdFast(colstack(oldX),Wtranslation,params)+linTerm);
    originalExpArg=expArg;
    expArg=expArg-max(expArg(:)); 
        scaleFactor=max(abs(expArg(:)));
    alphaNormalized=alpha*scaleFactor;
    fX=reshape(exp(expArg/alphaNormalized),size(oldX));
    K = (fX.^eta) .* (oldX.^(1-eta));
    if params.injective
        [currX,v,w] = sinkhornInjective(K,params);
    else
        [currX,v,w] = sinkhorn(K,params);
    end
    params.v = v; params.w = w;
    % calc objectives
    objs(ii)=-originalExpArg'*currX(:)-alphaNormalized*myEntropy(currX(:));
    assert(~isnan(objs(ii)) && ~isinf(objs(ii)))
    % update
    oldX = currX;
    % stop criteria
    if ii>1 && ~isnan(objs(ii)) && ~isnan(objs(ii-1))
        change = abs(objs(ii)-objs(ii-1))/abs(objs(ii));
    else
        change = inf;
    end
    running = (ii<maxIter) && (change > tol);
end
X = currX;
obj_opt=X(:)'*calcWProdFast(X(:),Wtranslation,params)-Wtranslation*k;
end

