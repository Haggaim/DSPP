function [lambMin,lambMax] = findLambdaMinAndMaxFast(params)
%compute maximal and minimal eigenvalues of W on the DS subspace
n = params.n;
m=params.k;


%find lambMin (for initial translation)

%preprocessing:steps not needed for each iteration
%create A, a full rank matrix whose columns define the linear constraints
preprocessed=[];

if params.injective
    [A_ds,~] = getDoublyStochasticConstraints(m+1,n);
    A=A_ds';
    A=A(:,1:n+m);  %because there are n+(m+1) constraints and one is redundant
    vDim=(m+1)*n;
else
    [A_ds,~] = getDoublyStochasticConstraints(n);
    A=A_ds';
    A=A(:,1:(2*n-1));
    vDim=n^2;
end
%create Ginv: inverse of Gram matrix
preprocessed.Ginv=inv(A'*A);
preprocessed.A=A;
%find lambaMin
opt.tol=10^-6;
funHandle=@(v) eigHelper(preprocessed,v,params);
[~,lambMaxValue]=eigs(funHandle,vDim,1,'lm',opt);
if lambMaxValue>0
    lambMax=lambMaxValue;
    opt.tol=10^(-6)*lambMax;
    Shift=1.1*lambMax;
    newFunHandle=@(v) Shift*projectAffinedDS(v,preprocessed)-eigHelper(preprocessed,v,params);
    [~,lambMinShifted]=eigs(newFunHandle,vDim,1,'lm',opt);
    lambMin=Shift-lambMinShifted;
else
    lambMin=lambMaxValue;
    opt.tol=10^(-4)*abs(lambMin);
    Shift=1.1*lambMin;
    newFunHandle=@(v)eigHelper(preprocessed,v,params)-Shift*projectAffinedDS(v,preprocessed);
    [~,lambMaxShifted]=eigs(newFunHandle,vDim,1,'lm',opt);
    lambMax=lambMaxShifted+Shift;
end
end
%------------------------------------------------------
%subfunctions
%------------------------------------------------------
function u=eigHelper(preprocessed,v,params)
%implement u=F'*W*F*v in the correct order
u=projectAffinedDS(v,preprocessed);
u = calcWProdFast(u,0,params);
u=projectAffinedDS(u,preprocessed);
end

function u=projectAffinedDS(v,preprocessed)
%project v onto the affine space defined by the DS matrices.
%n is the dimension of the affine constraints
Pv=preprocessed.A'*v;
Pv=preprocessed.Ginv*Pv;
Pv=preprocessed.A*Pv;
%project onto kernel
u=v-Pv;
end





