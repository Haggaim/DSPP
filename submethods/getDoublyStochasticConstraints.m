function [A,b] = getDoublyStochasticConstraints(n,k)
%represent the affine DS constraints Ax=b where x is the columnstack DS matrix X %
if nargin < 2
    k = n;
end
T = getTensorTranspose(n,k);
CRows = kron(ones(k,1)',eye(n,n)); % sums up rows
CCols = kron(ones(n,1)',eye(k,k)) * T; % sums up cols
A = [CRows;CCols];
b = ones(n + k,1);
end

function T = getTensorTranspose(n,m)
% Returns a matrix T such that vec(X')=T*vec(X) for any nxm matrix X.
%
% Example:
% n = 3;
% m = 4;
% X = zeros(n,m);
% X(:)=1:numel(X);
% T = getTensorTranspose(n,m);
% vec(X')
% T*vec(X)

[i,j] = meshgrid(1:n,1:m);
T = sparse(sub2ind([m n],j,i),sub2ind([n m],i,j),1);
end
