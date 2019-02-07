function [X,v,w,iter_num,dist_v,dist_w] = sinkhorn(K,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this code was adapted from the code released by the authors of
% "Entropic Metric Alignment for Correspondence Problems"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxSinkhornIter = getoptions(params,'maxSinkhornIter',1000 );
verbose = getoptions(params,'verbose',false);
n = size(K,1);
sinkhornTol = getoptions(params,'sinkhornTol',10^-6);
v = getoptions(params,'v',ones(n,1));
w = getoptions(params,'w',ones(n,1));
dist_v=inf(1,maxSinkhornIter);
dist_w=inf(1,maxSinkhornIter);
Kw = K*w;
for i=1:maxSinkhornIter % check this...
    if sum(Kw==0)>0
        error('division by zero')
    end
    v = 1./Kw;
    Kv = K'*v;
    if sum(Kv==0)>0
        error('division by zero')
    end
    w = 1./Kv;
    Kw = K*w;
    %compute distance from correct marginals
    dist_v(i)=log(max(w.*Kv)/min(w.*Kv));
    dist_w(i)=log(max(v.*Kw)/min(v.*Kw));
    % Stop when margins are close to 1
    %if norm((v.*Kw-1),'inf')<sinkhornTol && norm((w.*Kv-1),'inf')<sinkhornTol
    if i>1 
        ratio(i-1)=dist_w(i)/dist_w(i-1);
        if ratio(i-1)>0.9
            break;
        end
    end
    if (dist_v(i)<sinkhornTol && dist_w(i)<sinkhornTol) 
        break;
    end
end

Mv=repmat(v,[1,n]);
Mw=repmat(w',[n,1]);
X=Mv.*K.*Mw;
%X = diag(v)*K*diag(w);
if verbose
    fprintf('Sinkhorn projection ended after %d iterations...\n',i);
end
iter_num=i;

% if iter_num>1
%     figure; plot(ratio);
% else
%     figure; plot(zeros(1,2));
% end
% title('decay ratio of d');


end
