
function [X_proj,X_opt,final_obj,opt_obj] = solveDSpp(D1,D2,problemType)
%-------------------------------------------------------------------------
%input:  Two symmetric matrices: n by n matrix D1 and k by k matrix D2
%        where n>=k.
%        a choice of problemType between two alternatives:
%        'L1' minimizes \sum_{qrst} X_qr X_st |D1(q,s)-D2(r,t)| over 
%        "permutation matrices" X.
%        'GW' (default) minimizes -tr(X'D1XD2) over permutations X. Solving for this
%        objective is faster than for 'L1' since evaluating the objective needs O(n^3) actions instead of O(n^4) ;
%
%output: Solution of the convex relaxation X_convex and its objective
%        obj_convex.
%        Final integer solution X_final and its objective obj_final.
%---------------------------------------------------------------------------------
% Code written by Haggai Maron and Nadav Dym. Please cite
% "DS++: A flexible, scalable and provably tight relaxation for matching problems"
% questions regarding the code can be sent to nadavdym@gmail.com
%---------------------------------------------------------------------------------

n= size(D1,1);
k= size(D2,1);
params.n=n;
params.k=k;
params.problemType=problemType;
params.D1=D1;
params.D2=D2;
if strcmp(params.problemType,'L1')
   params.W=getL1W(D1,D2,params);
end
translationVec = linspace(0,1.1,10);
plotMatrices= getoptions(params,'plotMatrices',1);
params.injective=getoptions(params,'injective',~(n==k));
params.linTerm = getoptions(params,'linTerm',colstack(zeros(k,n)));
% find max translation

[lambMin,lambMax] = findLambdaMinAndMaxFast(params);

maximalTranslation=-1*lambMin;
minimalTranslation=-lambMax;
X_opt_arr = cell(1,length(translationVec));
for ii=1:length(translationVec)
    t=translationVec(ii);
    params.translation=t*minimalTranslation+(1-t)*maximalTranslation;
    [X_opt,local_opt_obj] = solveDsQuadprog(params);
    if ii==1  %this is the objective from the convex relaxation
        opt_obj=local_opt_obj;
    end
    X_opt_arr{ii} = X_opt;
    params.Xinit=X_opt;
end
% plot objectives
objs= zeros(1,numel(translationVec));
for ii = 1:numel(translationVec)
    X = X_opt_arr{ii};
    %         objs(ii) =  X(:)'*W*X(:);
    objs(ii) =  X(:)'*calcWProdFast(X(:),0,params);
end
%plot output matrices
if plotMatrices
    if params.injective %plot without extra  row as well
        figure;
        hold on;
        for ii=1:length(translationVec)
            subplot(2,ceil(length(translationVec)/2),ii);
            imagesc(X_opt_arr{ii}(2:k+1,:));
            colorbar;
            title(sprintf('Hysterestic iter no %d',ii ));
        end
    end
else
    figure;
    hold on;
    for ii=1:length(translationVec)
        subplot(2,ceil(length(translationVec)/2),ii);
        imagesc(X_opt_arr{ii});
        colorbar;
        title(sprintf('iter no %d',ii ));
    end
end
%project
final_obj = X_opt(:)'*calcWProdFast(X_opt(:),0,params);

if params.injective
    X_opt=X_opt(2:k+1,:);
end
[~,maxInd]=max(X_opt');
X_proj = zeros(k,n);
X_proj(sub2ind([k n],1:k,maxInd))=1;
%In paper we compute closest permutation using linear assignment, here we
%use simpler rounding so as not to require a solver for linear assignment


end