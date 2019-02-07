function out =  calcWProdFast(v,Wtranslation,params)
%fast calculation of out=(W+Wtranslation*I)*v
m=params.k;
n=params.n;

%% calculate W*x
switch params.problemType
    case 'L1'
        if params.injective
            Xaug=reshape(v,[m+1 n]);
            outMat=zeros(size(Xaug));
            X=Xaug(2:(m+1),:);
            x=X(:);
            subVec=params.W*x;
            subMat=reshape(subVec, [m n]);
            outMat(2:(m+1),:)=subMat;
            out2=outMat(:);
        else
            out2 = params.W*v;
        end
    case 'GW'
        if params.injective
            Xaug=reshape(v,[m+1 n]);
            outMat=zeros(size(Xaug));
            X=Xaug(2:(m+1),:);
            subMat=-2*params.D2*X*params.D1+ones(m,1)*(  ones(1,m)*X*params.D1.^2  );
            outMat(2:(m+1),:)=subMat; %+Wtranslation*X;
            out2=outMat(:);
        else
            out2 = -(colstack(params.D2*reshape(v,[params.n params.n])*params.D1));
        end
end


%% translate final solution

out = out2+Wtranslation*v;

end


