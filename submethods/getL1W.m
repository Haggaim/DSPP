function W = getL1W(D1,D2,params)
%compute the matrix W representing the quadratic form x'*W*x
% defined by \sum_qrst X_qr*X_st |D2(q,s)-D1(r,t)|

n = params.n;
k = params.k;
D2_4D = zeros(k,n,k,n);
D1_4D = zeros(k,n,k,n);

for q = 1:k
    for r = 1:n
        for s = 1:k
            for t = 1:n
                D2_4D(q,r,s,t)=D2(q,s);
                D1_4D(q,r,s,t)=D1(r,t);                
            end
        end
    end
end
D2_2D = reshape(D2_4D,[n*k n*k]);
D1_2D = reshape(D1_4D,[n*k n*k]);

differenceMAt_2D = abs(D2_2D-D1_2D);

W = differenceMAt_2D;

end