function out = myEntropy(x)
%compute entropy of the matrix x
if size(x,1)<size(x,2) x=x'; end
    
assert(sum(x<0)==0)

zeroIdx = (x ==0);
out = -x(~zeroIdx)'*log(x(~zeroIdx));



end