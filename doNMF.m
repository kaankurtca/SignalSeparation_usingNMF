function [B,W,err] = doNMF(V,niter,B,W)

F=size(V,1); T=size(V,2);
Ones=ones(F,T);
err=zeros(niter,1);
for i=1:niter
    W= W.* (B'*(V./(B*W+eps))) ./  (B'*Ones);
    B= B.* ((V./(B*W+eps))*W') ./  (Ones*W'); 
    
    error=((V-(B*W)).^2)./(F*T);    %MeanSquaredError is calculated in each iteration.
    err(i)=(sum(sum(error)));
end
sumB=sum(B);
B=B*diag(1./sumB);
W=diag(sumB)*W;

% Learning rule from course slides.

end