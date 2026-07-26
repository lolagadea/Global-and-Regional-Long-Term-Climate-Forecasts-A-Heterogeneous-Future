%April 2013
%Estimate an ar(p) model to forecast with the direct method
%Include the case p=0
function  [beta, y_hat, varres]=estima_arp_direct_h(y,p,h)
t=length(y);
if p==0
    %beta=ones(length(y),1)\y;
    [beta,t_nw, se,res, r2, varres, varBhat, y_hat]=ols_hac_forecast(y,ones(length(y),1));
else
X=embed(y,p); %build a matrix of lags
%beta=[ones(size(X,1)-h,1),X(1:end-h,2:end)]\X(h+1:end,1);
[beta,t_nw,se, res, r2, varres, varBhat, y_hat]=ols_hac_forecast(X(h+1:end,1),[ones(size(X,1)-h,1),X(1:end-h,2:end)]);
end
         
[T,k]=size([ones(size(X,1)-h,1),X(1:end-h,2:end)]);
if k>=T
    varres=res'*res;
end