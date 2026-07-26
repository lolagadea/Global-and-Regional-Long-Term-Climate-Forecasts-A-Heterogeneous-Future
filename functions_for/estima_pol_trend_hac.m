%October 2012
%Looking for a trend, revised January 2014
function [beta,t_ratio,y_hat,varres,varbeta]=estima_pol_trend_hac(y,k)
t=length(y); 
trend=1:1:t; trend=trend';
if k==0
    X=ones(t,1);
else
trend_k=[];
for i=1:k
    trend_k=[trend_k,trend.^i];
end
    X=[ones(t,1),trend_k];   
end
    [beta,t_ratio,~,~,~,varres,varbeta]=ols_hac_forecast(y,X);
    y_hat=X*beta;


 