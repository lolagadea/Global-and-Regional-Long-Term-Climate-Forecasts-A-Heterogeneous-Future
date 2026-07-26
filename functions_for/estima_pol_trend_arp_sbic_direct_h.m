%July 2019
%Looking for a trend
%Estimate a pol_trend_arp (k,p) model
function [beta,t_ratio,y_hat,r2, varres, varBhat]=estima_pol_trend_arp_sbic_direct_h(y,k,p,h)
t=length(y); 
trend=1:1:t; trend=trend';
TREND=[];

for i=0:k
    TREND=[TREND,trend.^i];
end
    X=TREND;
    if p>=1
    [beta,t_ratio,~,y_hat,r2,varres,varBhat]=arols_trend_direct_h(y,X,p,h);
    else 
    [beta,t_ratio,~,y_hat,r2, varres, varBhat]=estima_pol_trend(y,k);
    end
    

function [beta,t_ratio,res,y_hat,r2,varres,varBhat]=arols_trend_direct_h(y,xt,p,h)


X=embed(y,p); %build a matrix of lags

%beta=[ones(size(X,1)-h+1,1),X(1:end-h+1,2:end)]\X(h:end,1);
[beta,t_ratio, se,res, r2, varres, varBhat, y_hat]=ols_hac_forecast(X(h+1:end,1),[xt(p+h+1:end,:),X(1:end-h,2:end)]);

function [beta,t_ratio,res,y_hat,r2, varres, varBhat]=estima_pol_trend(y,k)
t=length(y);
trend=1:1:t; trend=trend';
TREND=[];
for i=0:k
    TREND=[TREND,trend.^i];
    X=[TREND];
end
[beta,t_ratio, se,res, r2, varres, varBhat, y_hat]=ols_hac_forecast(y,X);






          
