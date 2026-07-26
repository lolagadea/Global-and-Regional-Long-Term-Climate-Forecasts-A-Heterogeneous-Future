%October 2012
%Looking for a trend, revised January 2014
function [k,beta,t_ratio,y_hat,varres,varbeta]=estima_select_pol_trend_log_hac_sbic(Y,maxpol);
t=length(Y); 
trend=1:1:t; trend=trend';
ltrend=log(trend);
SBIC=zeros(maxpol,1);
    y=Y;
    t=length(y);
    TREND=[];   
    for i=1:maxpol
    TREND=[TREND,ltrend.^i];
    X=[ones(t,1),TREND];
    [beta,t_ratio,se,res]=ols_hac_forecast(y,X);
    SBIC(i)=sbic(res,i+1);
    end

[value,k]=min(SBIC); 
  [beta,t_ratio,se,res]=ols_hac_forecast(y,ones(t,1)); %model without trend
  sbic_c=sbic(res,1);
  if abs(sbic_c)<value
      k=0;
      X=[ones(t,1)];
    [beta,t_ratio,se,res,r2,varres,varbeta]=ols_hac_forecast(y,X);
    y_hat=X*beta;
    %plot(ltrend,y,'b',ltrend,y_hat,'r');
  else
ltrend_k=[];
for i=1:k
    ltrend_k=[ltrend_k,ltrend.^i];
end
    X=[ones(t,1),ltrend_k];
    [beta,t_ratio,se,res,r2,varres,varbeta]=ols_hac_forecast(y,X);
    y_hat=X*beta;
    %plot(ltrend,y,'b',ltrend,y_hat,'r');
  end
  
  %The classical SBIC
function[y]=sbic(res,k)
  T=length(res);
  sigma2_res = res'*res/(T);
  y=log(sigma2_res)+log(T)*k/T;
  
