%July 2019
%Estimate a pol_trend_arp (k,p) model
%Select k and p globally by SBIC
function [p,k,beta,t_ratio,y_hat,res]=estima_select_pol_trend_arp_hac_sbic(y,maxpol,pmax)
t=length(y); 
SBIC=NaN(maxpol,pmax);
for i=1:maxpol
    for j=1:pmax
        Xarp=compute_arp(y,j);
        Xpol=compute_trend(i,t); 
        Xpol=Xpol(j+1:end,:);
        z=y(j+1:t,1);
        X=[ones(t-j,1),Xpol,Xarp];
        [~,~, se,res, ~, ~, ~, ~]=ols_hac_forecast(z,X);
        SBIC(i,j)=sbic(res,i+j+1);
    end
end

[~,I] = min(SBIC(:));
[I_row, I_col] = ind2sub(size(SBIC),I);
p=I_col;
k=I_row;

z=y(p+1:t,1);
Xpol=compute_trend(k,t-p);
Xarp=compute_arp(y,p);
X=[ones(t-p,1),Xpol,Xarp];
[beta,t_ratio, se,res, ~, ~, ~, y_hat]=ols_hac_forecast(z,X);


function TREND=compute_trend(k,t)
trend=1:1:t; trend=trend';
TREND=[];
    for i=1:k
    TREND=[TREND,trend.^i];
    end
    
function x=compute_arp(y,p)
%This function produces a matrix with p lags of y
t=length(y);
x=zeros(t-p,p);
%let's now fill the columns on x...
for j=1:p
    x(:,j)=y(p-j+1:t-j);
end

%The classical SBIC
function[y]=sbic(res,k)
  T=length(res);
  sigma2_res = res'*res/(T);
  y=log(sigma2_res)+log(T)*k/T;

          
