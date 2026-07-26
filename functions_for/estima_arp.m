%July 2019
%Estimate a ar(p)
%Include the case p=0
function  [beta,t_ratio,p,y_hat,res]=estima_arp(y,p)
[beta, t_ratio, res,y_hat]=arols(y,p);

        
function[beta, t_ratio, res, y_hat]=arols(y,p)
t=length(y);
x=zeros(t-p,p);
%let's now fill the columns in x...
for j=1:p
    x(:,j)=y(p-j+1:t-j);
end
x=[ones(t-p,1) x];

[beta,t_ratio, se,res, ~, ~, ~, y_hat]=ols_hac_forecast(y(p+1:t),x);

