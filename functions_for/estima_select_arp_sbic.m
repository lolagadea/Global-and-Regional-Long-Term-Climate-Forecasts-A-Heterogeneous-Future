%July 2019
%Estimate a ar(p) model; p selected with SBIC
function  [beta,t_ratio,p,y_hat,res]=estima_select_arp_sbic(y,pmax)
t=length(y); 
p=select_p_sbic(y,pmax);
[beta, t_ratio, res, y_hat]=arols(y,p);
         
  
function[p]=select_p_sbic(y,pmax)
for i=1:pmax
   [~, ~,res]=arols(y,i);
   SBIC(i)=sbic(res,i+1);
end
[~,p]=min(SBIC);


function[beta, t_ratio, res, y_hat]=arols(y,p)
t=length(y);
x=zeros(t-p,p);
%let's now fill the columns in x...
for j=1:p
    x(:,j)=y(p-j+1:t-j);
end
x=[ones(t-p,1) x];
[beta,t_ratio,se, res, ~, ~, ~, y_hat]=ols_hac_forecast(y(p+1:t),x);


  %The classical SBIC
function[y]=sbic(res,k);
  T=length(res);
  sigma2_res = res'*res/(T);
  y=log(sigma2_res)+log(T)*k/T;
  

