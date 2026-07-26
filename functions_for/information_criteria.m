function[aic,sbic,hqc]=information_criteria(res,xt,k)
if isempty(res)==0
  T=length(res);
  sigma2_res = res'*res/(T);
else
    T=length(xt);
    sigma2_res=var(diff(xt,1));
end
  
  aic=log(sigma2_res)+2*(k+1)/(T);
  sbic=log(sigma2_res)+log(T)*k/T;
  hqc=log(sigma2_res)+2*log(log(T))*k/T;

  