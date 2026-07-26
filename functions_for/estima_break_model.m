%estima_break_model
function  [b,t_nw, se_nw, res, r2, varres, varBhat, y_hat]=estima_break_model(y,TB,model)
T=length(y);
constant=ones(T,1);
trend=(1:T)';
nb=length(TB);
for i=1:nb
    DU(:,i)=(trend>TB(i)); DT(:,i)=(trend>TB(i)).*(trend-TB(i)); 
end
    if model==1
  		reg=[constant,trend,DU];
    elseif model==2
    	reg=[constant,trend,DT];
    elseif model==3
    	reg=[constant,DU,trend,DT];
    end

 [b,t_nw, se_nw, res, r2, varres, varBhat, y_hat]=ols_hac_forecast(y,reg);
 