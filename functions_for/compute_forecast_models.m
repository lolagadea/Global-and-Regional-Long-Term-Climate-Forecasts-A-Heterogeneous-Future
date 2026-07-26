%Compupe forecast with different models

function [FORECAST,CI, RES]=compute_forecast_models(y,f)
warning('off')

t=length(y);
m=14;
pmax=12;
maxpol=12;

FORECAST=NaN(m,1);
CI=NaN(m,2);
RES=NaN(t,m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Forecast with models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MEAN Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%xf_mean=mean(y);
[b,t_nw, se_nw,res, r2, varres, varBhat, y_hat]=ols_hac_forecast(y,ones(t,1));

xf_mean=b;
FORECAST(1)=xf_mean;

ci_mean=xf_mean+norminv([0.05 0.95])*sqrt(varres); %conditional
CI(1,:)=ci_mean;
RES(:,1)=res;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LINEAR-TREND Model%%%%%%%%%%%%%%%%%%%%%%%%%%
k=1;
[beta,t_ratio,y_hat,varres,varbeta]=estima_pol_trend_hac(y,k);

res=y-y_hat;
trend_f=[1,t+f];

xf_linear_trend=beta'*trend_f';


FORECAST(2)=xf_linear_trend;

ci_linear_trend=xf_linear_trend+norminv([0.05 0.95])*sqrt(varres); %conditional
CI(2,:)=ci_linear_trend;
RES(:,2)=res;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%POL-TREND Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=estima_select_pol_trend_hac_sbic(y,maxpol);
[beta,t_ratio,y_hat,varres,varbeta]=estima_pol_trend_hac(y,k);
res=y-y_hat;
slope=compute_slope(beta,t);
if k>0
    trend_f=[];
for i=0:k
    trend_f=[trend_f,(t+f)^i];
end

xf_pol_trend=beta'*trend_f';
ci_pol_trend=xf_pol_trend+norminv([0.05 0.95])*sqrt(varres); %conditional
beta2=[beta(1),slope];
trend_f2=[1,(t+f)];
xf2_pol_trend=beta2*trend_f2'; %comppute_pol_trend_average_slope model
ci2_pol_trend=xf2_pol_trend+norminv([0.05 0.95])*sqrt(varres); %conditional
elseif k==0
    xf_pol_trend=beta;
    xf2_pol_trend=beta;
    ci_pol_trend=xf_pol_trend+norminv([0.05 0.95])*sqrt(varres);
    ci2_pol_trend=ci_pol_trend;
end

FORECAST(3)=xf_pol_trend;
FORECAST(4)=xf2_pol_trend;
CI(3,:)=ci_pol_trend;
CI(4,:)=ci2_pol_trend;
RES(:,3)=res;
RES(:,4)=res;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%POL-TREND-LOG Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=estima_select_pol_trend_log_hac_sbic(y,maxpol);
[beta,t_ratio,y_hat,varres,varbeta]=estima_pol_trend_log_hac(y,k);
res=y-y_hat;
if k>0
    trend_f=[];
for i=0:k
    trend_f=[trend_f,log((t+f))^i];
end

xf_pol_trend_log=beta'*trend_f';
elseif k==0
    xf_pol_trend_log=beta;
end

FORECAST(5)=xf_pol_trend_log;

ci_pol_trend_log=xf_pol_trend_log+norminv([0.05 0.95])*sqrt(varres); %conditional
CI(5,:)=ci_pol_trend_log;
RES(:,5)=res;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%STRUCT-BREAKS Model%%%%%%%%%%%%%%%%%%%%%%%%%%
model=3;  
criteria=2; %BIC if criteria=2, AIC if criteria=1
eps=0.15;	%chose a value from {0.01, 0.05, 0.10, 0.15, 0.25}
kmax=fix(12*(t/100)^(1/4));
[wald,cv,tb]=qfgls_simus(y,kmax,model, criteria,eps);
if wald>cv(2)
[beta_break,t_nw, se_nw, res, r2, varres, varBhat, y_hat]=estima_break_model(y,tb,model);
xf_struct_breaks=beta_break(1)+beta_break(2)+beta_break(4)*(t+f);
ci_struct_breaks=xf_struct_breaks+norminv([0.05 0.95])*sqrt(varres); %conditional
else 
k=1;
[beta,t_ratio,y_hat,varres,varbeta]=estima_pol_trend_hac(y,k);
res=y-y_hat;
trend_f=[1,t+f];
xf_struct_breaks=beta'*trend_f';
ci_struct_breaks=xf_struct_breaks+norminv([0.05 0.95])*sqrt(varres); %conditional
end
FORECAST(6)=xf_struct_breaks;
CI(6,:)=ci_struct_breaks;
RES(:,6)=res;


%%%%%%%%%%%%%%%%%%%%%%%%%Pol-trend+AR(p) Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[p,k]=estima_select_pol_trend_arp_hac_sbic(y,maxpol,pmax);
[beta,t_ratio,y_hat,r2, varres, varBhat]=estima_pol_trend_arp_sbic_direct_h(y,k,p,f);
res=y(p+f+1:end)-y_hat;

if p+f>=length(y)
    xf_pol_trend_arp=NaN;
    ci_pol_trend_arp=[NaN NaN];
else

if p>0 && k>0
z=[];    
for w=0:p-1
z=[z,y(t-w)];
end

trend_f=[];
for i=0:k
    trend_f=[trend_f,(t+f)^i];
end
Z=[trend_f';z'];
xf_pol_trend_arp=beta'*Z;

elseif p==0 && k==0
xf_pol_trend_arp=beta;

elseif p==0 && k>0
    trend_f=[];
for i=0:k
    trend_f=[trend_f,(t+f)^i];
end
xf_pol_trend_arp=beta'*trend_f';  

elseif p>0 && k==0    
z=[];    
for w=0:p-1
z=[z,y(t-w)];
end 
z=[1,z];
xf_pol_trend_arp=z*beta;
end

FORECAST(7)=xf_pol_trend_arp;
ci_pol_trend_arp=xf_pol_trend_arp+norminv([0.05 0.95])*sqrt(varres); %conditional
CI(7,:)=ci_pol_trend_arp;
end
RES(p+f+1:t,7)=res;

%%%%%%%%%%%%%%%%%%%%%%%%%Pol-trend+AR(p)-average-slope Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[p,k]=estima_select_pol_trend_arp_hac_sbic(y,maxpol,pmax);
[beta,t_ratio,y_hat,r2, varres, varBhat]=estima_pol_trend_arp_sbic_direct_h(y,k,p,f);
res=y(p+f+1:t)-y_hat;

if p+f>=length(y)
    xf_pol_trend_arp_av_sl=NaN;
    ci_pol_trend_arp_av_sl=[NaN NaN];
else

if p>0 && k>0
z=[];    
for w=0:p-1
z=[z,y(t-w)];
end
z=z';

slope=compute_slope(beta(1:k+1),t);
beta2=[beta(1);slope;beta(k+2:end)];
trend_f=[1,(t+f)];
Z=[trend_f';z];
xf_pol_trend_arp_av_sl=beta2'*Z; 

elseif p==0 && k==0
xf_pol_trend_arp_av_sl=beta;

elseif p==0 && k>0
slope=compute_slope(beta(1:k+1),t);
trend_f=[1,(t+f)];
xf_pol_trend_arp_av_sl=[beta(1),slope]*trend_f';

elseif p>0 && k==0    
z=1;    
for v=0:p-1
z=[z;y(t-v)];
end
    xf_pol_trend_arp_av_sl=beta'*z;
end

FORECAST(8)=xf_pol_trend_arp_av_sl;
ci_pol_trend_arp_av_sl=xf_pol_trend_arp_av_sl+norminv([0.05 0.95])*sqrt(varres); %conditional
CI(8,:)=ci_pol_trend_arp_av_sl;
end
RES(p+f+1:t,8)=res;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%AR(p) Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,~,p]=estima_select_arp_sbic(y,pmax);
[beta, y_hat, varres]=estima_arp_direct_h(y,p,f);
res=y(p+f+1:t)-y_hat;

if p+f>=length(y)
    xf_arp=NaN;
    ci_arp=[NaN NaN];
else

z=[];
if p>0
z=[];    
for w=0:p-1
z=[z,y(t-w)];
end 
z=[1,z];
xf_arp=z*beta;
elseif p==0
xf_arp=beta;
end


FORECAST(9)=xf_arp;
ci_arp=xf_arp+norminv([0.05 0.95])*sqrt(varres); %conditional
CI(9,:)=ci_arp;
end
RES(p+f+1:t,9)=res;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RW model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xf_rw=y(end);
res=my_diff(y,f);
varres=res'*res/(length(y)-1);

FORECAST(10)=xf_rw;
ci_rw=xf_rw+norminv([0.05 0.95])*sqrt(varres); %conditional
CI(10,:)=ci_rw;
RES(f+1:t,10)=res;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RW with dritf model%%%%%%%%%%%%%%%%%%%%%%%%
dy=diff(y,1);
reg=ones(t-1,1);
[alpha, ~, ~,~, ~, varres]=ols_hac_forecast(dy,reg);
res=my_diff(y,f)-alpha*f;
xf_rwd=y(t)+alpha*f;
varres=res'*res/(length(y)-2);

FORECAST(11)=xf_rwd;

ci_rwd=xf_rwd+norminv([0.05 0.95])*sqrt(varres); %conditional
CI(11,:)=ci_rwd;
RES(f+1:t,11)=res;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%IMA model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [whittle] = arfima_estimate(y,'FWHI',[0 1]);
% res=whittle.errors;

dy=diff(y,1);
mu=mean(dy);
ym=y-mean(y);
dym=diff(ym,1);
ima=armax(dym,[0 1]);
theta=-ima.c(2);
eps=y(t); 
for k=1:t-1
    eps=eps+theta^(k-1)*(theta-1)*y(t-k);
end
if theta<0.97 
eps=eps-mu/(1-theta);
xf_ima=mu*f+y(t)-theta*eps;
res=my_diff(y,f)-mu*f+theta*eps;
elseif theta>0.99
xf_ima=y(1)+(t+f)*mu;
res=my_diff(y,f)-mu*(t+f);
else
eps=eps-mu/(1-theta);
xf_ima=y(1)+mu*f+y(t)-theta*eps;
res=my_diff(y,f)-mu*f+theta*eps;
end
%varres=ima.Report.Fit.MSE;
varres=res'*res/(length(y)-2);

FORECAST(12)=xf_ima;
ci_ima=xf_ima+norminv([0.05 0.95])*sqrt(varres); %conditional
CI(12,:)=ci_ima;
RES(f+1:t,12)=res;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Arfima Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[d,ar,sigma2,res,y_hat]=estima_arfima(y,0);
xf= arfima_forecast(y,f,d,[],[],mean(y),sigma2);
xf_arfima=xf(end);
varres=sigma2;

FORECAST(13)=xf_arfima;
ci_arfima=xf_arfima+norminv([0.05 0.95])*sqrt(varres); %conditional
CI(13,:)=ci_arfima;
RES(:,13)=res;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%AR20(p) Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p=20;
[beta, y_hat, varres]=estima_arp_direct_h(y,p,f);
res=y(p+f+1:t)-y_hat;

if p+f>=length(y)
    xf_arp20=NaN;
    ci_arp20=[NaN NaN];
else

z=[];    
for w=0:p-1
z=[z,y(t-w)];
end 
z=[1,z];
xf_arp20=z*beta;

FORECAST(14)=xf_arp20;
ci_arp20=xf_arp20+norminv([0.05 0.95])*sqrt(varres); %conditional
CI(14,:)=ci_arp20;
end
RES(p+f+1:t,14)=res;

