function [Xf,TB,YEARS_breaks, SBIC]=struct_break_roll_forecast_h(Y,years,names,f,w)

%Model=1: Y{t}=a0+a1*DU+b0*t+e{t} where DU=1(t>TB)
%Model=2: Y{t}=a0+b0*t+b1*DT+e{t} where DT=1(t>TB)*(t-TB)
%Model=3: Y{t}=a0+a1*DU+b0*t+b1*DT+e{t} where DU=1(t>TB) and DT=1(t>TB)*(t-TB)
model=3;  
criteria=2; %BIC if criteria=2, AIC if criteria=1
eps=0.15;	%chose a value from {0.01, 0.05, 0.10, 0.15, 0.25}
kmax=fix(12*(w/100)^(1/4));

[t,n]=size(Y);

f2=strcat('for',num2str(f));
w2=strcat('w',num2str(w));

Xf=NaN(t,n);
Xf(1:w,:)=Y(1:w,:);


K=NaN(n,t-w-f+1);
SBIC=NaN(n,t-w-f+1);
TB=NaN(n,t-w-f+1);
YEARS_breaks=NaN(n,t-w-f+1);

k=1; %linear model if no breaks

for i=1:n  
for j=1:t-w+1-f
X=Y(j:w+j-1,:);
years_w=years(j:w+j-1);
x=X(:,i);
[wald,cv,tb]=qfgls_simus(x,kmax,model, criteria,eps);
if wald>cv(2)
[beta_break,t_nw, se_nw, res, r2, varres, varBhat, y_hat]=estima_break_model(x,tb,model);
trend=(w+f);
xf=beta_break(1)+beta_break(2)+beta_break(4)*trend;
[aic,sbic,hqc]=information_criteria(res,x,4);
TB(i,j)=tb;
YEARS_breaks(i,j)=years_w(tb);
else
    trend=(w+f);
    trend=[1,trend];
[beta,t_ratio,y_hat,varres,varbeta]=estima_pol_trend_hac(x,k);
res=x-y_hat;
[aic,sbic,hqc]=information_criteria(res,x,2);
xf=trend*beta;
end  
SBIC(i,j)=sbic;
Xf(w+j-1+f,i)=xf;
end
end


[m1,m2]=compute_subplots(n);
T=1:1:t;
f=figure;
for i=1:n
subplot(m1,m2,i),
plot(T,Y(:,i),'b',T,Xf(:,i),'r');
title(names(i,:));      
if t>140
set(gca,'XTick',[1 20 40 60 80 100 120 t],'FontSize',8);
set(gca,'XTickLabel', years([1 20 40 60 80 100 120 t]),'FontSize',8);
elseif t<70
set(gca,'XTick',[1 10 20 30 40 50 t],'FontSize',8);
set(gca,'XTickLabel', years([1 10 20 30 40 50 t]),'FontSize',8);
else
    disp('Define other ticks and labels for the figures')
end
end
legend('original','forecast');
set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
print(f,strcat('Figure_struct_breaks_roll_',f2,'_',w2), '-dpdf')
print(f,strcat('Figure_pol_struct_breaks_roll_',f2,'_',w2), '-deps')


save(strcat(f2,'_',w2,'_pol_trend_roll.dat'),'Xf','-ASCII');
