function [Xf, K, SBIC]=pol_trend_average_slope_roll_forecast_h(Y,years,names,f,w)

[t,n]=size(Y);

maxpol=12;

f2=strcat('for',num2str(f));
w2=strcat('w',num2str(w));

Xf=NaN(t,n);
Xf(1:w,:)=Y(1:w,:);

K=NaN(n,t-w-f+1);
SBIC=NaN(n,t-w-f+1);

for i=1:n  
    y=Y(:,i);
k=estima_select_pol_trend_hac_sbic(y,maxpol);

for j=1:t-w+1-f
X=Y(j:w+j-1,:);
x=X(:,i);
k=estima_select_pol_trend_hac_sbic(x,maxpol);
K(i,j)=k;
[beta,t_ratio,y_hat,varres,varbeta]=estima_pol_trend_hac(x,k);
res=x-y_hat;
[aic,sbic,hqc]=information_criteria(res,x,k+1);
SBIC(i,j)=sbic;
slope=compute_slope(beta,w);

if k>0
beta2=[beta(1),slope];
trend_f=[1,(w+f)];
xf=beta2*trend_f';

else xf=beta;
end
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
print(f,strcat('Figure_pol_trend_av_sl_roll_',f2,'_',w2), '-dpdf')
print(f,strcat('Figure_pol_trend_av_sl_roll_',f2,'_',w2), '-deps')


save(strcat(f2,'_',w2,'_pol_trend_av_sl_roll.dat'),'Xf','-ASCII');
