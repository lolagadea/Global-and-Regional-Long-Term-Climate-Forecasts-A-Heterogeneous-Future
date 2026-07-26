function [Xf, K, P, SBIC]=pol_trend_arp_average_slope_roll_forecast_h(Y,years,names,f,w)
[t,n]=size(Y);

maxpol=12;
pmax=12;

f2=strcat('for',num2str(f));
w2=strcat('w',num2str(w));

Xf=NaN(t,n);
Xf(1:w,:)=Y(1:w,:);

K=NaN(n,t-w-f+1);
P=NaN(n,t-w-f+1);
SBIC=NaN(n,t-w-f+1);

for i=1:n  
    y=Y(:,i);
[p,k]=estima_select_pol_trend_arp_hac_sbic(y,maxpol,pmax);
if p>0 && w<=f
    Xf(w+f:t,i)=NaN(t-w-f+1,1);
    K=NaN(n,t-w-f+1);
    P=NaN(n,t-w-f+1);
    SBIC=NaN(n,t-w-f+1);
    fprintf('The window is too small\n');
else
for j=1:t-w+1-f
X=Y(j:w+j-1,:);
x=X(:,i);
[p,k]=estima_select_pol_trend_arp_hac_sbic(x,maxpol,pmax);
K(i,j)=k;
P(i,j)=p;
t1=length(x);
[beta,t_ratio,y_hat,r2, varres, varBhat]=estima_pol_trend_arp_sbic_direct_h(x,k,p,f);
res=x(p+f+1:t1,1)-y_hat;
[aic,sbic,hqc]=information_criteria(res,x,k+p+1);
SBIC(i,j)=sbic;
if k>0 && p>0
z=[];    
for v=0:p-1
z=[z;y(t-v)];
end
slope=compute_slope(beta(1:k+1),t);
beta2=[beta(1);slope;beta(k+2:end)];
trend_f=[1,(t+f)];
Z=[trend_f';z];
xf=beta2'*Z;  %compute_pol_trend_arp_average_slope model
%xf=[beta(1),slope]*trend_f'/(1-sum(beta(k+2:end)));
elseif k==0 && p~=0
    z=1;    
for v=0:p-1
z=[z;y(t-v)];
end
    xf=beta'*z;
elseif k~=0 && p==0
    slope=compute_slope(beta(1:k+1),t);
    trend_f=[1,(t+f)];
    xf=[beta(1),slope]*trend_f';
elseif k==0 && p==0
    xf=beta;
end
Xf(w+j-1+f,i)=xf;
end   
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
print(f,strcat('Figure_pol_trend_arp_av_sl_roll_',f2,'_',w2), '-dpdf')
print(f,strcat('Figure_arp_av_sl_roll_',f2,'_',w2), '-deps')

save(strcat(f2,'_',w2,'_arp_av_sl_roll.dat'),'Xf','-ASCII');