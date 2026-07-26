
function [Xf, SBIC]=ima_roll_forecast_h(Y,years,names,f,w)

[t,n]=size(Y);

f2=strcat('for',num2str(f));
w2=strcat('w',num2str(w));

Xf=NaN(t,n);
Xf(1:w,:)=Y(1:w,:);

SBIC=NaN(n,t-w-f+1);

for i=1:n
for j=1:t-w+1-f 
X=Y(j:w+j-1,:);
x=X(:,i);
dx=diff(x,1);
mu=mean(dx);
xm=x-mean(x);
dxm=diff(xm,1);
ima=armax(dxm,[0 1]);
theta=-ima.c(2);
sigma2=ima.Report.Fit.MSE;
l=length(x);
SBIC(i,j)=log(sigma2)+log(l)*1/l;
eps=x(w); 
for k=1:w-1
    eps=eps+theta^(k-1)*(theta-1)*x(w-k);
end
if theta<0.97 
eps=eps-mu/(1-theta);
xf=mu*f+x(w)-theta*eps;
elseif theta>0.99
xf=x(1)+(w+f)*mu;
else
eps=eps-mu/(1-theta);
xf=x(1)+mu*f+x(w)-theta*eps;
end
Xf(w+j+f-1,i)=xf;
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
print(f,strcat('Figure_ima_roll_',f2,'_',w2), '-dpdf')
print(f,strcat('Figure_ima_roll_',f2,'_',w2), '-deps')

save(strcat(f2,'_',w2,'_ima_roll.dat'),'Xf','-ASCII');
