function [Xf, P, SBIC]=arp_roll_forecast_h(Y,years,names,f,w)

[t,n]=size(Y);

pmax=12;
f2=strcat('for',num2str(f));
w2=strcat('w',num2str(w));

Xf=NaN(t,n);
Xf(1:w,:)=Y(1:w,:);

P=NaN(n,t-w-f+1);
SBIC=NaN(n,t-w-f+1);

for i=1:n
y=Y(:,i);
[~,~,p]=estima_select_arp_sbic(y,pmax);

if p>0 && w<=f
    Xf(w+f:t,i)=NaN(t-w-f+1,1);
    SBIC=NaN(n,t-w-f+1);
    fprintf('The window is too small\n');
else

for j=1:t-w+1-f   
X=Y(j:w+j-1,:);
x=X(:,i);
t1=length(x);
[~,~,p]=estima_select_arp_sbic(x,pmax);
P(i,j)=p;
[beta, y_hat]=estima_arp_direct_h(x,p,f);
res=x(p+f+1:t1,1)-y_hat;
[aic,sbic,hqc]=information_criteria(res,x,p+1);
SBIC(i,j)=sbic;

z=[];
if p>0
for l=0:p-1
z=[z,x(w-l)];
end
z=[1,z];
xf=z*beta;
elseif p==0
xf=beta;
end

Xf(w+j+f-1,i)=xf;
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
print(f,strcat('Figure_arp_roll_',f2,'_',w2), '-dpdf')
print(f,strcat('Figure_arp_roll_',f2,'_',w2), '-deps')

save(strcat(f2,'_',w2,'_arp_roll.dat'),'Xf','-ASCII');

