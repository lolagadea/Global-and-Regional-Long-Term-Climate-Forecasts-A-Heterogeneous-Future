function [Xf, SBIC]=arfima_roll_forecast_h(Y,years,names,f,w)

[t,n]=size(Y);


f2=strcat('for',num2str(f));
w2=strcat('w',num2str(w));

Xf=NaN(t,n);
Xf(1:w,:)=Y(1:w,:);

SBIC=NaN(n,t-w-f+1);


for i=1:n
y=Y(:,i);
%ar_part=select_arfima(y);
ar_part=0; %estimation problems in some points
for j=1:t-w+1-f    
X=Y(j:w+j-1,:);
x=X(:,i);
[d,ar,sigma2, res]=estima_arfima(x,ar_part);
[aic,sbic,hqc]=information_criteria(res,x,2+length(ar));
SBIC(i,j)=sbic;
z=[];
xf= arfima_forecast(y,f,d,ar,[],mean(x),sigma2);
Xf(w+j+f-1,i)=xf(end);
end
end

Xf=Xf(1:t,:);

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
print(f,strcat('Figure_arfima_roll_',f2,'_',w2), '-dpdf')
print(f,strcat('Figure_arfima_roll_',f2,'_',w2), '-deps')

save(strcat(f2,'_',w2,'_arfima_roll.dat'),'Xf','-ASCII');