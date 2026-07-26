function [CI, ci_models]=compute_CI_comb_select_models(Y,RDOS,quantile,w,f,selected_models,Weig,gob_reg)
q=quantile;
Z=Y(:,q);
t=length(Z);
F=[1,10,25,50];
if gob_reg==0
W=[50,75,100]; 
elseif gob_reg==1
    W=25;
end
select_w=find(W==w);
select_f=find(F==f);
m=14;
m2=nansum(selected_models);

L=NaN(t-w-f+1,m2);
l=t-w-f+1;
f2=strcat('for',num2str(f));
w2=strcat('w',num2str(w));
Y=squeeze(RDOS.Xf(:,q,1:m,select_w,select_f));
Y=Y(:,selected_models==1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Combinated models with SBIC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y_comb_sbic=Y*Weig';
    l=t-w-f+1;
    y0=Z(w+f:t);
    for k=1:m2
    L(:,k)=(y0-Y(w+f:t,k));
    end
    FE=Y_comb_sbic(w+f:t,1)+L*Weig';
    var_fe=compute_cov_FE(L,Weig);
    CI=[norminv(0.025)*sqrt(var_fe),norminv(0.975)*sqrt(var_fe)];
    var_models=diag(cov(L));
    ci_models=[norminv(0.025)*sqrt(var_models),norminv(0.975)*sqrt(var_models)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function var_fe=compute_cov_FE(L,Weig)
[t,m]=size(L);
S=cov(L);
var_fe=Weig*diag(S);
if m>1    
pair=combinator(m,2,'c');
n=size(pair,1);
for i=1:n
    var_fe=var_fe+2*Weig(pair(i,1))*Weig(pair(i,2))*S(pair(i,1),pair(i,2));
end
end

