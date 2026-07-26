function beta=compute_betas_sum1_all(y,X)
global y X
m=size(X,2);
t=length(y);
%Optimization procedure
beta0=[ones(t,1),X]\y;
A=[];
b=[];
Aeq=[0,ones(1,m)]; beq=[1];
lB=[-Inf zeros(1,m)];
uB=[Inf ones(1,m)];
options=optimset('Algorithm','interior-point', 'Display','iter','LargeScale',...
'off','Hessian','bfgs','TolX',1e-5,'TolFun',1e-1000,'TolCon',1e-1000,'MaxFunEvals',...
50000,'MaxIter',5300,'FunValCheck','off','DiffMinChange',1e-8, 'DiffMaxChange', 0.1);
[beta,fval,exitflag,output,lambda,grad,hessian]=fmincon(@fun,beta0,A,b,Aeq,beq,lB,uB,[],options);
clear global
function fmin=fun(beta0)
global y X
t=length(y);
m=size(X,2);

if m==2
    res=beta0(1)*ones(t,1)+beta0(2)*X(:,1)+beta0(3)*X(:,2) - y;

elseif m==3
    res=beta0(1)*ones(t,1)+beta0(2)*X(:,1)+beta0(3)*X(:,2)+ beta0(4)*X(:,3) - y;

elseif m==4
    res=beta0(1)*ones(t,1)+beta0(2)*X(:,1)+beta0(3)*X(:,2)+ beta0(4)*X(:,3)+...
    beta0(5)*X(:,4) - y;

elseif m==5
    res=beta0(1)*ones(t,1)+beta0(2)*X(:,1)+beta0(3)*X(:,2)+ beta0(4)*X(:,3)+...
    beta0(5)*X(:,4)+ beta0(6)*X(:,5) - y;

elseif m==6
    res=beta0(1)*ones(t,1)+beta0(2)*X(:,1)+beta0(3)*X(:,2)+ beta0(4)*X(:,3)+...
    beta0(5)*X(:,4)+ beta0(6)*X(:,5)+ beta0(7)*X(:,6)- y;

elseif m==7
    res=beta0(1)*ones(t,1)+beta0(2)*X(:,1)+beta0(3)*X(:,2)+ beta0(4)*X(:,3)+...
    beta0(5)*X(:,4)+ beta0(6)*X(:,5)+ beta0(7)*X(:,6) +beta0(8)*X(:,7) - y;

elseif m==8
    res=beta0(1)*ones(t,1)+beta0(2)*X(:,1)+beta0(3)*X(:,2)+ beta0(4)*X(:,3)+...
    beta0(5)*X(:,4)+ beta0(6)*X(:,5)+ beta0(7)*X(:,6) +beta0(8)*X(:,7) +...
    beta0(9)*X(:,8) - y;

elseif m==9
    res=beta0(1)*ones(t,1)+beta0(2)*X(:,1)+beta0(3)*X(:,2)+ beta0(4)*X(:,3)+...
    beta0(5)*X(:,4)+ beta0(6)*X(:,5)+ beta0(7)*X(:,6) +beta0(8)*X(:,7) +...
    beta0(9)*X(:,8)+ beta0(10)*X(:,9) - y;

elseif m==10
    res=beta0(1)*ones(t,1)+beta0(2)*X(:,1)+beta0(3)*X(:,2)+ beta0(4)*X(:,3)+...
    beta0(5)*X(:,4)+ beta0(6)*X(:,5)+ beta0(7)*X(:,6) +beta0(8)*X(:,7) +...
    beta0(9)*X(:,8)+ beta0(10)*X(:,9)+ beta0(11)*X(:,10)- y;

elseif m==11
    res=beta0(1)*ones(t,1)+beta0(2)*X(:,1)+beta0(3)*X(:,2)+ beta0(4)*X(:,3)+...
    beta0(5)*X(:,4)+ beta0(6)*X(:,5)+ beta0(7)*X(:,6) +beta0(8)*X(:,7) +...
    beta0(9)*X(:,8)+ beta0(10)*X(:,9)+ beta0(11)*X(:,10)+ beta0(12)*X(:,11) - y;

elseif m==12
res=beta0(1)*ones(t,1)+beta0(2)*X(:,1)+beta0(3)*X(:,2)+ beta0(4)*X(:,3)+...
    beta0(5)*X(:,4)+ beta0(6)*X(:,5)+ beta0(7)*X(:,6) +beta0(8)*X(:,7) +...
    beta0(9)*X(:,8)+ beta0(10)*X(:,9)+ beta0(11)*X(:,10)+ beta0(12)*X(:,11)+beta0(13)*X(:,12) - y;

elseif m==13
res=beta0(1)*ones(t,1)+beta0(2)*X(:,1)+beta0(3)*X(:,2)+ beta0(4)*X(:,3)+...
    beta0(5)*X(:,4)+ beta0(6)*X(:,5)+ beta0(7)*X(:,6) +beta0(8)*X(:,7) +...
    beta0(9)*X(:,8)+ beta0(10)*X(:,9)+ beta0(11)*X(:,10)+ beta0(12)*X(:,11)+beta0(13)*X(:,12)+beta0(14)*X(:,13) - y;

end

fmin=res'*res;

