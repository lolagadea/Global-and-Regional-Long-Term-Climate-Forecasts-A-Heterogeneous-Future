%This is the main function of the procedure of Perron and Yabu: "Testing for Shifts in Trend with an Integrated or Stationary 
%Noise Component," Journal of Business and Economic Statistics 27 (2009), 369-396.

function [wald,cv,TB]=qfgls_simus(y,kmax,model, criteria,eps)
%***************************************
%define some matrix
if model==1
    VR=[0,0,1];
	v_t=[-4.30;-4.39;-4.39;-4.34;-4.32;-4.45;-4.42;-4.33;-4.27;-4.27];
    cv=[1.60,2.07,3.33;1.52,1.97,3.24;1.41,1.88,3.05;1.26,1.74,3.12;0.91,1.33,2.83];
elseif model==2
    VR=[0,0,1];
	v_t=[-4.27;-4.41;-4.51;-4.55;-4.56;-4.57;-4.51;-4.38;-4.26;-4.26];
    cv=[1.52,2.02,3.37;1.40,1.93,3.27;1.28,1.86,3.20;1.13,1.67,3.06;0.74,1.28,2.61];
elseif model==3
    VR=[0,1,0,0;0,0,0,1]; 
	v_t=[-4.38;-4.65;-4.78;-4.81;-4.90;-4.88;-4.75;-4.70;-4.41;-4.41];
    cv=[2.96,3.55,5.02;2.82,3.36,4.78;2.65,3.16,4.59;2.48,3.12,4.47;2.15,2.79,4.57];
end

T=length(y);
constant=ones(T,1);
trend=seqa(1,1,T);
vect1=zeros(fix((1-2*eps)*T)+2,1);
TBi=max([fix(eps*T);(kmax+2)]);

while TBi<=fix((1-eps)*T)
	lam1=TBi/T;
    DUi=(trend>TBi); DTi=(trend>TBi).*(trend-TBi); 
    if model==1
  		reg=[constant,trend,DUi];
    elseif model==2
    	reg=[constant,trend,DTi];
    elseif model==3
    	reg=[constant,DUi,trend,DTi];
    end
		
	%Estimation of alpha			
	khat=IC(y,reg,kmax,criteria);	
	u=(eye(T)-reg*inv(reg'*reg)*reg')*y; du=diff(u,1);
	depu=u; regu=lagn(u,1);
	j=1;
	while j<=khat-1
		regu=[regu,lagn(du,j)];
		j=j+1;
	end
	depu=trimr(depu,khat,0);
	regu=trimr(regu,khat,0);
	b=inv(regu'*regu)*(regu'*depu);
	ehat=depu-regu*b;
	VCV=(ehat'*ehat/rows(ehat))*inv(regu'*regu);
	ahat=b(1); vahat=VCV(1,1);
	tau1=(ahat-1)/sqrt(vahat);

	%Upper Biased Estimator
	if lam1<=0.1;                     t05=v_t(1);
	elseif (lam1>0.1)*(lam1<=0.2)==1; t05=v_t(2);
	elseif (lam1>0.2)*(lam1<=0.3)==1; t05=v_t(3);
	elseif (lam1>0.3)*(lam1<=0.4)==1; t05=v_t(4);
	elseif (lam1>0.4)*(lam1<=0.5)==1; t05=v_t(5);
	elseif (lam1>0.5)*(lam1<=0.6)==1; t05=v_t(6);
	elseif (lam1>0.6)*(lam1<=0.7)==1; t05=v_t(7);
	elseif (lam1>0.7)*(lam1<=0.8)==1; t05=v_t(8);
	elseif (lam1>0.8)*(lam1<=0.9)==1; t05=v_t(9);
	elseif lam1>0.9;                  t05=v_t(10);
	end
	IP=fix((khat+1)/2); r=cols(reg); k=10; c1=sqrt((1+r)*T);
	c2=((r+1)*T-t05*t05*(IP+T))/(t05*(t05+k)*(IP+T));
	if tau1>t05
		ctau=-tau1; 
	elseif (tau1<=t05)*(tau1>-k)==1
		ctau=IP*(tau1/T)-(r+1)/(tau1+c2*(tau1+k));
	elseif (tau1<=-k)*(tau1>-c1)==1
		ctau=IP*(tau1/T)-(r+1)/tau1;
	elseif tau1<=-c1
		ctau=0;
	end
	rhomd1 = ahat+ctau*sqrt(vahat);
	rhomd1= 1*(rhomd1>=1)+rhomd1*(abs(rhomd1)<1)-0.99*(rhomd1<=-1);
	amu=rhomd1;

	CR=(T^0.5)*abs(amu-1);	
	amus=amu*(CR>1)+1*(CR<=1);
	gdep=([(y(1)),(trimr(y,1,0)-amus*trimr(y,0,1))'])';
	greg=([(reg(1,:))',(trimr(reg,1,0)-amus*trimr(reg,0,1))'])';	
	b=inv(greg'*greg)*greg'*gdep;
	v=gdep-greg*b;
	

	if khat==1
		h0=v'*v/(rows(v));
	else
		ki=1;
		if amus==1
			while ki<=khat-1
				if ki==1
					regv=lagn(v,ki);
				else
					regv=[regv,lagn(v,ki)];
				end
				ki=ki+1;
			end			
			depv=trimr(v,khat-1,0);
			regv=trimr(regv,khat-1,0);		
			beta=inv(regv'*regv)*regv'*depv;
			e=depv-regv*beta;

            if model==1
				vbeta=zeros(3,khat-1);
				ki=1;
				while ki<=khat-1
					DUki=(trend>TBi-ki);
					regki=[constant,trend,DUki]; 
					gdepki=([(y(1)),(trimr(y,1,0)-amus*trimr(y,0,1))'])';
					gregki=([(reg(1,:))',(trimr(regki,1,0)-amus*trimr(regki,0,1))'])';	
					vbeta(:,ki)=inv(gregki'*gregki)*gregki'*gdepki;
					ki=ki+1;
				end			
				b(3)=b(3)-vbeta(3,:)*beta;
                h0=e'*e/(T-khat);
            elseif (model==2)
    			h0=(e'*e/(T-khat))/((1-beta'*ones(khat-1,1))^2);			
            elseif (model==3)
				vbeta=zeros(4,khat-1);
				ki=1;
				while ki<=khat-1
					DUki=(trend>TBi-ki);
					DTki=(trend>TBi-ki).*(trend-TBi);
					regki=[constant,DUki,trend,DTki]; 
					gdepki=([(y(1)),(trimr(y,1,0)-amus*trimr(y,0,1))'])';
					gregki=([(reg(1,:))',(trimr(regki,1,0)-amus*trimr(regki,0,1))'])';	
					vbeta(:,ki)=inv(gregki'*gregki)*gregki'*gdepki;
					ki=ki+1;
				end
                sige=e'*e/(T-khat);
    			h0=sige/((1-beta'*ones(khat-1,1))^2);			
				b(2)=sqrt(h0)*(b(2)-vbeta(2,:)*beta)/sqrt(sige);
            end   
		elseif abs(amus)<1
			[~,h0] = h0W(v);
		end
	end
	VCV=h0*inv(greg'*greg);
	vect1(TBi-fix(eps*T)+1)=(VR*b)'*inv(VR*VCV*VR')*(VR*b);
	TBi=TBi+1;
end

wald = log(sum(exp(vect1./2))/T);
TB = breakdate(y,model,kmax,eps);


% fprintf('*********************************************\n');
% if (model==1);
%     fprintf('Model 1: Structural change in intercept\n');
%     fprintf('Y{t}=a0+a1*DU+b0*t+e{t} where DU=1(t>TB)\n');
% elseif (model==2);
%     fprintf('Model 2: Structural change in slope\n');
%     fprintf('Y{t}=a0+b0*t+b1*DT+e{t} where DT=1(t>TB)*(t-TB)\n');
% elseif (model==3);
%     fprintf('Model 3: Structural change in both intercept and slope\n');
%     fprintf('Y{t}=a0+a1*DU+b0*t+b1*DT+e{t} where DU=1(t>TB) and DT=1(t>TB)*(t-TB)\n');
% end
% fprintf('\n');
% fprintf('THe EXP test Statistic (W-RQF)  =%6.4f\n', wald);	
    if eps==0.01; cv=cv(1,:);
    elseif eps==0.05; cv=cv(2,:);
    elseif eps==0.10; cv=cv(3,:);
    elseif eps==0.15; cv=cv(4,:);
    elseif eps==0.25; cv=cv(5,:);
    end
% display 'the critical value al 10%, 5% and 1% are:'
% fprintf('%6.4f %6.4f %6.4f\n',cv')
% fprintf('The estimate of the break date TB is: %6.4f\n', TB);
% fprintf('obtained by minimizing the sum of squared residuals from a regression of the relevant series on a constant, a time trend,\n');
% if (model==1);
%     fprintf('a level shift dummy.\n');
% elseif (model==2);
%     fprintf('a slope shift dummy.\n');
% elseif (model==3);
%     fprintf('a level shift and a slope shift dummies.\n');
% end
% fprintf('*********************************************\n');

%%%%%%%%%%%%%%%%%%FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%to take the n lag of x.
%x is a TxK matrix. return is also a TxK matrix with the first n elements 0.  
function y=lagn(x,n)
     if n>0
		y =[zeros(n,cols(x)) ; trimr(x,0,n)];
	else
		y =[trimr(x,abs(n),0) ; zeros(abs(n),cols(x))];
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%diff(x,k) is same as x{t}-x{t-k}
function y=diff(x,k)
	if k==0
		y = x;
	else
		y = [zeros(k,cols(x)) ; (trimr(x,k,0)-trimr(lagn(x,k),k,0))];
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AIC if criteria=1, BIC if criteria=2
function khat=IC(x,reg,kmax,criteria)
	ICV = zeros(kmax+1,1);
	dep = trimr(x,kmax,0);
	rego = trimr(reg,kmax,0);
	e=dep-rego*inv(rego'*rego)*rego'*dep;
	ICV(1) = log(e'*e/rows(e)); % in the case of k=0.
	i = 1;
	while i <= kmax
		reg = [reg,lagn(x,i)];
		rego = trimr(reg,kmax,0);
		[~,e,~] = olsqr2(dep,rego);
		if (criteria==1) %AIC%
			ICV(i+1) = log(e'*e/rows(e))+(2*i)/(rows(e));		
		elseif (criteria==2) %BIC%
			ICV(i+1) = log(e'*e/rows(e))+(log(rows(e))*i)/(rows(e));
		end
		i = i+1;
	end
	khat = minindc(ICV)-1;
	khat=khat*(khat>=1)+1*(khat==0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function is to take the autocovariance functions of x as
%ACF(x) = (R(0), R(1),...,R(T-1))'.
%Thus, R(1)=R(0), R(2)=R(1), and so on.
function R=ACV(x)
	T = rows(x);
	R = zeros(T,1); 
	xbar = mean(x); 
	i = 1;
	while i <= T
		j = i-1;
		R(i) = (x(1:T-j,1)-xbar)'*(x(1+j:T,1)-xbar)/T;
		i = i+1;	
	end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%---h(0)---%
%Andrews'(1991) method using a window:
%Quadratic Spectral window
function [m,h0]= h0W(x)
    T=rows(x);	
    b=olsqr(x(2:T),x(1:T-1));
	a=(4*(b^2))/((1-b)^4);
	R = ACV(x);
	lamda = zeros(T-1,1);
	s=seqa(1,1,T-1);

	m = 1.3221*((a*T)^(1/5));
	delta = (6*pi*s)/(5*m);
	i=1;
	while i<=T-1; 
		lamda(i) = 3*((sin(delta(i))/delta(i))-cos(delta(i)))/(delta(i)^2);
		i = i+1;
	end
	h0 = (R(1) + 2*lamda'*R(2:T,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TB= breakdate(x,model,kmax,eps)
    T=rows(x);
    constant=ones(T,1); 
    trend=seqa(1,1,T);
    TB0=max([fix(eps*T);(kmax+2)]);
    TBN=fix((1-eps)*T);

    vsum=ones(TBN-TB0+1,1)*1000;
    TBi=TB0;
    while TBi<=TBN
        DUi=(trend>TBi); DTi=(trend>TBi).*(trend-TBi); 
        if model==1
  		    reg=[constant,trend,DUi];
        elseif model==2
    	    reg=[constant,trend,DTi];
        elseif model==3
    	    reg=[constant,DUi,trend,DTi];
        end
		e=(eye(T)-reg*inv(reg'*reg)*reg')*x;
        vsum(TBi-TB0+1)=e'*e;
        TBi=TBi+1;
    end

    TB=TB0+minindc(vsum)-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b=olsqr(y,x)
b=inv(x'*x)*x'*y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ahat,r,fit]=olsqr2(y,x)
ahat=(x'*x)\(x'*y);
r=y-x*ahat;
fit=x*ahat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         AUXILIARY GAUSS TO MATLAB FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z=trimr(x,n1,n2)
% PURPOSE: return a matrix (or vector) x stripped of the specified rows.
% -----------------------------------------------------
% USAGE: z = trimr(x,n1,n2)
% where: x = input matrix (or vector) (n x k)
%       n1 = first n1 rows to strip
%       n2 = last  n2 rows to strip
% NOTE: modeled after Gauss trimr function
% -----------------------------------------------------
% RETURNS: z = x(n1+1:n-n2,:)
% -----------------------------------------------------

% written by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jpl%jpl.econ.utoledo.edu

  [n, ~] = size(x);
  if (n1+n2) >= n; 
     error('Attempting to trim too much in trimr');
  end;
  h1 = n1+1;   
  h2 = n-n2;
  z = x(h1:h2,:);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function seq=seqa(a,b,c)
% PURPOSE: produce a sequence of values
% -----------------------------------------------------
% USAGE: y = seqa(a,b,c)
%  where    a = initial value in sequence 
%           b = increment
%           c = number of values in the sequence  
% -----------------------------------------------------
% RETURNS: a sequence, (a:b:(a+b*(c-1)))' in MATLAB notation
% ----------------------------------------------------- 
% NOTE: a Gauss compatability function
% -----------------------------------------------------

% written by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jpl@jpl.econ.utoledo.edu
       
% seqa Gauss eqivalent of seqa(a,b,c)
seq=(a:b:(a+b*(c-1)))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = diagrv(x,v)
% PURPOSE: replaces main diagonal of a square matrix
% -----------------------------------------
% USAGE: y - diagrv(x,v)
% where: x = input matrix
%        v = vector to replace main diagonal
% -----------------------------------------
% RETURNS: y = matrix x with v placed on main diagonal
% -----------------------------------------
% NOTE: a Gauss compatability function
% -----------------------------------------------------

% written by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jpl@jpl.econ.utoledo.edu
  
[r,c] = size(x);
if r ~= c;
  error('x matrix not square')
end;
if length(v) ~= r;
  error('v is not conformable with x')
end;
y = x - diag(diag(x)) + diag(v);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c=minindc(A)
    [~,n]=size(A);
    c=zeros(n,1);
    for i=1:n
        [~,index]=min(A(:,i));
        c(i)=index;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function re=rev(x)
% REV 
% REV(x) reverses the elements of x
% It is equivalent to the Matlab functin flipud
d=size(x); li=d(1,1);
ind=li:-1:1;
re=x(ind,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m=rows(A)
[m,~]=size(A);
            
function n=cols(A)
[~,n]=size(A);