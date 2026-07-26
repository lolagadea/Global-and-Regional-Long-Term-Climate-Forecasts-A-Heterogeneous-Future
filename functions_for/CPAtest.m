function  [teststat,critval,pval,sign]=CPAtest(loss1,loss2,tau, alpha, choice)

% This function performs the asymptotic Conditional Predictive Ability Test
% INPUTS: loss1 and loss2, Tx1 vectors of losses over the out of sample period for the two models under consideration
%         tau, the forecast horizon
%         alpha, niminal risk level: 1%, 5%, 10%
%         choice: 1 if unconditional ; 2 if conditional
%
% OUTPUTS: teststat, the test statistic of the conditional predictive ability test
%          critval, the critical value of the test for a 5% level (the test is a chi square)
%          pval, the p-value of the test
%
% Raffaella Giacomini, 2003

lossdiff1=loss1-loss2;  
% loss diferential 

TT = size(lossdiff1,1);

if choice==1
    instruments=ones(TT,1);
    lossdiff=lossdiff1;
    T=TT;
else    
    instruments=[ones(TT-tau,1) lossdiff1(1:end-tau)];
    lossdiff=lossdiff1(tau+1:end);                                         %loss diferential
    T=TT-tau;
end      

% create the regressor matrix given by lossdiff*ht', where ht is the matrix of instruments
reg = -999*ones(size(instruments));
for jj = 1:size(instruments,2)
   reg(:,jj) = instruments(:,jj).*lossdiff;
end
if tau == 1
   % calculate the test stat as nR^2 from the regression of one on lossdiff*ht
   res.beta = reg\ones(T,1);   
   err = ones(T,1)-reg*res.beta;
   r2 = 1-mean(err.^2);
   teststat = T*r2;
   q = size(reg,2);
   critval = chi2inv(1-alpha,q);
   pval = 1 - cdf('chi2',abs(teststat),q);
else
   zbar = mean(reg)';
   nlags = tau-1;
   omega = NeweyWest(reg,nlags);
   teststat = T*zbar'*inv(omega)*zbar;
   q = size(reg,2);
   critval = chi2inv(1-alpha,q);
   pval = 1 - cdf('chi2',abs(teststat),q);
end

av_diff_loss=mean(loss1-loss2);

if av_diff_loss<0
    sign=0; %negative, the first is better
elseif av_diff_loss>0
    sign=1; %positive, the second is better
end

function omegahat = NeweyWest(Z,nlags)

% Returns the Newey-West estimator of the asymptotic variance matrix
% INPUTS: Z, a nxk matrix with rows the vector zt'
%         nlags, the number of lags
%
% OUTPUTS: omegahat, the Newey-West estimator of the covariance matrix


[n,k] = size(Z);

% de-mean the variables
Z = Z - ones(size(Z,1),1)*mean(Z);

gamma = -999*ones(nlags,k);
samplevar = Z'*Z/n; % sample variance
omegahat = samplevar;
if nlags > 0
   % sample autocovariances
   for ii = 1:nlags
      Zlag = [zeros(ii,k);Z(1:n-ii,:)];
      gamma = (Z'*Zlag +Zlag'*Z)/n;
      weights = 1 - (ii/(nlags+1));
      omegahat = omegahat + weights*gamma;
   end
end



   
   