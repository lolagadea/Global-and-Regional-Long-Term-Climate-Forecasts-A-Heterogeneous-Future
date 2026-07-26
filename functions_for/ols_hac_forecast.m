function[b,t_nw, se_nw, res, r2, varres, varBhat, y_hat]=ols_hac_forecast(y,X)
% Forecast-specific HAC routine.
% This version preserves the output interface used in the accepted-paper
% forecast codes. It is kept separate from the general ols_hac function to
% avoid changing results in the replication package.

 % calculate Bhat
  b = pinv(X)*y;[Q, R]=qr(X,0); b= R\(Q'*y);
  %bb=inv(X'*X)*X'*y;
  % Determine the size of the matrix of regressors
  [T, k] =size(X);
  L = round(4*((T/100)^(2/9)));
  %L=T^(1/3);
 
  % Generate residuals
  e = y - X * b;
 
  % Calculate the Newey-West autocorrelation consistent covariance
  % estimator. See p.464 of Greene (2000) for more information
  % Note that there was a typo in the 2000 edition that was corrected
  % in the 2003 edition.
  Q = 0;
  for l = 0:L
    w_l = 1-l/(L+1);
    for t = l+1:T
      if (l==0)   % This calculates the S_0 portion
        Q = Q  + e(t) ^2 * X(t, :)' * X(t,:);
      else        % This calculates the off-diagonal terms
        Q = Q + w_l * e(t) * e(t-l)* ...
          (X(t, :)' * X(t-l,:) + X(t-l, :)' * X(t,:));
      end
    end
  end
  Q = 1/(T-k) * Q;
  
  % Calculate Newey-White standard errors
  varBhat = T * inv(X' * X) * Q * inv(X' * X);

  % calculate standard errors and t-stats
  se_nw = sqrt(diag(varBhat));
  t_nw = b ./ se_nw;
  
  %standard errors normales:
   var_beta=inv(X'*X)*(e'*e)/(T-k);
   se_beta=diag(var_beta).^0.5;
   t_ratio=b./se_beta;
   varres=(e'*e)/(T-k);
  
   %Residuals, fitted values and r2
   res=e; 
   y_hat=X*b;
   r2=1-e'*e/((y-mean(y))'*(y-mean(y)));
  
   
    

