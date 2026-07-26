function [beta,t_ratio,R2,R2bar,F,varres,varbeta] = ols_hac(y,X)
% =========================================================================
% ols_hac.m
% -------------------------------------------------------------------------
% Purpose:
%   Estimate an OLS regression with HAC (heteroskedasticity and
%   autocorrelation consistent) standard errors using the
%   Newey-West estimator.
%
% Usage:
%   [beta,t_ratio,R2,R2bar,F,varres,varbeta] = ols_hac(y,X)
%
% Inputs:
%   y : dependent variable (T x 1)
%   X : regressor matrix   (T x k)
%
% Outputs:
%   beta     : OLS coefficient estimates
%   t_ratio  : HAC t-ratios
%   R2       : coefficient of determination
%   R2bar    : adjusted R-squared
%   F        : F statistic
%   varres   : residual variance
%   varbeta  : HAC covariance matrix of beta
%
% Notes:
%   - HAC covariance matrix is computed using the Newey-West estimator.
%   - The truncation lag follows the original implementation:
%
%         p = floor(4*(T/100)^(2/9))
%
%   - The original econometric implementation is preserved.
%
% =========================================================================

% Ensure column vector
y = y(:);

[T,k] = size(X);

% -------------------------------------------------------------------------
% 1. OLS estimation
% -------------------------------------------------------------------------

beta = (X' * X) \ (X' * y);

yhat = X * beta;
res  = y - yhat;

% -------------------------------------------------------------------------
% 2. Goodness-of-fit statistics
% -------------------------------------------------------------------------

SSR = res' * res;
SST = (y - mean(y))' * (y - mean(y));

R2 = 1 - SSR/SST;

R2bar = 1 - ((T-1)/(T-k)) * (1-R2);

varres = SSR / (T-k);

% -------------------------------------------------------------------------
% 3. HAC covariance matrix (Newey-West)
% -------------------------------------------------------------------------

% Original truncation lag
p = floor(4 * (T/100)^(2/9));

XX = (X' * X) / T;

S0 = zeros(k,k);

for t = 1:T
    S0 = S0 + (res(t)^2) * (X(t,:)' * X(t,:));
end

S0 = S0 / T;

S = S0;

for j = 1:p

    Sj = zeros(k,k);

    for t = j+1:T

        Sj = Sj + ...
            res(t) * res(t-j) * ...
            (X(t,:)' * X(t-j,:));

    end

    Sj = Sj / T;

    weight = 1 - j/(p+1);

    S = S + weight * (Sj + Sj');

end

varbeta = (1/T) * inv(XX) * S * inv(XX);

% -------------------------------------------------------------------------
% 4. HAC t-ratios
% -------------------------------------------------------------------------

se = sqrt(diag(varbeta));

t_ratio = beta ./ se;

% -------------------------------------------------------------------------
% 5. F statistic
% -------------------------------------------------------------------------

if k > 1

    R = [zeros(k-1,1), eye(k-1)];

    F = (R*beta)' * inv(R*varbeta*R') * (R*beta) / (k-1);

else

    F = NaN;

end

end