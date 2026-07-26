function [dfa,p] = calcula_adf_t_sbic(y,pmax)
% =========================================================================
% calcula_adf_t_sbic.m
% -------------------------------------------------------------------------
% Purpose:
%   Compute the Augmented Dickey-Fuller statistic with deterministic trend,
%   selecting the lag length by the Schwarz Bayesian Information Criterion
%   (SBIC).
%
% Usage:
%   [dfa,p] = calcula_adf_t_sbic(y,pmax)
%
% Inputs:
%   y     : univariate time series
%   pmax  : maximum number of ADF augmentation lags
%
% Outputs:
%   dfa   : ADF t-statistic for the lagged level term
%   p     : selected lag length
%
% Notes:
%   - The ADF regression includes:
%       constant, linear trend, lagged level, and p lagged differences.
%   - The original econometric logic is preserved.
% =========================================================================

% Ensure column vector
y = y(:);

sbic = NaN(pmax+1,1);

for i = 0:pmax

    [~,~,res] = ols_adf_t(y,i);

    T = length(res);
    sigma2_res = res' * res / T;

    % SBIC for ADF regression with constant, trend, lagged level, and i lags
    sbic(i+1) = log(sigma2_res) + (i+3) * log(T) / T;

end

[~,p] = min(sbic);

% Correction because MATLAB indices start at 1 and lag length starts at 0
p = p - 1;

[beta,sigma,~] = ols_adf_t(y,p);

% ADF statistic: t-ratio of the lagged level coefficient
dfa = beta(3) / sigma(3);

end


% =========================================================================
% Local function: ols_adf_t
% -------------------------------------------------------------------------
% Estimate ADF regression with deterministic trend.
% =========================================================================

function [beta,sigma,res] = ols_adf_t(y,p)

n = length(y);

dy = NaN(n-1,1);

for j = 1:n-1
    dy(j) = y(j+1) - y(j);
end

x = zeros(n-1-p,p);

% Lagged differences
for j = 1:p
    x(:,j) = dy(p-j+1:n-1-j);
end

trend = (1:n-1-p)';

% Regressors:
%   constant, trend, lagged level, lagged differences
x = [ones(n-1-p,1), trend, y(p+1:n-1), x];

[t,k] = size(x);

beta = ols(dy(p+1:n-1),x);

res = dy(p+1:n-1) - x*beta;

sigma2 = inv(x'*x) * (res'*res) / (t-k);

sigma = sqrt(diag(sigma2));

end


% =========================================================================
% Local function: ols
% -------------------------------------------------------------------------
% Compute OLS estimates using QR decomposition.
% =========================================================================

function beta = ols(y,X)

[Q,R] = qr(X,0);
beta = R \ (Q' * y);

end