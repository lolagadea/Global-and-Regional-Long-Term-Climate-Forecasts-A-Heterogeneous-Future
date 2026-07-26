% ========================================================================
% compute_weights_BIC_Australia.m
% ------------------------------------------------------------------------
% Compute BIC weights for Australia temperature characteristics.
% Replication package IJF - Regional Forecast module.
%
% Inputs:
%   03_Regional/Introduction/Results/QUANTILES_monthly_1960_2023.mat
%
% Outputs:
%   03_Regional/Forecast/Australia/Results_for/Weights_BIC_Australia.dat
%   03_Regional/Forecast/Australia/Results_for/BIC_Australia.dat
% ========================================================================

clear; clc;
clear global;
warning('off');

restoredefaultpath;
rehash toolboxcache;

% ------------------------------------------------------------------------
% 1. Portable paths
% ------------------------------------------------------------------------

THIS_DIR     = fileparts(mfilename('fullpath'));
FORECAST_DIR = fileparts(THIS_DIR);
REGIONAL_DIR = fileparts(FORECAST_DIR);
ROOT_DIR     = fileparts(REGIONAL_DIR);

FUN_DIR     = fullfile(ROOT_DIR,'functions');
FUN_FOR_DIR = fullfile(ROOT_DIR,'functions_for');

addpath(genpath(FUN_DIR));
addpath(genpath(FUN_FOR_DIR));

INTRO_RES_DIR = fullfile(REGIONAL_DIR,'Introduction','Results');
OUT_RES_DIR   = fullfile(THIS_DIR,'Results_for');

if ~isfolder(OUT_RES_DIR)
    mkdir(OUT_RES_DIR);
end

% ------------------------------------------------------------------------
% 2. Load regional quantiles
% ------------------------------------------------------------------------

QUANT_FILE = fullfile(INTRO_RES_DIR,'QUANTILES_monthly_1960_2023.mat');

assert(exist(QUANT_FILE,'file')==2, ...
    'QUANTILES_monthly_1960_2023.mat not found in Regional/Introduction/Results.');

load(QUANT_FILE,'QUANTILES_monthly');

Y = QUANTILES_monthly.Australia;

[t,n] = size(Y);

% ------------------------------------------------------------------------
% 3. Settings
% ------------------------------------------------------------------------

m      = 14;
pmax   = 12;
maxpol = 12;

SBIC = NaN(n,m);

% ------------------------------------------------------------------------
% 4. Compute BIC by characteristic and model
% ------------------------------------------------------------------------

for i = 1:n

    y = Y(:,i);

    % 1. Mean model
    y_hat = mean(y);
    res = y - y_hat;
    [aic,sbic,hqc] = information_criteria(res,y,0);
    SBIC(i,1) = sbic;

    % 2. Linear trend model
    k = 1;
    [beta,t_ratio,y_hat,varres,varbeta] = estima_pol_trend_hac(y,k);
    res = y - y_hat;
    [aic,sbic,hqc] = information_criteria(res,y,2);
    SBIC(i,2) = sbic;

    % 3. Polynomial trend model
    k = estima_select_pol_trend_hac_sbic(y,maxpol);
    [beta,t_ratio,y_hat] = estima_pol_trend_hac(y,k);
    res = y - y_hat;
    [aic,sbic,hqc] = information_criteria(res,y,k+1);
    SBIC(i,3) = sbic;
    SBIC(i,4) = sbic;      % slope average model

    % 5. Polynomial trend log model
    k = estima_select_pol_trend_log_hac_sbic(y,maxpol);
    [beta,t_ratio,y_hat] = estima_pol_trend_log_hac(y,k);
    res = y - y_hat;
    [aic,sbic,hqc] = information_criteria(res,y,k+1);
    SBIC(i,5) = sbic;

    % 6. Structural breaks model
    model = 3;
    criteria = 2;          % BIC if criteria = 2, AIC if criteria = 1
    eps = 0.15;
    kmax = fix(12*(t/100)^(1/4));

    [wald,cv,tb] = qfgls_simus(y,kmax,model,criteria,eps);

    if wald > cv(2)
        [beta_break,t_nw,se_nw,res,r2,varres,varBhat,y_hat] = ...
            estima_break_model(y,tb,model);
        [aic,sbic,hqc] = information_criteria(res,y,k+1);
    else
        sbic = SBIC(i,2);
    end

    SBIC(i,6) = sbic;

    % 7. Polynomial trend + AR(p) model
    [p,k,beta,t_ratio,y_hat,res] = ...
        estima_select_pol_trend_arp_hac_sbic(y,maxpol,pmax);

    [aic,sbic,hqc] = information_criteria(res,y,1+p+k);
    SBIC(i,7) = sbic;
    SBIC(i,8) = sbic;      % slope average model

    % 9. AR(p) model
    [beta,t_ratio,p,y_hat,res] = estima_select_arp_sbic(y,pmax);
    [aic,sbic,hqc] = information_criteria(res,y,1+p);
    SBIC(i,9) = sbic;

    % 10. Random walk model
    res = diff(y,1);
    [aic,sbic,hqc] = information_criteria(res,y,0);
    SBIC(i,10) = sbic;

    % 11. Random walk with drift model
    [b,t_nw,se_nw,res] = ols_hac_forecast(diff(y,1),ones(t-1,1));
    [aic,sbic,hqc] = information_criteria(res,y,1);
    SBIC(i,11) = sbic;

    % 12. IMA model
    x = y;
    dx = diff(x,1);
    mu = mean(dx);
    xm = x - mean(x);
    dxm = diff(xm,1);

    ima = armax(dxm,[0 1]);
    theta = -ima.c(2);
    sigma2 = ima.Report.Fit.MSE;

    SBIC(i,12) = log(sigma2) + log(t)*1/t;

    % 13. ARFIMA model
    ar_part = 0;
    [d,ar,sigma2,res,y_hat] = estima_arfima(y,ar_part);
    [aic,sbic,hqc] = information_criteria(res,y,2+length(ar));
    SBIC(i,13) = sbic;

    % 14. AR(20) model
    p = 20;
    [beta,t_ratio,p,y_hat,res] = estima_arp(y,p);
    [aic,sbic,hqc] = information_criteria(res,y,1+p);
    SBIC(i,14) = sbic;

end

% ------------------------------------------------------------------------
% 5. Compute BIC weights
% ------------------------------------------------------------------------

W = NaN(n,m);

for i = 1:n
    for j = 1:m
        W(i,j) = exp(-0.5*SBIC(i,j)) / sum(exp(-0.5*SBIC(i,:)));
    end
end

% ------------------------------------------------------------------------
% 6. Save outputs
% ------------------------------------------------------------------------

save(fullfile(OUT_RES_DIR,'Weights_BIC_Australia.dat'),'W','-ASCII');
save(fullfile(OUT_RES_DIR,'BIC_Australia.dat'),'SBIC','-ASCII');