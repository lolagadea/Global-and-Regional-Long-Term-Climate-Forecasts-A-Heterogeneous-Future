% ========================================================================
% compute_weights_BIC.m
% ------------------------------------------------------------------------
% Compute BIC-based model weights for the global forecast competition.
%
% Outputs:
%   Results_for/Weights_BIC.dat
%   Results_for/BIC.dat
% ========================================================================

clear; clc;
warning('off');

% ------------------------------------------------------------------------
% 1. Portable paths
% ------------------------------------------------------------------------

THIS_DIR = fileparts(mfilename('fullpath'));
FORECAST_DIR = fileparts(THIS_DIR);
GLOBE_DIR = fileparts(FORECAST_DIR);
ROOT_DIR = fileparts(GLOBE_DIR);

FUN_GENERAL_DIR  = fullfile(ROOT_DIR, 'functions');
FUN_FORECAST_DIR = fullfile(ROOT_DIR, 'functions_for');

RES_DIR = fullfile(THIS_DIR, 'Results_for');

if ~isfolder(RES_DIR), mkdir(RES_DIR); end

addpath(genpath(FUN_GENERAL_DIR));
addpath(genpath(FUN_FORECAST_DIR));

% ------------------------------------------------------------------------
% 2. Load data
% ------------------------------------------------------------------------

DATA_FILE = fullfile(ROOT_DIR, ...
    '02_Globe', 'Introduction', 'Results', ...
    'QUANTILES_monthly_Globe_1880_1960.mat');

assert(isfile(DATA_FILE), ...
    'Input file QUANTILES_monthly_Globe_1880_1960.mat not found.');

load(DATA_FILE, 'QUANTILES_monthly');

Y = QUANTILES_monthly.Globe;

[t,n] = size(Y);

m = 14;
pmax = 12;
maxpol = 12;

SBIC = NaN(n,m);

% ------------------------------------------------------------------------
% 3. Compute BIC for each model
% ------------------------------------------------------------------------

for i = 1:n
    
    y = Y(:,i);
    
    % Mean model
    y_hat = mean(y);
    res = y - y_hat;
    [~,sbic,~] = information_criteria(res,y,0);
    SBIC(i,1) = sbic;
    
    % Linear trend model
    k = 1;
    [~,~,y_hat,~,~] = estima_pol_trend_hac(y,k);
    res = y - y_hat;
    [~,sbic,~] = information_criteria(res,y,2);
    SBIC(i,2) = sbic;
    
    % Polynomial trend model
    k = estima_select_pol_trend_hac_sbic(y,maxpol);
    [~,~,y_hat] = estima_pol_trend_hac(y,k);
    res = y - y_hat;
    [~,sbic,~] = information_criteria(res,y,k+1);
    SBIC(i,3) = sbic;
    SBIC(i,4) = sbic; % average-slope version
    
    % Polynomial trend in logs
    k = estima_select_pol_trend_log_hac_sbic(y,maxpol);
    [~,~,y_hat] = estima_pol_trend_log_hac(y,k);
    res = y - y_hat;
    [~,sbic,~] = information_criteria(res,y,k+1);
    SBIC(i,5) = sbic;
    
    % Structural break model
    model = 3;
    criteria = 2;
    eps = 0.15;
    kmax = fix(12*(t/100)^(1/4));
    
    [wald,cv,tb] = qfgls_simus(y,kmax,model,criteria,eps);
    
    if wald > cv(2)
        [~,~,~,res,~,~,~,~] = estima_break_model(y,tb,model);
        [~,sbic,~] = information_criteria(res,y,k+1);
    else
        sbic = SBIC(i,2);
    end
    
    SBIC(i,6) = sbic;
    
    % Polynomial trend + AR(p)
    [p,k,~,~,~,res] = estima_select_pol_trend_arp_hac_sbic(y,maxpol,pmax);
    [~,sbic,~] = information_criteria(res,y,1+p+k);
    SBIC(i,7) = sbic;
    SBIC(i,8) = sbic; % average-slope version
    
    % AR(p)
    [~,~,p,~,res] = estima_select_arp_sbic(y,pmax);
    [~,sbic,~] = information_criteria(res,y,1+p);
    SBIC(i,9) = sbic;
    
    % Random walk
    res = diff(y,1);
    [~,sbic,~] = information_criteria(res,y,0);
    SBIC(i,10) = sbic;
    
    % Random walk with drift
    dx = diff(y,1);
    reg = ones(t-1,1);
    b = ols_hac(dx,reg);
    res = dx - reg*b;
    [~,sbic,~] = information_criteria(res,y,1);
    SBIC(i,11) = sbic;
    
    % IMA model
    x = y;
    dx = diff(x,1);
    xm = x - mean(x);
    dxm = diff(xm,1);
    ima = armax(dxm,[0 1]);
    sigma2 = ima.Report.Fit.MSE;
    SBIC(i,12) = log(sigma2) + log(t)*1/t;
    
    % ARFIMA model
    ar_part = 0;
    [~,ar,~,res,~] = estima_arfima(y,ar_part);
    [~,sbic,~] = information_criteria(res,y,2+length(ar));
    SBIC(i,13) = sbic;
    
    % AR(20) model
    p = 20;
    [~,~,p,~,res] = estima_arp(y,p);
    [~,sbic,~] = information_criteria(res,y,1+p);
    SBIC(i,14) = sbic;
    
end

% ------------------------------------------------------------------------
% 4. Compute BIC weights
% ------------------------------------------------------------------------

W = NaN(n,m);

for i = 1:n
    for j = 1:m
        W(i,j) = exp(-0.5*SBIC(i,j)) / sum(exp(-0.5*SBIC(i,:)));
    end
end

% ------------------------------------------------------------------------
% 5. Save outputs
% ------------------------------------------------------------------------

save(fullfile(RES_DIR,'Weights_BIC.dat'), 'W', '-ASCII');
save(fullfile(RES_DIR,'BIC.dat'), 'SBIC', '-ASCII');

fprintf('\nBIC weights completed.\n');
fprintf('Results saved in:\n%s\n', RES_DIR);