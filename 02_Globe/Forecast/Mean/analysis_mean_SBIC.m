% ========================================================================
% analysis_mean_SBIC.m
% ------------------------------------------------------------------------
% Analysis of SBIC selection and long-run forecasts for the global mean.
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

FUN_GENERAL_DIR  = fullfile(ROOT_DIR,'functions');
FUN_FORECAST_DIR = fullfile(ROOT_DIR,'functions_for');

FIG_DIR = fullfile(THIS_DIR,'Figures');
TAB_DIR = fullfile(THIS_DIR,'Tables');
RES_DIR = fullfile(THIS_DIR,'Results');


if ~isfolder(FIG_DIR), mkdir(FIG_DIR); end
if ~isfolder(TAB_DIR), mkdir(TAB_DIR); end
if ~isfolder(RES_DIR), mkdir(RES_DIR); end

addpath(genpath(FUN_GENERAL_DIR));
addpath(genpath(FUN_FORECAST_DIR));

% ------------------------------------------------------------------------
% 2. Load data
% ------------------------------------------------------------------------

DATA_FILE = fullfile(ROOT_DIR, ...
    '02_Globe','Introduction','Results', ...
    'QUANTILES_monthly_Globe_1880_2023.mat');

assert(isfile(DATA_FILE), ...
    'Input file QUANTILES_monthly_Globe_1880_2023.mat not found.');

load(DATA_FILE,'QUANTILES_monthly');

% ------------------------------------------------------------------------
% 3. Global mean series
% ------------------------------------------------------------------------

y = QUANTILES_monthly.Globe(:,1);
t = length(y);

m      = 14;
pmax   = 12;
maxpol = 12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Select model with SBIC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MEAN Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y_hat = mean(y);
res   = y-y_hat;
[aic,sbic,hqc] = information_criteria(res,y,0);
SBIC(1) = sbic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LINEAR-TREND Model%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1;
[beta,t_ratio,y_hat,varres,varbeta] = estima_pol_trend_hac(y,k);
res = y-y_hat;
[aic,sbic,hqc] = information_criteria(res,y,2);
SBIC(2) = sbic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%POL-TREND Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = estima_select_pol_trend_hac_sbic(y,maxpol);
[beta,t_ratio,y_hat] = estima_pol_trend_hac(y,k);
res = y-y_hat;
[aic,sbic,hqc] = information_criteria(res,y,k+1);
SBIC(3) = sbic;
SBIC(4) = sbic; % slope average model

%%%%%%%%%%%%%%%%%%%%%%%%%%%POL-TREND-LOG Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = estima_select_pol_trend_log_hac_sbic(y,maxpol);
[beta,t_ratio,y_hat] = estima_pol_trend_log_hac(y,k);
res = y-y_hat;
[aic,sbic,hqc] = information_criteria(res,y,k+1);
SBIC(5) = sbic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%STRUCT-BREAK Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = 3;
criteria = 2;
eps = 0.15;
kmax = fix(12*(t/100)^(1/4));

[wald,cv,tb] = qfgls_simus(y,kmax,model,criteria,eps);

if wald > cv(2)
    [beta_break,t_nw,se_nw,res,r2,varres,varBhat,y_hat] = ...
        estima_break_model(y,tb,model);
    [aic,sbic,hqc] = information_criteria(res,y,k+1);
else
    sbic = SBIC(2);
end

SBIC(6) = sbic;

%%%%%%%%%%%%%%%%%%%%%%%%%Pol-trend+AR(p) Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[p,k,beta,t_ratio,y_hat,res] = ...
    estima_select_pol_trend_arp_hac_sbic(y,maxpol,pmax);

[aic,sbic,hqc] = information_criteria(res,y,1+p+k);
SBIC(7) = sbic;
SBIC(8) = sbic; % slope average model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%AR(p) Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[beta,t_ratio,p,y_hat,res] = estima_select_arp_sbic(y,pmax);
[aic,sbic,hqc] = information_criteria(res,y,1+p);
SBIC(9) = sbic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RW model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res = diff(y,1);
[aic,sbic,hqc] = information_criteria(res,y,0);
SBIC(10) = sbic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RW + drift model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[b,t_nw,se_nw,res] = ...
    ols_hac_forecast(diff(y,1),ones(t-1,1));

[aic,sbic,hqc] = information_criteria(res,y,1);
SBIC(11) = sbic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%IMA model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = y;
dx = diff(x,1);
mu = mean(dx);

xm = x-mean(x);
dxm = diff(xm,1);

ima = armax(dxm,[0 1]);

theta = -ima.c(2);
sigma2 = ima.Report.Fit.MSE;

SBIC(12) = log(sigma2)+log(t)*1/t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Arfima Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ar_part = 0;
[d,ar,sigma2,res,y_hat] = estima_arfima(y,ar_part);
[aic,sbic,hqc] = information_criteria(res,y,2+length(ar));
SBIC(13) = sbic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%AR20(p) Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 20;
[beta,t_ratio,p,y_hat,res] = estima_arp(y,p);
[aic,sbic,hqc] = information_criteria(res,y,1+p);
SBIC(14) = sbic;

% ------------------------------------------------------------------------
% 4. BIC weights and LaTeX table
% ------------------------------------------------------------------------

W_BIC = exp(-0.5*SBIC) ./ sum(exp(-0.5*SBIC));

model_labels = { ...
    'mean', ...
    'linear-trend', ...
    'pol-trend (k=5)', ...
    'pol-trend-av-sl (k=5)', ...
    'pol-trend-log (k=8)', ...
    'struct-breaks', ...
    'pol-trend-arp (k=3, p=1)', ...
    'pol-trend-arp-av-sl (k=3, p=1)', ...
    'arp (p=2)', ...
    'rw', ...
    'rwd', ...
    'ima', ...
    'arfima', ...
    'arp20'};

TAB_FILE = fullfile(TAB_DIR,'Table_mean_SBIC_Globe_1880_2023.tex');

fid = fopen(TAB_FILE,'w');

fprintf(fid,'\\begin{tabular}{lrr}\n');
fprintf(fid,'\\toprule\n');
fprintf(fid,'Models & BIC & weights \\\\\n');
fprintf(fid,'\\midrule\n');

for ii = 1:length(model_labels)
    fprintf(fid,'%s & %.4f & %.4f \\\\\n', ...
        model_labels{ii}, SBIC(ii), W_BIC(ii));
end

fprintf(fid,'\\bottomrule\n');
fprintf(fid,'\\end{tabular}\n');

fclose(fid);

% ------------------------------------------------------------------------
% 5. Forecast table for selected models
% ------------------------------------------------------------------------

F = [1,10,25,50,100];

FORECAST_table = NaN(length(F),3);
CI_table       = NaN(length(F),3,2);

% Benchmark model: linear trend
k_linear = 1;
[beta_linear,~,y_hat_linear,varres_linear,~] = ...
    estima_pol_trend_hac(y,k_linear);

% Selected model: polynomial trend
k_pol = estima_select_pol_trend_hac_sbic(y,maxpol);
[beta_pol,~,y_hat_pol,varres_pol,~] = ...
    estima_pol_trend_hac(y,k_pol);

slope_pol = compute_slope(beta_pol,t);

for ii = 1:length(F)
    
    f = F(ii);
    
    % Benchmark: linear trend
    x_linear = [1,t+f];
    xf_linear = beta_linear' * x_linear';
    ci_linear = xf_linear + norminv([0.05 0.95])*sqrt(varres_linear);
    
    % Selected model: polynomial trend
    x_pol = [];
    for jj = 0:k_pol
        x_pol = [x_pol,(t+f)^jj];
    end
    
    xf_pol = beta_pol' * x_pol';
    ci_pol = xf_pol + norminv([0.05 0.95])*sqrt(varres_pol);
    
    % Selected model: polynomial trend with average slope
    beta_avsl = [beta_pol(1),slope_pol];
    x_avsl = [1,t+f];
    
    xf_avsl = beta_avsl * x_avsl';
    ci_avsl = xf_avsl + norminv([0.05 0.95])*sqrt(varres_pol);
    
    FORECAST_table(ii,:) = [xf_linear,xf_pol,xf_avsl];
    CI_table(ii,1,:) = ci_linear;
    CI_table(ii,2,:) = ci_pol;
    CI_table(ii,3,:) = ci_avsl;
    
end

TAB_FILE2 = fullfile(TAB_DIR,'Table_mean_forecast_Globe_1880_2023.tex');

fid = fopen(TAB_FILE2,'w');

fprintf(fid,'\\begin{tabular}{lrrr}\n');
fprintf(fid,'\\toprule\n');
fprintf(fid,'Horizon & Benchmark model: linear trend & Selected model: pol-trend & Selected model: pol-trend-av-sl \\\\\n');
fprintf(fid,'\\midrule\n');

for ii = 1:length(F)
    
    fprintf(fid,'%d & %.2f & %.2f & %.2f \\\\\n', ...
        F(ii), FORECAST_table(ii,1), FORECAST_table(ii,2), FORECAST_table(ii,3));
    
    fprintf(fid,' & (%.2f,%.2f) & (%.2f,%.2f) & (%.2f,%.2f) \\\\\n', ...
        CI_table(ii,1,1), CI_table(ii,1,2), ...
        CI_table(ii,2,1), CI_table(ii,2,2), ...
        CI_table(ii,3,1), CI_table(ii,3,2));
    
end

fprintf(fid,'\\bottomrule\n');
fprintf(fid,'\\end{tabular}\n');

fclose(fid);

% ------------------------------------------------------------------------
% 6. Save local RDOS structure (compatibility)
% ------------------------------------------------------------------------

RDOS.SBIC    = SBIC;
RDOS.weights = W_BIC;
RDOS.models  = model_labels;
RDOS.forecast_table = FORECAST_table;
RDOS.forecast_ci    = CI_table;
RDOS.forecast_H     = F;

save(fullfile(RES_DIR,'RDOS_for_mean_Globe_1880_2023.mat'),'RDOS');

fprintf('\nMean SBIC analysis completed.\n');
fprintf('Tables saved in:\n%s\n', TAB_DIR);
