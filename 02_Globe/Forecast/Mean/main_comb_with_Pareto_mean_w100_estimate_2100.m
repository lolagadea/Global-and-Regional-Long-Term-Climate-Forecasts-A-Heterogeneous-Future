% ========================================================================
% main_comb_with_Pareto_mean_w100_estimate_2100.m
% ------------------------------------------------------------------------
% Forecast combinations for the global mean up to year 2100.
% Model selection is based on Pareto-superior models for h=25, w=100.
% ========================================================================

clear; clc;
clear global;
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

COMP_RES_DIR = fullfile(FORECAST_DIR, ...
    'run_forecast_competition','Results_for');

TAB_DIR = fullfile(THIS_DIR,'Tables');
RES_DIR = fullfile(THIS_DIR,'Results');

if ~isfolder(TAB_DIR), mkdir(TAB_DIR); end
if ~isfolder(RES_DIR), mkdir(RES_DIR); end

addpath(genpath(FUN_GENERAL_DIR));
addpath(genpath(FUN_FORECAST_DIR));
addpath(genpath(COMP_RES_DIR));

% ------------------------------------------------------------------------
% 2. Load inputs
% ------------------------------------------------------------------------

DATA_FILE = fullfile(ROOT_DIR, ...
    '02_Globe','Introduction','Results', ...
    'QUANTILES_monthly_Globe_1880_2023.mat');

BIC_FILE = fullfile(COMP_RES_DIR,'BIC.dat');

assert(isfile(DATA_FILE), ...
    'Input file QUANTILES_monthly_Globe_1880_2023.mat not found.');

assert(isfile(BIC_FILE), ...
    'BIC.dat not found.');

load(DATA_FILE,'QUANTILES_monthly');

QUANTILES = QUANTILES_monthly.Globe;
BIC_all = load(BIC_FILE);

% ------------------------------------------------------------------------
% 3. Settings
% ------------------------------------------------------------------------

Q = 1;              % global mean
q = Q;

W = 100;
select_w = 3;      % w=100 is the third rolling window
f_select = 25;     % Pareto selection horizon
select_f = 3;      % h=25 is the third horizon in [1,10,25,50]

f100 = 77;         % forecast horizon from 2023 to 2100

mod = 14;
nmodels_table = 17;

FORECAST_table = NaN(nmodels_table,1);
CI_fit_table   = NaN(nmodels_table,1,2);
CI_for_table   = NaN(nmodels_table,1,2);

SELECTED = NaN(1,14);

% ------------------------------------------------------------------------
% 4. Load Pareto selection for h=25, w=100
% ------------------------------------------------------------------------

f2 = strcat('for',num2str(f_select));
w2 = strcat('w',num2str(W));

GW_FILE = fullfile(COMP_RES_DIR, ...
    strcat('GW_rdos_',f2,'_',w2,'_all_quantiles.mat'));

assert(isfile(GW_FILE), ...
    ['Missing GW results file: ', GW_FILE]);

load(GW_FILE,'RDOS');

Z = QUANTILES_monthly.Globe(:,q);
t = length(Z);

selected_models = RDOS.GW.models(q,1:mod);
SELECTED(:,:) = selected_models;

Y = RDOS.Xf(:,q,:,select_w,select_f);
Y = squeeze(Y(:,:,1:mod));
Y = Y(:,selected_models==1);

m = size(Y,2);

% ------------------------------------------------------------------------
% 5. Combination weights using selected models
% ------------------------------------------------------------------------

% BIC weights
BIC = BIC_all(q,selected_models==1);

Weig_sbic = NaN(1,m);

for k = 1:m
    Weig_sbic(k) = exp(-0.5*BIC(k)) / sum(exp(-0.5*BIC));
end

% Simple mean weights
Weig_mean = (1/m)*ones(1,m);

% Estimated weights constrained to sum to one
X = Y;

if size(X,2) == 1
    BETA = 1;
else
    beta = compute_betas_sum1_all(Z(W+f_select:t,1),X(W+f_select:t,:));
    BETA = beta(2:end);
end

Weig_betas = BETA';

clear global;

% ------------------------------------------------------------------------
% 6. Forecast to 2100 and confidence intervals
% ------------------------------------------------------------------------

[FORECAST,CI_fit,RES] = compute_forecast_models(Z,f100);

FORECAST2 = FORECAST(selected_models==1);
RES = RES(:,selected_models==1);

% Combination 1: BIC weights
FOR_comb1 = Weig_sbic*FORECAST2;

[ci_for_sbic,ci_models] = ...
    compute_CI_comb_select_models(QUANTILES,RDOS,q,W,f_select, ...
    selected_models,Weig_sbic,0);

CI_for_comb1 = FOR_comb1 + ci_for_sbic;
CI_for = FORECAST2 + ci_models;

var_comb1 = compute_cov_comb(RES,Weig_sbic);
CI_fit_comb1 = FOR_comb1 + ...
    [norminv(0.05)*sqrt(var_comb1), ...
     norminv(0.95)*sqrt(var_comb1)];

% Combination 2: simple mean
FOR_comb2 = mean(FORECAST2);

ci_for_mean = ...
    compute_CI_comb_select_models(QUANTILES,RDOS,q,W,f_select, ...
    selected_models,Weig_mean,0);

CI_for_comb2 = FOR_comb2 + ci_for_mean;

var_comb2 = compute_cov_comb(RES,Weig_mean);
CI_fit_comb2 = FOR_comb2 + ...
    [norminv(0.05)*sqrt(var_comb2), ...
     norminv(0.95)*sqrt(var_comb2)];

% Combination 3: estimated weights
FOR_comb3 = FORECAST2'*BETA;

ci_for_betas = ...
    compute_CI_comb_select_models(QUANTILES,RDOS,q,W,f_select, ...
    selected_models,Weig_betas,0);

CI_for_comb3 = FOR_comb3 + ci_for_betas;

var_comb3 = compute_cov_comb(RES,Weig_betas);
CI_fit_comb3 = FOR_comb3 + ...
    [norminv(0.05)*sqrt(var_comb3), ...
     norminv(0.95)*sqrt(var_comb3)];

% ------------------------------------------------------------------------
% 7. Store results
% ------------------------------------------------------------------------

FORECAST_table(:,1) = [FORECAST;FOR_comb1;FOR_comb2;FOR_comb3];

CI_fit_table(:,1,1) = ...
    [CI_fit(:,1);CI_fit_comb1(1);CI_fit_comb2(1);CI_fit_comb3(1)];

CI_fit_table(:,1,2) = ...
    [CI_fit(:,2);CI_fit_comb1(2);CI_fit_comb2(2);CI_fit_comb3(2)];

pos = find(selected_models==1);

CI_for_table(pos',1,1) = CI_for(:,1);
CI_for_table(mod+1:nmodels_table,1,1) = ...
    [CI_for_comb1(1);CI_for_comb2(1);CI_for_comb3(1)];

CI_for_table(pos',1,2) = CI_for(:,2);
CI_for_table(mod+1:nmodels_table,1,2) = ...
    [CI_for_comb1(2);CI_for_comb2(2);CI_for_comb3(2)];

F_comb = [FOR_comb1;FOR_comb2;FOR_comb3];
CI_for_comb = [CI_for_comb1;CI_for_comb2;CI_for_comb3];
CI_fit_comb = [CI_fit_comb1;CI_fit_comb2;CI_fit_comb3];

disp('The 2100 forecast values of combined models are:')
F_comb

disp('Their CI with forecast errors are:')
CI_for_comb

disp('Their CI with fitted errors are:')
CI_fit_comb

% ------------------------------------------------------------------------
% 8. LaTeX tables
% ------------------------------------------------------------------------

models = {'mean','linear-trend','pol-trend','pol-trend-av-sl', ...
    'pol-trend-log','struct-breaks','pol-trend-arp', ...
    'pol-trend-arp-av-sl','arp','rw','rwd','ima','arfima','arp20', ...
    'combined 1','combined 2','combined 3'};

% Forecast-error intervals
TAB_FILE_FOR = fullfile(TAB_DIR, ...
    'Table_forecast_Pareto_mean_ci_for_w100_2100_Globe_1880_2023.tex');

fid = fopen(TAB_FILE_FOR,'w');

fprintf(fid,'\\begin{tabular}{lc}\n');
fprintf(fid,'\\toprule\n');
fprintf(fid,'Models & 2100 forecast \\\\\n');
fprintf(fid,'\\midrule\n');

for i = 1:length(models)
    
    fprintf(fid,'%s & %5.2f \\\\\n',models{i},FORECAST_table(i,1));
    
    fprintf(fid,' & (%5.2f,%5.2f) \\\\\n', ...
        CI_for_table(i,1,1),CI_for_table(i,1,2));
    
end

fprintf(fid,'\\bottomrule\n');
fprintf(fid,'\\end{tabular}\n');

fclose(fid);

% Fitted-error intervals
TAB_FILE_FIT = fullfile(TAB_DIR, ...
    'Table_forecast_Pareto_mean_ci_fit_w100_2100_Globe_1880_2023.tex');

fid = fopen(TAB_FILE_FIT,'w');

fprintf(fid,'\\begin{tabular}{lc}\n');
fprintf(fid,'\\toprule\n');
fprintf(fid,'Models & 2100 forecast \\\\\n');
fprintf(fid,'\\midrule\n');

for i = 1:length(models)
    
    fprintf(fid,'%s & %5.2f \\\\\n',models{i},FORECAST_table(i,1));
    
    fprintf(fid,' & (%5.2f,%5.2f) \\\\\n', ...
        CI_fit_table(i,1,1),CI_fit_table(i,1,2));
    
end

fprintf(fid,'\\bottomrule\n');
fprintf(fid,'\\end{tabular}\n');

fclose(fid);

% ------------------------------------------------------------------------
% 9. Save results
% ------------------------------------------------------------------------

RDOS_FILE = fullfile(RES_DIR,'RDOS_for_mean_Globe_1880_2023.mat');

if isfile(RDOS_FILE)
    Sload = load(RDOS_FILE);
    
    if isfield(Sload,'RDOS')
        RDOS_mean = Sload.RDOS;
    else
        RDOS_mean = struct();
    end
else
    RDOS_mean = struct();
end

RDOS_mean.comb_Pareto_2100.H_selection = f_select;
RDOS_mean.comb_Pareto_2100.H_forecast  = f100;
RDOS_mean.comb_Pareto_2100.W = W;
RDOS_mean.comb_Pareto_2100.models = models;
RDOS_mean.comb_Pareto_2100.forecast = FORECAST_table;
RDOS_mean.comb_Pareto_2100.ci.fit = CI_fit_table;
RDOS_mean.comb_Pareto_2100.ci.for = CI_for_table;
RDOS_mean.comb_Pareto_2100.selected = SELECTED;

RDOS = RDOS_mean;

save(RDOS_FILE,'RDOS');

fprintf('\nPareto-combined 2100 forecast for the global mean completed.\n');
fprintf('Tables saved in:\n%s\n', TAB_DIR);
fprintf('Results saved in:\n%s\n', RES_DIR);