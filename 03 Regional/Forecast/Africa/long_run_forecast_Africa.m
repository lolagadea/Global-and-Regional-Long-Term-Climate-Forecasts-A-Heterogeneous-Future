% ========================================================================
% long_run_forecast_Africa.m
% ------------------------------------------------------------------------
% Long-run density forecast for 2100 using the three Pareto methods.
% Region: Africa.
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

RES_DIR = fullfile(THIS_DIR,'Results_for');
TAB_DIR = fullfile(THIS_DIR,'Tables');

if ~isfolder(RES_DIR), mkdir(RES_DIR); end
if ~isfolder(TAB_DIR), mkdir(TAB_DIR); end

% ------------------------------------------------------------------------
% 2. Load inputs
% ------------------------------------------------------------------------

DATA_FILE = fullfile(INTRO_RES_DIR,'QUANTILES_monthly_1960_2023.mat');
RDOS_FILE = fullfile(RES_DIR,'RDOS_Africa_1960_2023.mat');
BIC_FILE  = fullfile(RES_DIR,'BIC_Africa.dat');

MODEL_FILE1 = fullfile(RES_DIR,'selected_the model_allQ.mat');
MODEL_FILE2 = fullfile(RES_DIR,'selected_the model2_allQ.mat');
ONE_FILE    = fullfile(RES_DIR,'RDOS_one_model.mat');

assert(isfile(DATA_FILE),   'QUANTILES_monthly_1960_2023.mat not found.');
assert(isfile(RDOS_FILE),   'RDOS_Africa_1960_2023.mat not found.');
assert(isfile(BIC_FILE),    'BIC_Africa.dat not found.');
assert(isfile(MODEL_FILE1), 'selected_the model_allQ.mat not found.');
assert(isfile(MODEL_FILE2), 'selected_the model2_allQ.mat not found.');
assert(isfile(ONE_FILE),    'RDOS_one_model.mat not found.');

load(DATA_FILE,'QUANTILES_monthly');
load(RDOS_FILE,'RDOS');

BIC = load(BIC_FILE);

QUANTILES = QUANTILES_monthly.Africa;

% ------------------------------------------------------------------------
% 3. Settings
% ------------------------------------------------------------------------

F = [1,10,25];
W = 25;

select_w = 1;
select_f = 3;      % h = 25
f = F(select_f);
w = W(select_w);

f2 = strcat('for',num2str(f));
w2 = strcat('w',num2str(w));

Q = 9:19;          % q05-q95
mod = 14;

f100 = 77;         % forecast horizon corresponding to year 2100

GW_FILE = fullfile(RES_DIR, ...
    strcat('GW_rdos_',f2,'_',w2,'_all_quantiles.mat'));

assert(isfile(GW_FILE), 'GW file not found: %s', GW_FILE);

load(GW_FILE,'RDOS');

% ------------------------------------------------------------------------
% 4. Initialize output
% ------------------------------------------------------------------------

FOR_comb.m0.for    = NaN(length(Q),1);
FOR_comb.m0.ci.fit = NaN(length(Q),2);
FOR_comb.m0.ci.for = NaN(length(Q),2);

FOR_comb.m1.for    = NaN(length(Q),1);
FOR_comb.m1.ci.fit = NaN(length(Q),2);
FOR_comb.m1.ci.for = NaN(length(Q),2);

FOR_comb.m2.for    = NaN(length(Q),1);
FOR_comb.m2.ci.fit = NaN(length(Q),2);
FOR_comb.m2.ci.for = NaN(length(Q),2);

% ========================================================================
% Method 0
% ========================================================================

for q = 1:length(Q)

    selected_models = RDOS.GW.models(Q(q),1:mod);

    Z = QUANTILES(:,Q(q));

    Y = RDOS.Xf(:,Q(q),:,select_w,select_f);
    Y = squeeze(Y(:,:,1:mod));
    Y = Y(:,selected_models==1);

    m = size(Y,2);

    BIC2 = BIC(Q(q),selected_models==1);

    Weig_sbic = NaN(1,m);

    for k = 1:m
        Weig_sbic(k) = exp(-1/2*BIC2(k)) / ...
            sum(exp(-1/2*BIC2(1,:)));
    end

    [FORECAST,ci_fit,RES] = compute_forecast_models(Z,f100); %#ok<ASGLU>

    FORECAST2 = FORECAST(selected_models==1);
    RES = RES(:,selected_models==1);

    FOR_comb.m0.for(q) = Weig_sbic*FORECAST2;

    [ci_for_sbic,ci_models] = compute_CI_comb_select_models( ...
        QUANTILES,RDOS,Q(q),w,f,selected_models,Weig_sbic,1); %#ok<ASGLU>

    FOR_comb.m0.ci.for(q,:) = FOR_comb.m0.for(q) + ci_for_sbic;

    var_comb = compute_cov_comb(RES,Weig_sbic);

    FOR_comb.m0.ci.fit(q,:) = FOR_comb.m0.for(q) + ...
        [norminv(0.05)*sqrt(var_comb), norminv(0.95)*sqrt(var_comb)];

end

% ========================================================================
% Method 1
% ========================================================================

load(MODEL_FILE1,'MODELS');
load(MODEL_FILE2,'MODELS2');

selected_models_loss  = (MODELS(:,select_w,select_f))';
selected_models_loss2 = (MODELS2(:,select_w,select_f))';

for q = 1:length(Q)

    Z = QUANTILES(:,Q(q));

    Y = RDOS.Xf(:,Q(q),:,select_w,select_f);
    Y = squeeze(Y(:,:,1:mod));

    if length(find(selected_models_loss==1)) == 0
        disp('Warning: there are no selected models; using selected_models_loss2.')
        selected_models = selected_models_loss2;
    else
        selected_models = selected_models_loss;
    end

    Y = Y(:,selected_models==1);

    m = size(Y,2);

    BIC2 = BIC(Q(q),selected_models==1);

    Weig_sbic = NaN(1,m);

    for k = 1:m
        Weig_sbic(k) = exp(-1/2*BIC2(k)) / ...
            sum(exp(-1/2*BIC2(1,:)));
    end

    [FORECAST,ci_fit,RES] = compute_forecast_models(Z,f100); %#ok<ASGLU>

    FORECAST2 = FORECAST(selected_models==1);
    RES = RES(:,selected_models==1);

    FOR_comb.m1.for(q) = Weig_sbic*FORECAST2;

    [ci_for_sbic,ci_models] = compute_CI_comb_select_models( ...
        QUANTILES,RDOS,Q(q),w,f,selected_models,Weig_sbic,1); %#ok<ASGLU>

    FOR_comb.m1.ci.for(q,:) = FOR_comb.m1.for(q) + ci_for_sbic;

    var_comb = compute_cov_comb(RES,Weig_sbic);

    FOR_comb.m1.ci.fit(q,:) = FOR_comb.m1.for(q) + ...
        [norminv(0.05)*sqrt(var_comb), norminv(0.95)*sqrt(var_comb)];

end

% ========================================================================
% Method 2
% ========================================================================

load(ONE_FILE,'RDOS_one_model');

for q = 1:length(Q)

    Z = QUANTILES(:,Q(q));

    Y = RDOS.Xf(:,Q(q),:,select_w,select_f);
    Y = squeeze(Y(:,:,1:mod));

    selected_models = RDOS_one_model.select(:,select_w,select_f);

    Y = Y(:,selected_models==1);

    m = size(Y,2);

    BIC2 = BIC(Q(q),selected_models==1);

    Weig_sbic = NaN(1,m);

    for k = 1:m
        Weig_sbic(k) = exp(-1/2*BIC2(k)) / ...
            sum(exp(-1/2*BIC2(1,:)));
    end

    [FORECAST,ci_fit,RES] = compute_forecast_models(Z,f100); %#ok<ASGLU>

    FORECAST2 = FORECAST(selected_models==1);
    RES = RES(:,selected_models==1);

    FOR_comb.m2.for(q) = Weig_sbic*FORECAST2;

    [ci_for_sbic,ci_models] = compute_CI_comb_select_models( ...
        QUANTILES,RDOS,Q(q),w,f,selected_models,Weig_sbic,1); %#ok<ASGLU>

    FOR_comb.m2.ci.for(q,:) = FOR_comb.m2.for(q) + ci_for_sbic;

    var_comb = compute_cov_comb(RES,Weig_sbic);

    FOR_comb.m2.ci.fit(q,:) = FOR_comb.m2.for(q) + ...
        [norminv(0.05)*sqrt(var_comb), norminv(0.95)*sqrt(var_comb)];

end

% ------------------------------------------------------------------------
% 5. Save output
% ------------------------------------------------------------------------

OUT_FILE = fullfile(RES_DIR,'for_2100_Africa_1960_2023.mat');

save(OUT_FILE,'FOR_comb');

fprintf('\nLong-run forecast completed for Africa.\n');
fprintf('Saved in:\n%s\n',OUT_FILE);

clear global;