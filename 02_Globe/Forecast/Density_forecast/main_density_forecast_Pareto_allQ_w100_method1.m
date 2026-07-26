% ========================================================================
% main_density_forecast_Pareto_allQ_w100_method1.m
% ------------------------------------------------------------------------
% Density forecasts for all selected quantiles using Pareto-superior models.
%
% Method 1:
% Models are selected using the Pareto-superiority criterion applied jointly
% across quantiles. The loss functions are aggregated over quantiles and
% the Giacomini-White test is used to identify the Pareto-superior set.
%
% Selected models are combined using BIC weights.
%
% Globe, CRU data, 1880--2023
% Window: w = 100
% Horizons: h = 1, 10, 25
%
% This version is portable and writes outputs to:
%   Density_forecast/Figures
%   Density_forecast/Tables
%   Density_forecast/Results
% ========================================================================

clear; clc;
warning('off');

% ------------------------------------------------------------------------
% 1. Portable paths
% ------------------------------------------------------------------------

THIS_DIR     = fileparts(mfilename('fullpath'));
FORECAST_DIR = fileparts(THIS_DIR);
GLOBE_DIR    = fileparts(FORECAST_DIR);
ROOT_DIR     = fileparts(GLOBE_DIR);

FUN_GENERAL_DIR  = fullfile(ROOT_DIR,'functions');
FUN_FORECAST_DIR = fullfile(ROOT_DIR,'functions_for');

INTRO_RES_DIR = fullfile(GLOBE_DIR,'Introduction','Results');
COMP_RES_DIR  = fullfile(FORECAST_DIR,'run_forecast_competition','Results_for');

FIG_DIR = fullfile(THIS_DIR,'Figures');
TAB_DIR = fullfile(THIS_DIR,'Tables');
RES_DIR = fullfile(THIS_DIR,'Results');

if ~exist(FIG_DIR,'dir'), mkdir(FIG_DIR); end
if ~exist(TAB_DIR,'dir'), mkdir(TAB_DIR); end
if ~exist(RES_DIR,'dir'), mkdir(RES_DIR); end

addpath(genpath(FUN_GENERAL_DIR));
addpath(genpath(FUN_FORECAST_DIR));
addpath(genpath(COMP_RES_DIR));

% ------------------------------------------------------------------------
% 2. Load inputs
% ------------------------------------------------------------------------

FILE_QUANT = fullfile(INTRO_RES_DIR,'QUANTILES_monthly_Globe_1880_2023.mat');
FILE_RDOS  = fullfile(COMP_RES_DIR,'RDOS_Globe_1880_2023.mat');
FILE_SEL1  = fullfile(COMP_RES_DIR,'selected_the_model_allQ.mat');
FILE_SEL2  = fullfile(COMP_RES_DIR,'selected_the_model2_allQ.mat');
FILE_ONE   = fullfile(COMP_RES_DIR,'RDOS_one_model.mat');
FILE_BIC   = fullfile(COMP_RES_DIR,'BIC.dat');

assert(exist(FILE_QUANT,'file')==2, 'Missing input file: %s', FILE_QUANT);
assert(exist(FILE_RDOS,'file')==2,  'Missing input file: %s', FILE_RDOS);
assert(exist(FILE_SEL1,'file')==2,  'Missing input file: %s', FILE_SEL1);
assert(exist(FILE_SEL2,'file')==2,  'Missing input file: %s', FILE_SEL2);
assert(exist(FILE_ONE,'file')==2,   'Missing input file: %s', FILE_ONE);
assert(exist(FILE_BIC,'file')==2,   'Missing input file: %s', FILE_BIC);

load(FILE_QUANT,'QUANTILES_monthly');
load(FILE_RDOS,'RDOS');
load(FILE_SEL1);   % expected variables: MODELS
load(FILE_SEL2);   % expected variables: MODELS2
load(FILE_ONE,'RDOS_one_model');

QUANTILES = QUANTILES_monthly.Globe;
RDOS_run_comp = RDOS; %#ok<NASGU>

[t,~] = size(QUANTILES);

% ------------------------------------------------------------------------
% 3. Settings
% ------------------------------------------------------------------------

models = {'mean','linear-trend','pol-trend','pol-trend-av-sl', ...
          'pol-trend-log','struct-breaks','pol-trend-arp', ...
          'pol-trend-arp-av-sl','arp','rw','rwd','ima','arfima','arp20'};

F = [1,10,25];
W = 100;

select_w = 3;      % w = 100 in the original forecast competition arrays
mod = 14;          % number of individual models before combinations
Q = 9:19;          % q05, q10, ..., q95

RDOS_FOR = NaN(length(Q),3,length(W),length(F));
RDOS_FIT = NaN(length(Q),3,length(W),length(F));

% ------------------------------------------------------------------------
% 4. Density forecasts: Method 1
% ------------------------------------------------------------------------

for j = 1:length(W)

    w = W(j);

    for i = 1:length(F)

        f = F(i);

        selected_models_loss  = (MODELS(:,select_w,i))';
        selected_models_loss2 = (MODELS2(:,select_w,i))';

        % Kept for traceability with the old code, although method 1 uses
        % selected_models_loss and BIC weights.
        selected_models_inters = RDOS_one_model.select(:,select_w,i); %#ok<NASGU>
        selected_models_first  = RDOS_one_model.firstbest(:,select_w,i); %#ok<NASGU>
        selected_models_second = RDOS_one_model.secondbest(:,select_w,i); %#ok<NASGU>

        if ~any(selected_models_loss == 1)
            warning('No selected models in MODELS. Using MODELS2 instead.');
            selected_models_loss = selected_models_loss2;
        end

        for q = 1:length(Q)

            Z = QUANTILES(:,Q(q));

            Y = RDOS.Xf(:,Q(q),:,select_w,i);
            Y = squeeze(Y(:,:,1:mod));   % remove combined models

            Y = Y(:,selected_models_loss == 1);
            m = size(Y,2);

            % ------------------------------------------------------------
            % BIC weights for selected Pareto-superior models
            % ------------------------------------------------------------

            BIC_all = load(FILE_BIC);
            BIC_q = BIC_all(Q(q),selected_models_loss == 1);

            Weig_sbic = NaN(1,m);
            denom = sum(exp(-0.5 * BIC_q));

            for k = 1:m
                Weig_sbic(k) = exp(-0.5 * BIC_q(k)) / denom;
            end

            % ------------------------------------------------------------
            % Forecast and confidence intervals
            % ------------------------------------------------------------

            C1 = NaN(t,1); %#ok<NASGU>
            C1(1:w,:) = Z(1:w,:);
            C1(w+f:t) = Y(w+f:t,:) * Weig_sbic';

            [FORECAST,CI_fit,RES] = compute_forecast_models(Z,f); %#ok<ASGLU>

            FORECAST2 = FORECAST(selected_models_loss == 1);
            RES = RES(:,selected_models_loss == 1);
            RES = rmmissing(RES);

            FOR_comb = Weig_sbic * FORECAST2;

            var_comb = compute_cov_comb(RES,Weig_sbic);
            CI_fit_comb = FOR_comb + ...
                [norminv(0.05)*sqrt(var_comb), norminv(0.95)*sqrt(var_comb)];

            n_selected = sum(selected_models_loss == 1);

            if n_selected == 1
                varres = RES' * RES / t;
                CI_for_comb = FOR_comb + norminv([0.05 0.95]) * sqrt(varres);
            elseif n_selected > 1
                [ci_for_sbic, ~] = compute_CI_comb_select_models( ...
                    QUANTILES,RDOS,Q(q),w,f,selected_models_loss,Weig_sbic,0);
                CI_for_comb = FOR_comb + ci_for_sbic;
            else
                error('No selected models available for q=%d, f=%d, w=%d.', Q(q), f, w);
            end

            RDOS_FOR(q,1,j,i)   = FOR_comb;
            RDOS_FOR(q,2:3,j,i) = CI_for_comb;

            RDOS_FIT(q,1,j,i)   = FOR_comb;
            RDOS_FIT(q,2:3,j,i) = CI_fit_comb;

        end
    end
end

% ------------------------------------------------------------------------
% 5. Export LaTeX tables
% ------------------------------------------------------------------------

rows = {'q05','q10','q20','q30','q40','q50','q60','q70','q80','q90','q95'};

tableCaptionFor = ['Long-term density forecast with Pareto-superior ' ...
    'models for all quantiles (the Globe, CRU data, 1880--2023, w=100, method 1)'];
tableLabelFor = 'tab-Pareto-all-quantiles-method1-Globe-1880-2023';

fileFor = fullfile(TAB_DIR,'Table_onemodel_allquantiles_w100_method1_for_Globe_1880_2023.tex');
write_density_table_method1(fileFor, tableCaptionFor, tableLabelFor, rows, F, RDOS_FOR);

tableCaptionFit = ['Long-term density forecast with Pareto-superior ' ...
    'models for all quantiles, fitted intervals (the Globe, CRU data, 1880--2023, w=100, method 1)'];
tableLabelFit = 'tab-Pareto-all-quantiles-method1-fit-Globe-1880-2023';

fileFit = fullfile(TAB_DIR,'Table_onemodel_allquantiles_w100_method1_fit_Globe_1880_2023.tex');
write_density_table_method1(fileFit, tableCaptionFit, tableLabelFit, rows, F, RDOS_FIT);

% ------------------------------------------------------------------------
% 6. Display selected models
% ------------------------------------------------------------------------

for i = 1:length(F)
    selected_models_loss = (MODELS(:,select_w,i))';

    fprintf('\nThe selected models with method 1 for h=%d are:\n', F(i));
    disp(models(selected_models_loss == 1)');
end

% ------------------------------------------------------------------------
% 7. Save results
% ------------------------------------------------------------------------

save(fullfile(RES_DIR,'RDOS_FOR_allQ_method1_w100_Globe_1880_2023.mat'), ...
    'RDOS_FOR','F','W','Q','rows');

save(fullfile(RES_DIR,'RDOS_FIT_allQ_method1_w100_Globe_1880_2023.mat'), ...
    'RDOS_FIT','F','W','Q','rows');

fprintf('\nMethod 1 density forecast completed successfully.\n');
fprintf('Tables saved in:  %s\n', TAB_DIR);
fprintf('Results saved in: %s\n', RES_DIR);

% ========================================================================
% Local helper
% ========================================================================

function write_density_table_method1(file, tableCaption, tableLabel, rows, F, RDOS_OBJ)

fid = fopen(file,'w');
assert(fid > 0, 'Cannot open output file: %s', file);

fprintf(fid, ['\\begin{table}[h!]\\caption{%s}\\label{%s}' ...
    '\\begin{center}\\scalebox{0.8}{\\begin{tabular}{l|ccc}' ...
    '\\hline \\hline \n'], tableCaption, tableLabel);

fprintf(fid, ' & h=1 & h=10 & h=25 \\\\ \\hline \n');

for i = 1:length(rows)

    lineForecast = '';
    lineCI = '';

    for j = 1:length(F)
        lineForecast = [lineForecast, sprintf('%5.2f&', RDOS_OBJ(i,1,1,j))]; %#ok<AGROW>
        lineCI = [lineCI, sprintf('(%5.2f,%5.2f)&', ...
            RDOS_OBJ(i,2,1,j), RDOS_OBJ(i,3,1,j))]; %#ok<AGROW>
    end

    fprintf(fid, '%s&%s \\\\ \n', rows{i}, lineForecast(1:end-1));
    fprintf(fid, '&%s \\\\ \n', lineCI(1:end-1));

end

fprintf(fid, '\\hline \\hline \\end{tabular}}\\end{center}\\end{table}\n');
fclose(fid);

end
