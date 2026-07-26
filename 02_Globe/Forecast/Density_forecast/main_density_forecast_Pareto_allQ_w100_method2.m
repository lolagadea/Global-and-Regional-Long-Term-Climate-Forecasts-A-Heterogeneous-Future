% ========================================================================
% main_density_forecast_Pareto_allQ_w100_method2.m
% ------------------------------------------------------------------------
% Density forecasts for all selected quantiles using Pareto-superior models.
%
% Method 2:
% Models are selected jointly across quantiles/characteristics using the
% Pareto-superiority criterion. When no model is Pareto superior for all
% quantiles simultaneously, the models selected are those that are Pareto
% superior for the largest number of quantiles.
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
clear global;
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
ARFIMA_DIR        = fullfile(ROOT_DIR,'ARFIMA');

INTRO_RES_DIR = fullfile(GLOBE_DIR,'Introduction','Results');
COMP_RES_DIR  = fullfile(FORECAST_DIR,'run_forecast_competition','Results_for');

FIG_DIR = fullfile(THIS_DIR,'Figures');
TAB_DIR = fullfile(THIS_DIR,'Tables');
RES_DIR = fullfile(THIS_DIR,'Results');

if ~exist(FIG_DIR,'dir'), mkdir(FIG_DIR); end
if ~exist(TAB_DIR,'dir'), mkdir(TAB_DIR); end
if ~exist(RES_DIR,'dir'), mkdir(RES_DIR); end

if exist(FUN_GENERAL_DIR,'dir'),  addpath(genpath(FUN_GENERAL_DIR));  end
if exist(FUN_FORECAST_DIR,'dir'), addpath(genpath(FUN_FORECAST_DIR)); end
if exist(ARFIMA_DIR,'dir'),       addpath(genpath(ARFIMA_DIR));       end
if exist(COMP_RES_DIR,'dir'),     addpath(genpath(COMP_RES_DIR));     end

% ------------------------------------------------------------------------
% 2. Load inputs
% ------------------------------------------------------------------------

FILE_Q      = fullfile(INTRO_RES_DIR,'QUANTILES_monthly_Globe_1880_2023.mat');
FILE_RDOS   = fullfile(COMP_RES_DIR,'RDOS_Globe_1880_2023.mat');
FILE_SEL1   = fullfile(COMP_RES_DIR,'selected_the_model_allQ.mat');
FILE_SEL2   = fullfile(COMP_RES_DIR,'selected_the_model2_allQ.mat');
FILE_ONE    = fullfile(COMP_RES_DIR,'RDOS_one_model.mat');
FILE_BIC    = fullfile(COMP_RES_DIR,'BIC.dat');

assert(exist(FILE_Q,'file')==2,    'Input not found: %s', FILE_Q);
assert(exist(FILE_RDOS,'file')==2, 'Input not found: %s', FILE_RDOS);
assert(exist(FILE_SEL1,'file')==2, 'Input not found: %s', FILE_SEL1);
assert(exist(FILE_SEL2,'file')==2, 'Input not found: %s', FILE_SEL2);
assert(exist(FILE_ONE,'file')==2,  'Input not found: %s', FILE_ONE);
assert(exist(FILE_BIC,'file')==2,  'Input not found: %s', FILE_BIC);

load(FILE_Q,'QUANTILES_monthly');
QUANTILES = QUANTILES_monthly.Globe;
[t,~] = size(QUANTILES);

load(FILE_RDOS,'RDOS');
RDOS_run_comp = RDOS; 

load(FILE_SEL1);  % expected variables: MODELS, possibly others
load(FILE_SEL2);  % expected variables: MODELS2, possibly others
load(FILE_ONE,'RDOS_one_model');

BIC_all = load(FILE_BIC);

models = {'mean','linear-trend','pol-trend','pol-trend-av-sl', ...
          'pol-trend-log','struct-breaks','pol-trend-arp', ...
          'pol-trend-arp-av-sl','arp','rw','rwd','ima','arfima','arp20'};

% ------------------------------------------------------------------------
% 3. Settings
% ------------------------------------------------------------------------

F = [1 10 25];
W = 100;
select_w = 3;      % position of w=100 in the original competition arrays
mod = 14;          % number of individual models, excluding combinations
Q = 9:19;          % q05, q10, ..., q95

RDOS_FOR = NaN(length(Q),3,length(W),length(F));
RDOS_FIT = NaN(length(Q),3,length(W),length(F));

% ------------------------------------------------------------------------
% 4. Forecast loop
% ------------------------------------------------------------------------

for j = 1:length(W)

    w = W(j);

    for i = 1:length(F)

        f = F(i);

        selected_models_inters = RDOS_one_model.select(:,select_w,i);      % method 2
        selected_models_loss   = (MODELS(:,select_w,i))';                 
        selected_models_loss2  = (MODELS2(:,select_w,i))';                %#ok<NASGU>
        selected_models_first  = RDOS_one_model.firstbest(:,select_w,i);  %#ok<NASGU>
        selected_models_second = RDOS_one_model.secondbest(:,select_w,i); %#ok<NASGU>

        for q = 1:length(Q)

            fprintf('Method 2 | w=%d | h=%d | quantile index=%d\n', w, f, Q(q));

            Z = QUANTILES(:,Q(q));

            Y = RDOS.Xf(:,Q(q),:,select_w,i);
            Y = squeeze(Y(:,:,1:mod));              % remove combined models
            Y = Y(:,selected_models_inters==1);

            m = size(Y,2);

            % ------------------------------------------------------------
            % Combine selected models using BIC weights
            % ------------------------------------------------------------

            BIC = BIC_all(Q(q),selected_models_inters==1);

            Weig_sbic = NaN(1,m);
            for k = 1:m
                Weig_sbic(k) = exp(-0.5*BIC(k)) / sum(exp(-0.5*BIC(1,:)));
            end

            C1 = NaN(t,1); %#ok<NASGU>
            C1(1:w,:) = Z(1:w,:);
            C1(w+f:t) = Y(w+f:t,:) * Weig_sbic';

            % ------------------------------------------------------------
            % Forecast and confidence intervals
            % ------------------------------------------------------------

            [FORECAST,CI_fit,RES] = compute_forecast_models(Z,f); %#ok<ASGLU>

            FORECAST2 = FORECAST(selected_models_inters==1);

            RES = RES(:,selected_models_inters==1);
            RES = rmmissing(RES);

            FOR_comb = Weig_sbic * FORECAST2;

            var_comb = compute_cov_comb(RES,Weig_sbic);
            CI_fit_comb = FOR_comb + ...
                [norminv(0.05)*sqrt(var_comb), norminv(0.95)*sqrt(var_comb)];

            if sum(selected_models_inters==1)==0

                warning('There are no selected models for q=%d, h=%d.', Q(q), f);
                CI_for_comb = [NaN NaN];

            elseif sum(selected_models_inters==1)==1

                varres = RES' * RES / t;
                CI_for_comb = FOR_comb + norminv([0.05 0.95]) * sqrt(varres);

            else

                [ci_for_sbic, ci_models] = compute_CI_comb_select_models( ...
                    QUANTILES,RDOS,Q(q),w,f,selected_models_inters,Weig_sbic,0); %#ok<ASGLU>

                CI_for_comb = FOR_comb + ci_for_sbic;

            end

            RDOS_FOR(q,1,j,i)   = FOR_comb;
            RDOS_FOR(q,2:3,j,i) = CI_for_comb;

            RDOS_FIT(q,1,j,i)   = FOR_comb;
            RDOS_FIT(q,2:3,j,i) = CI_fit_comb;

        end
    end
end

% ------------------------------------------------------------------------
% 5. Tables
% ------------------------------------------------------------------------

rows = {'q05','q10','q20','q30','q40','q50','q60','q70','q80','q90','q95'};

caption_for = ['Long-term density forecast with Pareto-superior models ' ...
    'for all quantiles (the Globe, CRU data, 1880--2023, w=100, method 2)'];
label_for = 'tab-Pareto-one-model-all-quantiles-method2-Globe-1880-2023';

write_density_table( ...
    fullfile(TAB_DIR,'Table_onemodel_allquantiles_w100_method2_for_Globe_1880_2023.tex'), ...
    RDOS_FOR, F, rows, caption_for, label_for);

caption_fit = ['Long-term density forecast with Pareto-superior models ' ...
    'for all quantiles, fit uncertainty (the Globe, CRU data, 1880--2023, w=100, method 2)'];
label_fit = 'tab-Pareto-one-model-all-quantiles-method2-fit-Globe-1880-2023';

write_density_table( ...
    fullfile(TAB_DIR,'Table_onemodel_allquantiles_w100_method2_fit_Globe_1880_2023.tex'), ...
    RDOS_FIT, F, rows, caption_fit, label_fit);

% ------------------------------------------------------------------------
% 6. Display selected models
% ------------------------------------------------------------------------

for i = 1:length(F)

    f = F(i);
    selected_models_inters = RDOS_one_model.select(:,select_w,i);

    fprintf('\nThe selected models with method 2 for h=%d are:\n', f);
    disp(models(selected_models_inters==1)');

end

% ------------------------------------------------------------------------
% 7. Save results
% ------------------------------------------------------------------------

save(fullfile(RES_DIR,'RDOS_FOR_allQ_method2_w100_Globe_1880_2023.mat'), ...
    'RDOS_FOR');

save(fullfile(RES_DIR,'RDOS_FIT_allQ_method2_w100_Globe_1880_2023.mat'), ...
    'RDOS_FIT');

fprintf('\nMethod 2 completed successfully.\n');
fprintf('Tables saved in:  %s\n', TAB_DIR);
fprintf('Results saved in: %s\n', RES_DIR);

% ========================================================================
% Local helper
% ========================================================================

function write_density_table(outfile, RDOS_OBJ, F, rows, tableCaption, tableLabel)
%WRITE_DENSITY_TABLE Write LaTeX table for density forecasts.

fid = fopen(outfile,'w');
assert(fid~=-1, 'Could not open file for writing: %s', outfile);

fprintf(fid, ['\\begin{table}[h!]\\caption{%s}\\label{%s}' ...
    '\\begin{center}\\scalebox{0.8}{\\begin{tabular}{l|ccc}' ...
    '\\hline \\hline \n'], tableCaption, tableLabel);

fprintf(fid, ' & h=1 & h=10 & h=25 \\\\ \\hline \n');

for r = 1:length(rows)

    point_line = '';
    ci_line    = '';

    for c = 1:length(F)

        point_line = [point_line num2str(RDOS_OBJ(r,1,1,c),'%5.2f') '&']; %#ok<AGROW>

        ci_line = [ci_line ' (' ...
            num2str(RDOS_OBJ(r,2,1,c),'%5.2f') ',' ...
            num2str(RDOS_OBJ(r,3,1,c),'%5.2f') ')' '&']; %#ok<AGROW>

    end

    fprintf(fid, '%s & %s \\\\ \n', rows{r}, point_line(1:end-1));
    fprintf(fid, ' & %s \\\\ \n', ci_line(1:end-1));

end

fprintf(fid, '\\hline \\hline \\end{tabular}}\\end{center}\\end{table}\n');
fclose(fid);

end
