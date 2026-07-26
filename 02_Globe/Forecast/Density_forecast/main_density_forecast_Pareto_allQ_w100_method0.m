% ========================================================================
% main_density_forecast_Pareto_allQ_w100_method0.m
% ------------------------------------------------------------------------
% Density forecasts for all selected quantiles using Pareto-superior models.
%
% Method 0:
% Pareto-superior models are selected independently for each quantile.
%
% Globe, CRU / CRUTEM monthly quantiles, 1880--2023
% Window: w = 100
% Horizons: h = 1, 10, 25
%
% Outputs
%   - Tables  : LaTeX tables with forecast and fitted confidence intervals
%   - Results : MAT files with RDOS_FOR and RDOS_FIT
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
ARFIMA_DIR       = fullfile(ROOT_DIR,'ARFIMA');

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
addpath(genpath(ARFIMA_DIR));
addpath(genpath(COMP_RES_DIR));

% ------------------------------------------------------------------------
% 2. Load inputs
% ------------------------------------------------------------------------

load(fullfile(INTRO_RES_DIR,'QUANTILES_monthly_Globe_1880_2023.mat'), ...
    'QUANTILES_monthly');
QUANTILES = QUANTILES_monthly.Globe;
[t,~] = size(QUANTILES);

load(fullfile(COMP_RES_DIR,'RDOS_Globe_1880_2023.mat'),'RDOS');
RDOS_run_comp = RDOS; %#ok<NASGU>

load(fullfile(COMP_RES_DIR,'selected_the_model_allQ.mat'));
load(fullfile(COMP_RES_DIR,'selected_the_model2_allQ.mat'));
load(fullfile(COMP_RES_DIR,'RDOS_one_model.mat'));

% ------------------------------------------------------------------------
% 3. Settings
% ------------------------------------------------------------------------

models = {'mean','linear-trend','pol-trend','pol-trend-av-sl', ...
    'pol-trend-log','struct-breaks','pol-trend-arp', ...
    'pol-trend-arp-av-sl','arp','rw','rwd','ima','arfima','arp20'}; %#ok<NASGU>

F = [1,10,25];       % Forecast horizons
W = 100;             % Rolling window used in the paper
select_w = 3;        % Position of w=100 in RDOS.Xf from the competition code
mod = 14;            % Number of individual models before combinations
Q = 9:19;            % q05, q10, ..., q95

RDOS_FOR = NaN(length(Q),3,length(W),length(F));
RDOS_FIT = NaN(length(Q),3,length(W),length(F));

% ------------------------------------------------------------------------
% 4. Density forecast: method 0
% ------------------------------------------------------------------------

for j = 1:length(W)

    w = W(j);
    w2 = strcat('w',num2str(w));

    for i = 1:length(F)

        f = F(i);
        f2 = strcat('for',num2str(f));

        load(fullfile(COMP_RES_DIR, ...
            strcat('GW_rdos_',f2,'_',w2,'_all_quantiles.mat')),'RDOS');

        for q = 1:length(Q)

            selected_models = RDOS.GW.models(Q(q),1:mod);

            Z = QUANTILES(:,Q(q));
            Y = RDOS.Xf(:,Q(q),:,select_w,i);
            Y = squeeze(Y(:,:,1:mod)); % Remove combined models
            Y = Y(:,selected_models == 1);
            m = size(Y,2);

            % ------------------------------------------------------------
            % Combined model using BIC weights
            % ------------------------------------------------------------
            BIC_file = fullfile(COMP_RES_DIR,'BIC.dat');
            BIC = load(BIC_file);
            BIC = BIC(Q(q),selected_models == 1);

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
            [FORECAST,~,RES] = compute_forecast_models(Z,f);
            FORECAST2 = FORECAST(selected_models == 1);
            RES = RES(:,selected_models == 1);
            RES = rmmissing(RES);

            FOR_comb = Weig_sbic * FORECAST2;
            var_comb = compute_cov_comb(RES,Weig_sbic);
            CI_fit_comb = FOR_comb + ...
                [norminv(0.05)*sqrt(var_comb), norminv(0.95)*sqrt(var_comb)];

            n_selected = length(find(selected_models == 1));

            if n_selected == 0
                warning('There are no selected models for q=%d, h=%d, w=%d.', ...
                    Q(q), f, w);
                CI_for_comb = [NaN NaN];

            elseif n_selected == 1
                varres = RES' * RES / t;
                CI_for_comb = FOR_comb + norminv([0.05 0.95]) * sqrt(varres);

            else
                [ci_for_sbic, ~] = compute_CI_comb_select_models( ...
                    QUANTILES,RDOS,Q(q),w,f,selected_models,Weig_sbic,0);
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
% 5. Save tables
% ------------------------------------------------------------------------

rows = {'q05','q10','q20','q30','q40','q50','q60','q70','q80','q90','q95'};

caption_for = ['Long-term density forecast with Pareto-superior models ', ...
    'for all quantiles (Globe, CRU data, 1880--2023, w=100, method 0)'];
label_for = 'tab:Pareto-all-quantiles-method0-Globe-1880-2023-forecast';
file_for = fullfile(TAB_DIR, ...
    'Table_onemodel_allquantiles_w100_method0_for_Globe_1880_2023.tex');
write_density_table_method0(file_for,caption_for,label_for,rows,F,RDOS_FOR);

caption_fit = ['Fitted confidence intervals for the long-term density ', ...
    'forecast with Pareto-superior models for all quantiles ', ...
    '(Globe, CRU data, 1880--2023, w=100, method 0)'];
label_fit = 'tab:Pareto-all-quantiles-method0-Globe-1880-2023-fit';
file_fit = fullfile(TAB_DIR, ...
    'Table_onemodel_allquantiles_w100_method0_fit_Globe_1880_2023.tex');
write_density_table_method0(file_fit,caption_fit,label_fit,rows,F,RDOS_FIT);

% ------------------------------------------------------------------------
% 6. Save MAT results
% ------------------------------------------------------------------------

save(fullfile(RES_DIR,'RDOS_FOR_allQ_method0_w100_Globe_1880_2023.mat'), ...
    'RDOS_FOR');
save(fullfile(RES_DIR,'RDOS_Fit_allQ_method0_w100_Globe_1880_2023.mat'), ...
    'RDOS_FIT');

fprintf('Density forecast method 0 completed.\n');
fprintf('Tables saved in:  %s\n', TAB_DIR);
fprintf('Results saved in: %s\n', RES_DIR);

% ========================================================================
% Local helper
% ========================================================================

function write_density_table_method0(filename,caption,label,rows,F,RDOS_OUT)
% WRITE_DENSITY_TABLE_METHOD0
% Write LaTeX table for density forecasts and confidence intervals.

fid = fopen(filename,'w');
assert(fid > 0, 'Could not open output file: %s', filename);

fprintf(fid,'\\begin{table}[h!]\n');
fprintf(fid,'\\caption{%s}\n', caption);
fprintf(fid,'\\label{%s}\n', label);
fprintf(fid,'\\begin{center}\n');
fprintf(fid,'\\scalebox{0.8}{%%\n');
fprintf(fid,'\\begin{tabular}{l|ccc} \\hline \\hline\n');
fprintf(fid,' & h=%d & h=%d & h=%d \\\\ \\hline\n', F(1), F(2), F(3));

for i = 1:length(rows)

    point_line = '';
    ci_line = '';

    for j = 1:length(F)
        point_line = [point_line, sprintf('%5.2f&', RDOS_OUT(i,1,1,j))]; %#ok<AGROW>
        ci_line = [ci_line, sprintf('(%5.2f,%5.2f)&', ...
            RDOS_OUT(i,2,1,j), RDOS_OUT(i,3,1,j))]; %#ok<AGROW>
    end

    fprintf(fid,'%s & %s \\\\ \n', rows{i}, point_line(1:end-1));
    fprintf(fid,' & %s \\\\ \n', ci_line(1:end-1));
end

fprintf(fid,'\\hline \\hline\n');
fprintf(fid,'\\end{tabular}}\n');
fprintf(fid,'\\end{center}\n');
fprintf(fid,'\\end{table}\n');

fclose(fid);
end
