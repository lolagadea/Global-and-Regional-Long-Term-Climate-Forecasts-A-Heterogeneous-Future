% ========================================================================
% main_density_forecast_Pareto_allQ_w25_method2_Europe.m
% ------------------------------------------------------------------------
% Density forecasts using Pareto-superior models, all selected quantiles.
% Region: Europe. Window: w = 25. Method 2.
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

FIG_DIR = fullfile(THIS_DIR,'Figures');
RES_DIR = fullfile(THIS_DIR,'Results_for');
TAB_DIR = fullfile(THIS_DIR,'Tables');

if ~isfolder(FIG_DIR), mkdir(FIG_DIR); end
if ~isfolder(RES_DIR), mkdir(RES_DIR); end
if ~isfolder(TAB_DIR), mkdir(TAB_DIR); end

% ------------------------------------------------------------------------
% 2. Load inputs
% ------------------------------------------------------------------------

DATA_FILE = fullfile(INTRO_RES_DIR,'QUANTILES_monthly_1960_2023.mat');
RDOS_FILE = fullfile(RES_DIR,'RDOS_Europe_1960_2023.mat');
BIC_FILE  = fullfile(RES_DIR,'BIC_Europe.dat');

MODEL_FILE1 = fullfile(RES_DIR,'selected_the model_allQ.mat');
MODEL_FILE2 = fullfile(RES_DIR,'selected_the model2_allQ.mat');
ONE_FILE    = fullfile(RES_DIR,'RDOS_one_model.mat');

assert(isfile(DATA_FILE),   'QUANTILES_monthly_1960_2023.mat not found.');
assert(isfile(RDOS_FILE),   'RDOS_Europe_1960_2023.mat not found.');
assert(isfile(BIC_FILE),    'BIC_Europe.dat not found.');
assert(isfile(MODEL_FILE1), 'selected_the model_allQ.mat not found.');
assert(isfile(MODEL_FILE2), 'selected_the model2_allQ.mat not found.');
assert(isfile(ONE_FILE),    'RDOS_one_model.mat not found.');

load(DATA_FILE,'QUANTILES_monthly');
load(RDOS_FILE,'RDOS');
load(MODEL_FILE1,'MODELS');
load(MODEL_FILE2,'MODELS2'); %#ok<NASGU>
load(ONE_FILE,'RDOS_one_model');

QUANTILES = QUANTILES_monthly.Europe;
RDOS_run_comp = RDOS; %#ok<NASGU>

[t,~] = size(QUANTILES);

BIC_all = load(BIC_FILE);

% ------------------------------------------------------------------------
% 3. Settings
% ------------------------------------------------------------------------

models = { ...
    'mean', ...
    'linear-trend', ...
    'pol-trend', ...
    'pol-trend-av-sl', ...
    'pol-trend-log', ...
    'struct-breaks', ...
    'pol-trend-arp', ...
    'pol-trend-arp-av-sl', ...
    'arp', ...
    'rw', ...
    'rwd', ...
    'ima', ...
    'arfima', ...
    'arp20'};

F = [1,10,25];
W = 25;

select_w = 1;
mod = 14;
Q = 9:19;      % q05-q95

SELECTED = NaN(length(F),mod,length(Q)); %#ok<NASGU>

RDOS_FOR = NaN(length(Q),3,length(W),length(F));
RDOS_FIT = NaN(length(Q),3,length(W),length(F));

% ------------------------------------------------------------------------
% 4. Density forecast: Pareto-superior models, method 2
% ------------------------------------------------------------------------

for j = 1:length(W)

    w = W(j);

    for i = 1:length(F)

        f = F(i);

        selected_models_inters = RDOS_one_model.select(:,select_w,i);   % method 2
        selected_models_loss   = (MODELS(:,select_w,i))';              %#ok<NASGU>
        selected_models_first  = RDOS_one_model.firstbest(:,select_w,i); %#ok<NASGU>
        selected_models_second = RDOS_one_model.secondbest(:,select_w,i); %#ok<NASGU>

        for q = 1:length(Q)

            fprintf('Method 2: Europe, h=%d, quantile index=%d\n',f,Q(q));

            Z = QUANTILES(:,Q(q));

            Y = RDOS.Xf(:,Q(q),:,select_w,i);
            Y = squeeze(Y(:,:,1:mod));      % remove combined models
            Y = Y(:,selected_models_inters==1);

            m = size(Y,2);

            % ------------------------------------------------------------
            % Combined model with BIC weights
            % ------------------------------------------------------------

            BIC = BIC_all(Q(q),selected_models_inters==1);

            Weig_sbic = NaN(1,m);

            for k = 1:m
                Weig_sbic(k) = exp(-1/2*BIC(k)) / ...
                    sum(exp(-1/2*BIC(1,:)));
            end

            C1 = NaN(t,1); %#ok<NASGU>
            C1(1:w,:) = Z(1:w,:);
            C1(w+f:t) = Y(w+f:t,:)*Weig_sbic';

            % ------------------------------------------------------------
            % Forecast and confidence intervals
            % ------------------------------------------------------------

            [FORECAST,CI_fit,RES] = compute_forecast_models(Z,f); %#ok<ASGLU>

            FORECAST2 = FORECAST(selected_models_inters==1);

            RES = RES(:,selected_models_inters==1);
            RES = rmmissing(RES);

            FOR_comb = Weig_sbic*FORECAST2;

            var_comb = compute_cov_comb(RES,Weig_sbic);

            CI_fit_comb = FOR_comb + ...
                [norminv(0.05)*sqrt(var_comb), norminv(0.95)*sqrt(var_comb)];

            if length(find(selected_models_inters==1)) == 0

                disp('Warning: there are no selected models.')

            elseif length(find(selected_models_inters==1)) == 1

                varres = RES'*RES/t;
                CI_for_comb = FOR_comb + norminv([0.05 0.95])*sqrt(varres);

            elseif length(find(selected_models_inters==1)) > 1

                [ci_for_sbic, ci_models] = ...
                    compute_CI_comb_select_models( ...
                    QUANTILES,RDOS,Q(q),w,f,selected_models_inters,Weig_sbic,1); %#ok<ASGLU>

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
% 5. Table: forecast intervals
% ------------------------------------------------------------------------

rows = {'q05','q10','q20','q30','q40','q50','q60','q70','q80','q90','q95'};

tableCaption = ['Long-term density forecast with Pareto superior model ' ...
    'for all quantiles (Europe, CRU data, 1960--2023, w=25, method2)'];

tableLabel = 'tab-Pareto-one-model-all-quantiles-method2-Europe-1960-2023';

file = fullfile(TAB_DIR, ...
    'Table_onemodel_allquantiles_w25_method2_for_Europe_1960_2023.tex');

fid = fopen(file,'w');

tt = [' ' '\begin{table}[h!]\caption{' tableCaption '}\label{' ...
    tableLabel '}\begin{center}\scalebox{0.8}{\begin{tabular}{l|ccc} ' ...
    '\hline \hline ' newline];

fprintf(fid,tt);
fprintf(fid,['  & h=1 & h=10 & h=25 \\ \hline ' newline]);

for i = 1:length(rows)

    t1 = [];
    t2 = [];

    for j = 1:length(F)
        t1 = [t1 num2str(RDOS_FOR(i,1,1,j),'%5.2f') '&']; %#ok<AGROW>
        t2 = [t2 ' (' num2str(RDOS_FOR(i,2,1,j),'%5.2f') ',' ...
            num2str(RDOS_FOR(i,3,1,j),'%5.2f') ')' '&']; %#ok<AGROW>
    end

    tt = [rows{i} '&' t1(1:end-1) '\\ ' newline ...
        '&' t2(1:end-1) ' \\' newline];

    fprintf(fid,tt);

end

fprintf(fid,['\hline \hline \end{tabular}}\end{center}\end{table}' newline]);
fclose(fid);

% ------------------------------------------------------------------------
% 6. Table: fitted intervals
% ------------------------------------------------------------------------

file = fullfile(TAB_DIR, ...
    'Table_onemodel_allquantiles_w25_method2_fit_Europe_1960_2023.tex');

fid = fopen(file,'w');

tt = [' ' '\begin{table}[h!]\caption{' tableCaption '}\label{' ...
    tableLabel '}\begin{center}\scalebox{0.8}{\begin{tabular}{l|ccc} ' ...
    '\hline \hline ' newline];

fprintf(fid,tt);
fprintf(fid,['  & h=1 & h=10 & h=25 \\ \hline ' newline]);

for i = 1:length(rows)

    t1 = [];
    t2 = [];

    for j = 1:length(F)
        t1 = [t1 num2str(RDOS_FIT(i,1,1,j),'%5.2f') '&']; %#ok<AGROW>
        t2 = [t2 ' (' num2str(RDOS_FIT(i,2,1,j),'%5.2f') ',' ...
            num2str(RDOS_FIT(i,3,1,j),'%5.2f') ')' '&']; %#ok<AGROW>
    end

    tt = [rows{i} '&' t1(1:end-1) '\\ ' newline ...
        '&' t2(1:end-1) ' \\' newline];

    fprintf(fid,tt);

end

fprintf(fid,['\hline \hline \end{tabular}}\end{center}\end{table}' newline]);
fclose(fid);

% ------------------------------------------------------------------------
% 7. Display selected models
% ------------------------------------------------------------------------

for i = 1:length(F)

    selected_models_inters = RDOS_one_model.select(:,select_w,i);

    fprintf('\nThe selected models with method 2, h=%d, are:\n',F(i));
    disp(models(selected_models_inters==1));

end

% ------------------------------------------------------------------------
% 8. Save outputs
% ------------------------------------------------------------------------

save(fullfile(RES_DIR,'RDOS_FOR_allQ_method2_w25_Europe_1960_2023.mat'), ...
    'RDOS_FOR');

save(fullfile(RES_DIR,'RDOS_FIT_allQ_method2_w25_Europe_1960_2023.mat'), ...
    'RDOS_FIT');

fprintf('\nDensity forecast method2 completed for Europe.\n');

clear global;