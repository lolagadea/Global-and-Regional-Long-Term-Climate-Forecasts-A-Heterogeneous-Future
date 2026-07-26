% ========================================================================
% main_comb_with_Pareto_q95_w100.m
% ------------------------------------------------------------------------
% Forecast combinations for the global q95 using Pareto-superior models.
% Window: w = 100.
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

FIG_DIR = fullfile(THIS_DIR,'Figures');
TAB_DIR = fullfile(THIS_DIR,'Tables');
RES_DIR = fullfile(THIS_DIR,'Results');

if ~isfolder(FIG_DIR), mkdir(FIG_DIR); end
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

RDOS_COMP_FILE = fullfile(COMP_RES_DIR,'RDOS_Globe_1880_2023.mat');
BIC_FILE       = fullfile(COMP_RES_DIR,'BIC.dat');

assert(isfile(DATA_FILE), ...
    'Input file QUANTILES_monthly_Globe_1880_2023.mat not found.');

assert(isfile(RDOS_COMP_FILE), ...
    'RDOS_Globe_1880_2023.mat not found.');

assert(isfile(BIC_FILE), ...
    'BIC.dat not found.');

load(DATA_FILE,'QUANTILES_monthly');
load(RDOS_COMP_FILE,'RDOS');

QUANTILES = QUANTILES_monthly.Globe;
RDOS_run_comp = RDOS; 
BIC_all = load(BIC_FILE);

% ------------------------------------------------------------------------
% 3. Settings
% ------------------------------------------------------------------------

F = [1,10,25];

W = 100;
select_w = 3;     % w=100 is the third window in the forecast competition

Q = 19;
names_q = {'q95'};

mod = 14;
nmodels_table = 17;

nW=1;

RDOS_df.for = NaN(3,3,length(F),nW,length(Q));
RDOS_df.fit = NaN(3,3,length(F),nW,length(Q));

FORECAST_table = NaN(nmodels_table,length(F)*length(Q));
CI_fit_table   = NaN(nmodels_table,length(F)*length(Q),2);
CI_for_table   = NaN(nmodels_table,length(F)*length(Q),2);

SELECTED = NaN(1,14,length(F));

ifig = 1;

% ------------------------------------------------------------------------
% 4. Forecast combinations
% ------------------------------------------------------------------------

for v = 1:length(Q)
    
    a = 1;
    q = Q(v);
    
    Z = QUANTILES_monthly.Globe(:,q);
    t = length(Z);
    
    w = W;
    w2 = strcat('w',num2str(w));
    
    figure(ifig);
    
    for i = 1:length(F)
        
        f = F(i);
        f2 = strcat('for',num2str(f));
        
        GW_FILE = fullfile(COMP_RES_DIR, ...
            strcat('GW_rdos_',f2,'_',w2,'_all_quantiles.mat'));
        
        assert(isfile(GW_FILE), ...
            ['Missing GW results file: ', GW_FILE]);
        
        load(GW_FILE,'RDOS');
        
        selected_models = RDOS.GW.models(q,1:14);
        SELECTED(:,:,i) = selected_models;
        
        Y = RDOS.Xf(:,q,:,select_w,i);
        Y = squeeze(Y(:,:,1:mod));
        
        Y = Y(:,selected_models==1);
        m = size(Y,2);
        
        % ----------------------------------------------------------------
        % Combined models
        % ----------------------------------------------------------------
        
        % BIC weights
        BIC = BIC_all(q,selected_models==1);
        Weig_sbic = NaN(1,m);
        
        for k = 1:m
            Weig_sbic(k) = exp(-0.5*BIC(k)) / sum(exp(-0.5*BIC));
        end
        
        C1 = NaN(t,1);
        C1(1:w,:) = Z(1:w,:);
        C1(w+f:t) = Y(w+f:t,:)*Weig_sbic';
        
        % Simple mean
        C2 = NaN(t,1);
        C2(1:w,:) = Z(1:w,:);
        C2(w+f:t) = mean(Y(w+f:t,:),2);
        Weig_mean = (1/m)*ones(1,m);
        
        % Estimated weights constrained to sum to one
        C3 = NaN(t,1);
        C3(1:w,:) = Z(1:w,:);
        
        X = Y;
        
        if size(X,2) == 1
            BETA = 1;
        else
            beta = compute_betas_sum1_all(Z(w+f:t,1),X(w+f:t,:));
            BETA = beta(2:end);
        end
        
        C3(w+f:t) = Y(w+f:t,:)*BETA;
        Weig_betas = BETA';
        
        clear global;
        
        % ----------------------------------------------------------------
        % Forecasts and confidence intervals
        % ----------------------------------------------------------------
        
        [FORECAST,CI_fit,RES] = compute_forecast_models(Z,f);
        
        FORECAST2 = FORECAST(selected_models==1);
        RES = RES(:,selected_models==1);
        
        % Combination 1: BIC weights
        FOR_comb1 = Weig_sbic*FORECAST2;
        
        [ci_for_sbic,ci_models] = ...
            compute_CI_comb_select_models(QUANTILES,RDOS,q,w,f, ...
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
            compute_CI_comb_select_models(QUANTILES,RDOS,q,w,f, ...
            selected_models,Weig_mean,0);
        
        CI_for_comb2 = FOR_comb2 + ci_for_mean;
        
        var_comb2 = compute_cov_comb(RES,Weig_mean);
        CI_fit_comb2 = FOR_comb2 + ...
            [norminv(0.05)*sqrt(var_comb2), ...
             norminv(0.95)*sqrt(var_comb2)];
        
        % Combination 3: estimated weights
        FOR_comb3 = FORECAST2'*BETA;
        
        ci_for_betas = ...
            compute_CI_comb_select_models(QUANTILES,RDOS,q,w,f, ...
            selected_models,Weig_betas,0);
        
        CI_for_comb3 = FOR_comb3 + ci_for_betas;
        
        var_comb3 = compute_cov_comb(RES,Weig_betas);
        CI_fit_comb3 = FOR_comb3 + ...
            [norminv(0.05)*sqrt(var_comb3), ...
             norminv(0.95)*sqrt(var_comb3)];
        
        % ----------------------------------------------------------------
        % Store tables
        % ----------------------------------------------------------------
        
        FORECAST_table(:,v+(i-1)) = ...
            [FORECAST;FOR_comb1;FOR_comb2;FOR_comb3];
        
        CI_fit_table(:,v+(i-1),1) = ...
            [CI_fit(:,1);CI_fit_comb1(1);CI_fit_comb2(1);CI_fit_comb3(1)];
        
        CI_fit_table(:,v+(i-1),2) = ...
            [CI_fit(:,2);CI_fit_comb1(2);CI_fit_comb2(2);CI_fit_comb3(2)];
        
        pos = find(selected_models==1);
        
        CI_for_table(pos',i,1) = CI_for(:,1);
        CI_for_table(mod+1:nmodels_table,i,1) = ...
            [CI_for_comb1(1);CI_for_comb2(1);CI_for_comb3(1)];
        
        CI_for_table(pos',i,2) = CI_for(:,2);
        CI_for_table(mod+1:nmodels_table,i,2) = ...
            [CI_for_comb1(2);CI_for_comb2(2);CI_for_comb3(2)];
        
        F_comb = [FOR_comb1;FOR_comb2;FOR_comb3];
        CI_for_comb = [CI_for_comb1;CI_for_comb2;CI_for_comb3];
        CI_fit_comb = [CI_fit_comb1;CI_fit_comb2;CI_fit_comb3];
        
        % ----------------------------------------------------------------
        % Figure
        % ----------------------------------------------------------------
        
        subplot(length(F),1,a);
        
        years = 1880:2023+f;
        mcomb = length(F_comb);
        V = NaN(t+f,mcomb);
        
        for k = 1:mcomb
            V(t+f,k) = F_comb(k,1);
        end
        
        C = {'r','g','y'};
        
        plot([Z;NaN(f,1)],'Color','b');
        hold on;
        
        for k = 1:mcomb
            plot(V(:,k), ...
                'Marker','o', ...
                'MarkerFaceColor',C{k}, ...
                'Color',C{k});
            hold on;
        end
        
        hold off;
        
        legend('original','comb1','comb2','comb3','Location','best');
        
        if f == 1
            xt = [1 20 40 60 80 100 120 t+f];
        else
            xt = [1 20 40 60 80 100 120 t t+f];
        end
        
        set(gca,'XTick',xt,'FontSize',12);
        set(gca,'XTickLabel',years(xt),'FontSize',10);
        
        title(strcat('Forecast horizon',{' '}, ...
            'h=',num2str(f),{' '}, ...
            'and window',{' '}, ...
            'w=',num2str(w)));
        
        a = a+1;
        
        % ----------------------------------------------------------------
        % Display values and save temporary results
        % ----------------------------------------------------------------
        
        disp('The forecast values of combined models are:')
        F_comb
        
        disp('Their CI with forecast errors are:')
        CI_for_comb
        
        disp('Their CI with fitted errors are:')
        CI_fit_comb
        
        RDOS_df.fit(:,:,i,1,v) = [F_comb,CI_fit_comb];
        RDOS_df.for(:,:,i,1,v) = [F_comb,CI_for_comb];
        
        clear Y Weig_sbic Weig_mean Weig_betas BETA;
        
    end
    
    set(gcf,'PaperSize',[29.7 21.0], ...
        'PaperPosition',[0 0 29.7 21.0]);
    
    figfile = fullfile(FIG_DIR, ...
        'Figure_for_q95_comb_Pareto_Globe_1880_2023_w100');
    
    print(gcf,figfile,'-dpdf');
    print(gcf,figfile,'-dpng');
    print(gcf,figfile,'-deps');
    
    ifig = ifig+1;
    
end

% ------------------------------------------------------------------------
% 5. Remove non-selected individual models from the table
% ------------------------------------------------------------------------

S = [SELECTED(:,:,1)',SELECTED(:,:,2)',SELECTED(:,:,3)'];

for i = 1:14
    for j = 1:3
        if isnan(S(i,j))
            FORECAST_table(i,j) = NaN;
            CI_fit_table(i,j,1) = NaN;
            CI_fit_table(i,j,2) = NaN;
        end
    end
end

% ------------------------------------------------------------------------
% 6. LaTeX tables
% ------------------------------------------------------------------------

models = {'mean','linear-trend','pol-trend','pol-trend-av-sl', ...
    'pol-trend-log','struct-breaks','pol-trend-arp', ...
    'pol-trend-arp-av-sl','arp','rw','rwd','ima','arfima','arp20', ...
    'combined 1','combined 2','combined 3'};

% Forecast-error intervals
TAB_FILE_FOR = fullfile(TAB_DIR, ...
    'Table_forecast_Pareto_q95_ci_for_w100_Globe_1880_2023.tex');

fid = fopen(TAB_FILE_FOR,'w');

fprintf(fid,'\\begin{tabular}{lccc}\n');
fprintf(fid,'\\toprule\n');
fprintf(fid,'Models/horizon & h=1 & h=10 & h=25 \\\\\n');
fprintf(fid,'\\midrule\n');

for i = 1:length(models)
    
    fprintf(fid,'%s',models{i});
    
    for j = 1:length(F)
        fprintf(fid,' & %5.2f',FORECAST_table(i,j));
    end
    
    fprintf(fid,' \\\\\n');
    
    fprintf(fid,'');
    
    for j = 1:length(F)
        fprintf(fid,' & (%5.2f,%5.2f)', ...
            CI_for_table(i,j,1),CI_for_table(i,j,2));
    end
    
    fprintf(fid,' \\\\\n');
    
end

fprintf(fid,'\\bottomrule\n');
fprintf(fid,'\\end{tabular}\n');

fclose(fid);

% Fitted-error intervals
TAB_FILE_FIT = fullfile(TAB_DIR, ...
    'Table_forecast_Pareto_q95_ci_fit_w100_Globe_1880_2023.tex');

fid = fopen(TAB_FILE_FIT,'w');

fprintf(fid,'\\begin{tabular}{lccc}\n');
fprintf(fid,'\\toprule\n');
fprintf(fid,'Models/horizon & h=1 & h=10 & h=25 \\\\\n');
fprintf(fid,'\\midrule\n');

for i = 1:length(models)
    
    fprintf(fid,'%s',models{i});
    
    for j = 1:length(F)
        fprintf(fid,' & %5.2f',FORECAST_table(i,j));
    end
    
    fprintf(fid,' \\\\\n');
    
    fprintf(fid,'');
    
    for j = 1:length(F)
        fprintf(fid,' & (%5.2f,%5.2f)', ...
            CI_fit_table(i,j,1),CI_fit_table(i,j,2));
    end
    
    fprintf(fid,' \\\\\n');
    
end

fprintf(fid,'\\bottomrule\n');
fprintf(fid,'\\end{tabular}\n');

fclose(fid);

% ------------------------------------------------------------------------
% 7. Save results
% ------------------------------------------------------------------------

RDOS_FILE = fullfile(RES_DIR,'RDOS_for_q95_Globe_1880_2023.mat');

if isfile(RDOS_FILE)
    Sload = load(RDOS_FILE);
    
    if isfield(Sload,'RDOS')
        RDOS_q95 = Sload.RDOS;
    else
        RDOS_q95 = struct();
    end
else
    RDOS_q95 = struct();
end

RDOS_q95.comb_Pareto.H = F;
RDOS_q95.comb_Pareto.W = W;
RDOS_q95.comb_Pareto.models = models;
RDOS_q95.comb_Pareto.forecast = FORECAST_table;
RDOS_q95.comb_Pareto.ci.fit = CI_fit_table;
RDOS_q95.comb_Pareto.ci.for = CI_for_table;
RDOS_q95.comb_Pareto.selected = SELECTED;
RDOS_q95.comb_Pareto.RDOS_df = RDOS_df;

RDOS = RDOS_q95;

save(RDOS_FILE,'RDOS');

fprintf('\nPareto-combined forecasts for the global q95 completed.\n');
fprintf('Tables saved in:\n%s\n', TAB_DIR);
fprintf('Figures saved in:\n%s\n', FIG_DIR);
fprintf('Results saved in:\n%s\n', RES_DIR);