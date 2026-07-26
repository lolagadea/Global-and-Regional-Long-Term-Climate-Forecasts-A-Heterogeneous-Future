% ========================================================================
% main_density_all_quantiles_forecast.m
% ------------------------------------------------------------------------
% Joint model selection using the selected quantiles q09--q19.
%
% Inputs:
%   Results_for/GW_rdos_for*_w*_all_quantiles.mat
%   02_Globe/Introduction/Results/QUANTILES_monthly_Globe_1880_2023.mat
%
% Outputs:
%   Results_for/selected_the_model_allQ.mat
%   Results_for/selected_the_model2_allQ.mat
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
FIG_MODEL_DIR = fullfile(THIS_DIR, 'Figures_model_selection');

if ~isfolder(RES_DIR), mkdir(RES_DIR); end
if ~isfolder(FIG_MODEL_DIR), mkdir(FIG_MODEL_DIR); end

addpath(genpath(FUN_GENERAL_DIR));
addpath(genpath(FUN_FORECAST_DIR));
addpath(genpath(RES_DIR));

% ------------------------------------------------------------------------
% 2. Load data
% ------------------------------------------------------------------------

DATA_FILE = fullfile(ROOT_DIR, ...
    '02_Globe', 'Introduction', 'Results', ...
    'QUANTILES_monthly_Globe_1880_2023.mat');

assert(isfile(DATA_FILE), ...
    'Input file QUANTILES_monthly_Globe_1880_2023.mat not found.');

load(DATA_FILE, 'QUANTILES_monthly');

% ------------------------------------------------------------------------
% 3. Settings
% ------------------------------------------------------------------------

Q = 1:19; %#ok<NASGU>
q = 9:19;          % selected quantiles

W = [50,75,100];
F = [1,10,25];

mod = 14;

Z = QUANTILES_monthly.Globe(:,q);
t = size(Z,1);

RMSE   = NaN(mod,length(W),length(F));
BEST   = NaN(length(W),length(F));
MODELS = NaN(mod,length(W),length(F));
MODELS2 = NaN(mod,length(W),length(F));

% ------------------------------------------------------------------------
% 4. Select the best model by aggregate RMSE across selected quantiles
% ------------------------------------------------------------------------

for j = 1:length(W)
    
    w = W(j);
    w2 = strcat('w',num2str(w));
    
    for i = 1:length(F)
        
        f = F(i);
        f2 = strcat('for',num2str(f));
        
        GW_FILE = fullfile(RES_DIR, ...
            strcat('GW_rdos_',f2,'_',w2,'_all_quantiles.mat'));
        
        assert(isfile(GW_FILE), ...
            ['Missing GW results file: ', GW_FILE]);
        
        load(GW_FILE, 'RDOS');
        
        RMSE(:,j,i) = sum(RDOS.GW.RMSE(q,:));
        
        [~,pos] = min(RMSE(:,j,i));
        BEST(j,i) = pos;
        
    end
end

% ------------------------------------------------------------------------
% 5. GW test using aggregate loss across selected quantiles
% ------------------------------------------------------------------------

pairs = combinator(mod,2,'c');
k = mod*(mod-1)/2;

for j = 1:length(W)
    
    w = W(j);
    
    for i = 1:length(F)
        
        f = F(i);
        f2 = strcat('for',num2str(f));
        w2 = strcat('w',num2str(w));
        
        GW_FILE = fullfile(RES_DIR, ...
            strcat('GW_rdos_',f2,'_',w2,'_all_quantiles.mat'));
        
        assert(isfile(GW_FILE), ...
            ['Missing GW results file: ', GW_FILE]);
        
        load(GW_FILE, 'RDOS');
        
        GW_test = NaN(k,5);
        RMSE_pairs = NaN(k,2);
        teststat = NaN(k,1);
        pval = NaN(k,1);
        sign = NaN(k,1);
        rmse = NaN(mod,1);
        RMSE2 = NaN(mod,1); %#ok<NASGU>
        
        if t-w-f < f
            ff = t-w-f+1;
        else
            ff = f;
        end
        
        Y = RDOS.Xf(:,q,1:mod,j,i);
        y0 = Z(w+f:t,:);
        
        L = zeros(t-w-f+1,mod);
        
        for m = 1:mod
            for v = 1:length(q)
                L(:,m) = L(:,m) + (y0(:,v) - Y(w+f:t,v,m)).^2;
            end
            
            rmse(m) = sqrt(sum(L(:,m))/(t-w));
            RMSE2(m) = rmse(m);
        end
        
        for c = 1:size(pairs,1)
            
            L1 = L(:,pairs(c,1));
            L2 = L(:,pairs(c,2));
            
            if isequal(L1,L2)
                teststat(c) = NaN;
                pval(c) = NaN;
                sign(c) = NaN;
                
            elseif any(isnan(L1)) || any(isnan(L2))
                teststat(c) = NaN;
                pval(c) = NaN;
                sign(c) = NaN;
                
            elseif isempty(L1) || isempty(L2)
                teststat(c) = NaN;
                pval(c) = NaN;
                sign(c) = NaN;
                
            else
                [teststat(c),~,pval(c),sign(c)] = ...
                    CPAtest(L1,L2,ff,0.05,1);
            end
            
            RMSE_pairs(c,:) = [RMSE(pairs(c,1),j,i), ...
                               RMSE(pairs(c,2),j,i)];
        end
        
        GW_test = [RMSE_pairs, teststat, pval, sign];
        
        % If the test is significant:
        %   sign = 0 means first model is better
        %   sign = 1 means second model is better
        
        A = NaN(k,1);
        
        for v = 1:k
            if GW_test(v,4) <= 0.05
                A(v) = GW_test(v,5);
            end
        end
        
        M = [];
        
        if sum(isnan(A)) == k
            
            models = NaN(mod,1);
            
        else
            
            B = NaN(mod,mod);
            b = array_to_triang_sup_matrix(A,mod);
            B(:,:) = b;
            
            B2 = B;
            
            for mi = 1:mod
                for mj = 1:mod
                    if B(mi,mj) == 1
                        B2(mj,mi) = 0;
                    elseif B(mi,mj) == 0
                        B2(mj,mi) = 1;
                    end
                end
            end
            
            BEAT = NaN(mod,2);
            
            for m = 1:mod
                BEAT(m,1) = length(find(B2(m,:) == 0));
                BEAT(m,2) = length(find(B2(m,:) == 1));
            end
            
            % First-best models
            models = NaN(mod,1);
            
            for m = 1:mod
                if BEAT(m,1) > 0 && BEAT(m,2) == 0
                    models(m) = 1;
                end
            end
            
            % Second-best models
            M = find_second_bestModel(BEAT);
            
        end
        
        MODELS(:,j,i) = models;
        
        if ~isempty(M)
            MODELS2(M',j,i) = 1;
        end
        
    end
end

% ------------------------------------------------------------------------
% 6. Save outputs
% ------------------------------------------------------------------------

save(fullfile(RES_DIR,'selected_the_model_allQ.mat'),  'MODELS');
save(fullfile(RES_DIR,'selected_the_model2_allQ.mat'), 'MODELS2');

fprintf('\nJoint selected-quantiles model selection completed.\n');
fprintf('Results saved in:\n%s\n', RES_DIR);