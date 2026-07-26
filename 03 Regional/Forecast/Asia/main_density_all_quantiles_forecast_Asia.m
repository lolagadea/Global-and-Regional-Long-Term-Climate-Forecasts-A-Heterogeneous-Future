% ========================================================================
% main_density_all_quantiles_forecast_Asia.m
% ------------------------------------------------------------------------
% Selection of Pareto-superior models for the density forecast exercise
% using all selected quantiles for Asia.
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
RES_DIR       = fullfile(THIS_DIR,'Results_for');

if ~isfolder(RES_DIR), mkdir(RES_DIR); end

% ------------------------------------------------------------------------
% 2. Load data and forecast results
% ------------------------------------------------------------------------

DATA_FILE = fullfile(INTRO_RES_DIR,'QUANTILES_monthly_1960_2023.mat');
RDOS_FILE = fullfile(RES_DIR,'RDOS_Asia_1960_2023.mat');

assert(isfile(DATA_FILE), ...
    'QUANTILES_monthly_1960_2023.mat not found.');

assert(isfile(RDOS_FILE), ...
    'RDOS_Asia_1960_2023.mat not found.');

load(DATA_FILE,'QUANTILES_monthly');
load(RDOS_FILE,'RDOS');

% ------------------------------------------------------------------------
% 3. Settings
% ------------------------------------------------------------------------

Q = 1:19;
q = 9:19;          % selected quantiles: q05-q95

W = 25;
F = [1,10,25];

mod = 14;          % individual models only

Z = QUANTILES_monthly.Asia(:,q);
t = size(Z,1);

RMSE   = NaN(mod,length(W),length(F));
BEST   = NaN(length(W),length(F));
MODELS = NaN(mod,length(W),length(F));
MODELS2 = NaN(mod,length(W),length(F));

% ------------------------------------------------------------------------
% 4. Select best model by aggregate RMSE over selected quantiles
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
            'Required GW file not found: %s', GW_FILE);

        GW = load(GW_FILE,'RDOS');

        RMSE(:,j,i) = sum(GW.RDOS.GW.RMSE(q,:))';

        [value,pos] = min(RMSE(:,j,i)); %#ok<ASGLU>
        BEST(j,i) = pos;

    end

end

% ------------------------------------------------------------------------
% 5. GW test using all selected quantiles jointly
% ------------------------------------------------------------------------

pairs = combinator(mod,2,'c');
k = mod*(mod-1)/2;

for j = 1:length(W)

    w = W(j);

    for i = 1:length(F)

        f = F(i);

        fprintf('Density GW tests: Asia, window = %d, horizon = %d\n',w,f);

        GW_test = NaN(k,5);

        if t-w-f < f
            ff = t-w-f+1;
        else
            ff = f;
        end

        Y  = RDOS.Xf(:,q,1:mod,j,i);
        y0 = Z(w+f:t,:);

        L = zeros(t-w-f+1,mod);
        rmse = NaN(1,mod);
        RMSE2 = NaN(1,mod);

        for m = 1:mod
            for v = 1:length(q)
                L(:,m) = L(:,m) + (y0(:,v) - Y(w+f:t,v,m)).^2;
                rmse(m) = sqrt(sum(L(:,m))/(t-w));
                RMSE2(m) = rmse(m);
            end
        end

        RMSE_pairs = NaN(k,2);
        teststat   = NaN(k,1);
        pval       = NaN(k,1);
        signGW     = NaN(k,1);

        for c = 1:length(pairs)

            if L(:,pairs(c,1)) == L(:,pairs(c,2))
                teststat(c) = NaN;
                pval(c)     = NaN;
                signGW(c)   = NaN;

            elseif sum(isnan(L(:,pairs(c,1)))) > 0 || ...
                   sum(isnan(L(:,pairs(c,2)))) > 0
                teststat(c) = NaN;
                pval(c)     = NaN;
                signGW(c)   = NaN;

            elseif isempty(L(:,pairs(c,1))) == 1 || ...
                   isempty(L(:,pairs(c,2))) == 1
                teststat(c) = NaN;
                pval(c)     = NaN;
                signGW(c)   = NaN;

            else
                [teststat(c),~,pval(c),signGW(c)] = ...
                    CPAtest(L(:,pairs(c,1)),L(:,pairs(c,2)),ff,0.05,1);
            end

            RMSE_pairs(c,:) = ...
                [RMSE(pairs(c,1),j,i), RMSE(pairs(c,2),j,i)];

        end

        GW_test = [RMSE_pairs,teststat,pval,signGW];

        % If the test is significant, p-value <= 0.05:
        % sign = 0 means the first model is better;
        % sign = 1 means the second model is better.

        A = NaN(k,1);

        for v = 1:k
            if GW_test(v,4) <= 0.05
                A(v) = GW_test(v,5);
            end
        end

        models = NaN(mod,1);
        M = [];

        if sum(isnan(A)) < k

            B = NaN(mod,mod);
            b = array_to_triang_sup_matrix(A,mod);
            B(:,:) = b;

            B2 = B;

            % Fill the lower triangular matrix.
            % B2(i,j)=0 means that model i is better than model j.
            % B2(i,j)=1 means that model j is better than model i.

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
                BEAT(m,1) = length(find(B2(m,:)==0));
                BEAT(m,2) = length(find(B2(m,:)==1));
            end

            % First-best models
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

save(fullfile(RES_DIR,'selected_the model_allQ.mat'),'MODELS');
save(fullfile(RES_DIR,'selected_the model2_allQ.mat'),'MODELS2');

save(fullfile(RES_DIR,'selected_model_allQ_Asia.mat'), ...
    'MODELS','MODELS2','RMSE','BEST','W','F','q','mod');

clear global;