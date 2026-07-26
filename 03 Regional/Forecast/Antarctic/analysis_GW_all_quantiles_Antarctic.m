% ========================================================================
% analysis_GW_all_quantiles_Antarctic.m
% ------------------------------------------------------------------------
% Giacomini-White pairwise forecast comparison tests for Antarctic.
% Replication package IJF - Regional Forecast module.
%
% Inputs:
%   03_Regional/Introduction/Results/QUANTILES_monthly_1960_2023.mat
%   03_Regional/Forecast/Antarctic/Results/RDOS_Antarctic_1960_2023.mat
%
% Outputs:
%   03_Regional/Forecast/Antarctic/Results_for/GW_rdos_for*_w25_all_quantiles.mat
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
OUT_RES_DIR   = fullfile(THIS_DIR,'Results_for');

if ~isfolder(OUT_RES_DIR)
    mkdir(OUT_RES_DIR);
end

% ------------------------------------------------------------------------
% 2. Load data and forecasts
% ------------------------------------------------------------------------

QUANT_FILE = fullfile(INTRO_RES_DIR,'QUANTILES_monthly_1960_2023.mat');
RDOS_FILE  = fullfile(OUT_RES_DIR,'RDOS_Antarctic_1960_2023.mat');

assert(exist(QUANT_FILE,'file')==2, ...
    'QUANTILES_monthly_1960_2023.mat not found in Regional/Introduction/Results.');

assert(exist(RDOS_FILE,'file')==2, ...
    'RDOS_Antarctic_1960_2023.mat not found in Antarctic/Results.');

load(QUANT_FILE,'QUANTILES_monthly','name_labels');
load(RDOS_FILE,'RDOS');

Z = QUANTILES_monthly.ANC;

S = 1:19;      % fixed order of temperature characteristics
Z = Z(:,S);

[t,n] = size(Z);
names = name_labels;

m = 14;        % individual models only; combinations are excluded

% ------------------------------------------------------------------------
% 3. Forecast horizons and windows
% ------------------------------------------------------------------------

F = [1, 10, 25];
W = 25;

% ------------------------------------------------------------------------
% 4. GW tests
% ------------------------------------------------------------------------

for wi = 1:length(W)

    w  = W(wi);
    w2 = strcat('w',num2str(w));

    for fi = 1:length(F)

        f  = F(fi);
        f2 = strcat('for',num2str(f));

        fprintf('GW tests: Antarctic, window = %d, horizon = %d\n', w, f);

        if w + f < t

            Y = RDOS.Xf(:,S,:,wi,fi);
            Y = Y(:,:,1:m);      % remove combined models

            k = m*(m-1)/2;
            GW_test = NaN(k,5,n);
            RMSE2   = NaN(n,m);

            if t-w-f < f
                ff = t-w-f+1;
            else
                ff = f;
            end

            pairs = combinator(m,2,'c');

            for i = 1:n

                y0 = Z(w+f:t,i);

                L = NaN(length(y0),m);
                rmse = NaN(1,m);

                for j = 1:m
                    L(:,j) = (y0 - Y(w+f:t,i,j)).^2;
                    rmse(j) = sqrt(sum(L(:,j))/(t-w));
                    RMSE2(i,j) = rmse(j);
                end

                RMSE = NaN(k,2);
                teststat = NaN(k,1);
                pval     = NaN(k,1);
                signGW   = NaN(k,1);

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

                    else
                        [teststat(c),~,pval(c),signGW(c)] = ...
                            CPAtest(L(:,pairs(c,1)),L(:,pairs(c,2)),ff,0.05,1);
                    end

                    RMSE(c,:) = [rmse(pairs(c,1)), rmse(pairs(c,2))];

                end

                GW_test(:,:,i) = [RMSE,teststat,pval,signGW];

            end

            % ------------------------------------------------------------
            % Build GW dominance matrices
            % ------------------------------------------------------------

            A = NaN(k,n);

            for i = 1:n
                for j = 1:k
                    if GW_test(j,4,i) <= 0.05
                        A(j,i) = GW_test(j,5,i);
                    end
                end
            end

            B = NaN(m,m,n);

            for i = 1:n
                a = A(:,i);
                b = array_to_triang_sup_matrix(a,m);
                B(:,:,i) = b;
            end

            B2 = B;
            BEAT = NaN(m,2,n);

            for i = 1:n

                for mi = 1:m
                    for mj = 1:m
                        if B(mi,mj,i) == 1
                            B2(mj,mi,i) = 0;
                        elseif B(mi,mj,i) == 0
                            B2(mj,mi,i) = 1;
                        end
                    end
                end

                for kk = 1:m
                    BEAT(kk,1,i) = length(find(B2(kk,:,i)==0));
                    BEAT(kk,2,i) = length(find(B2(kk,:,i)==1));
                end

            end

            models  = NaN(n,m);
            models2 = NaN(n,m);

            for i = 1:n

                for kk = 1:m
                    if BEAT(kk,1,i) > 0 && BEAT(kk,2,i) == 0
                        models(i,kk) = 1;
                    end
                end

                M = find_second_bestModel(BEAT(:,:,i));
                models2(i,M') = 1;

            end

            RDOS.GW.test        = GW_test;
            RDOS.GW.test_matrix = B;
            RDOS.GW.models      = models;
            RDOS.GW.models2     = models2;
            RDOS.GW.RMSE        = RMSE2;

        else

            RDOS.GW.warning = sprintf( ...
                'Window and horizon incompatible with sample size: t=%d, w=%d, h=%d', ...
                t, w, f);

        end

        % ----------------------------------------------------------------
        % 5. Save GW results
        % ----------------------------------------------------------------

        save(fullfile(OUT_RES_DIR, ...
            strcat('GW_rdos_',f2,'_',w2,'_all_quantiles.mat')),'RDOS');

        clear Y

    end
end

clear global;