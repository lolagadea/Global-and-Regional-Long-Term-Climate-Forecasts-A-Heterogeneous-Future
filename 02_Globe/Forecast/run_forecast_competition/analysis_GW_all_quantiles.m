% ========================================================================
% analysis_GW_all_quantiles.m
% ------------------------------------------------------------------------
% Forecast Pareto competition for all global temperature characteristics.
% Computes Giacomini-White pairwise forecast comparison tests.
%
% Inputs:
%   Results_for/RDOS_Globe_1880_2023.mat
%   02_Globe/Introduction/Results/QUANTILES_monthly_Globe_1880_2023.mat
%
% Outputs:
%   Results_for/GW_rdos_for*_w*_all_quantiles.mat
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

addpath(genpath(FUN_GENERAL_DIR));
addpath(genpath(FUN_FORECAST_DIR));
addpath(genpath(RES_DIR));

% ------------------------------------------------------------------------
% 2. Load inputs
% ------------------------------------------------------------------------

DATA_FILE = fullfile(ROOT_DIR, ...
    '02_Globe', 'Introduction', 'Results', ...
    'QUANTILES_monthly_Globe_1880_2023.mat');

RDOS_FILE = fullfile(RES_DIR, 'RDOS_Globe_1880_2023.mat');

assert(isfile(DATA_FILE), 'Input quantiles file not found.');
assert(isfile(RDOS_FILE), 'RDOS_Globe_1880_2023.mat not found.');

load(DATA_FILE, 'QUANTILES_monthly', 'name_labels');
load(RDOS_FILE, 'RDOS');

% ------------------------------------------------------------------------
% 3. Settings
% ------------------------------------------------------------------------

F = [1, 10, 25, 50];
W = [50, 75, 100];

Z = QUANTILES_monthly.Globe;

S = 1:19;              % characteristics used in the forecast competition
Z = Z(:,S);

[t,n] = size(Z);
names = name_labels; %#ok<NASGU>

m = 14;                % number of individual models, excluding combinations

% ------------------------------------------------------------------------
% 4. GW tests and Pareto model selection
% ------------------------------------------------------------------------

for wi = 1:length(W)
    
    w = W(wi);
    w2 = strcat('w',num2str(w));
    
    for fi = 1:length(F)
        
        f = F(fi);
        f2 = strcat('for',num2str(f));
        
        if w + f < t
            
            Y = RDOS.Xf(:,S,:,wi,fi);
            Y = Y(:,:,1:m); % remove combined models
            
            npairs = m*(m-1)/2;
            GW_test = NaN(npairs,5,n);
            RMSE2 = NaN(n,m);
            
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
                
                teststat = NaN(npairs,1);
                pval = NaN(npairs,1);
                sign = NaN(npairs,1);
                RMSE = NaN(npairs,2);
                
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
                        
                    else
                        [teststat(c),~,pval(c),sign(c)] = ...
                            CPAtest(L1,L2,ff,0.05,1);
                    end
                    
                    RMSE(c,:) = [rmse(pairs(c,1)), rmse(pairs(c,2))];
                    
                end
                
                GW_test(:,:,i) = [RMSE, teststat, pval, sign];
                
            end
            
            A = NaN(npairs,n);
            
            for i = 1:n
                for j = 1:npairs
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
                
                for k = 1:m
                    BEAT(k,1,i) = length(find(B2(k,:,i) == 0));
                    BEAT(k,2,i) = length(find(B2(k,:,i) == 1));
                end
                
            end
            
            models = NaN(n,m);
            models2 = NaN(n,m);
            
            for i = 1:n
                for k = 1:m
                    if BEAT(k,1,i) > 0 && BEAT(k,2,i) == 0
                        models(i,k) = 1;
                    end
                end
                
                M = find_second_bestModel(BEAT(:,:,i));
                models2(i,M') = 1;
            end
            
            RDOS.GW.test = GW_test;
            RDOS.GW.test_matrix = B;
            RDOS.GW.models = models;
            RDOS.GW.models2 = models2;
            RDOS.GW.RMSE = RMSE2;
            
        else
            
            RDOS = 'w and h are incompatibles with t';
            
        end
        
        % ----------------------------------------------------------------
        % 5. Save GW results
        % ----------------------------------------------------------------
        
        OUT_FILE = fullfile(RES_DIR, ...
            strcat('GW_rdos_',f2,'_',w2,'_all_quantiles.mat'));
        
        save(OUT_FILE,'RDOS');
        
        clear Y
        
    end
end

fprintf('\nGW analysis completed.\n');
fprintf('Results saved in:\n%s\n', RES_DIR);