% ========================================================================
% main_real_time_mean_portable.m
% ------------------------------------------------------------------------
% Short-run real-time forecasting exercise for the global mean temperature.
%
% This script reproduces the real-time one-step-ahead forecasting exercise
% for the global mean temperature. Forecasts are generated recursively using
% the same forecasting model set as in the main forecast competition, and
% are then compared with individual benchmarks and forecast combinations.
%
% The script is portable: all paths are defined relative to the location of
% this file. It assumes the script is located in:
%
%   CODES_IJF/04_Short_run_forecast_real_time/Globe/
%
% Required inputs:
%   - 02_Globe/Introduction/Results/QUANTILES_monthly_Globe_1880_2023.mat
%   - 02_Globe/Forecast/run_forecast_competition/Results_for/Weights_BIC.dat
%   - 02_Globe/Forecast/run_forecast_competition/Results_for/BIC.dat
%   - 02_Globe/Forecast/run_forecast_competition/Results_for/GW_rdos_for1_w100_all_quantiles.mat
%
% Outputs:
%   - Figures saved in Figures/
%   - Numerical results saved in Results/real_time_mean_results_Globe.mat
% ========================================================================

clear; clc;
warning('off');

% ------------------------------------------------------------------------
% 1. Portable paths
% ------------------------------------------------------------------------
THIS_DIR = fileparts(mfilename('fullpath'));

% Script location:
% CODES_IJF/04_Short_run_forecast_real_time

SHORT_RUN_DIR = THIS_DIR;

% Root folder:
% CODES_IJF
ROOT_DIR = fileparts(fileparts(SHORT_RUN_DIR));

FUN_FORECAST_DIR = fullfile(ROOT_DIR, 'functions_for');
FUN_GENERAL_DIR  = fullfile(ROOT_DIR, 'functions');
ARFIMA_DIR       = fullfile(ROOT_DIR, 'ARFIMA');

GLOBE_DIR = fullfile(ROOT_DIR, '02_Globe');

GLOBE_INTRO_RES = fullfile( ...
    GLOBE_DIR, 'Introduction', 'Results');

GLOBE_FORECAST_RES = fullfile( ...
    GLOBE_DIR, 'Forecast', ...
    'run_forecast_competition', 'Results_for');

FIG_DIR = fullfile(THIS_DIR, 'Figures');
RES_DIR = fullfile(THIS_DIR, 'Results');

% ------------------------------------------------------------------------
% Add paths
% ------------------------------------------------------------------------
add_if_exists(FUN_FORECAST_DIR);
add_if_exists(FUN_GENERAL_DIR);
add_if_exists(ARFIMA_DIR);

addpath(genpath(GLOBE_FORECAST_RES));

% ------------------------------------------------------------------------
% Input files
% ------------------------------------------------------------------------
QUANTILES_FILE = fullfile( ...
    GLOBE_INTRO_RES, ...
    'QUANTILES_monthly_Globe_1880_2023.mat');

WEIGHTS_FILE = fullfile( ...
    GLOBE_FORECAST_RES, ...
    'Weights_BIC.dat');

BIC_FILE = fullfile( ...
    GLOBE_FORECAST_RES, ...
    'BIC.dat');

GW_RDOS_FILE = fullfile( ...
    GLOBE_FORECAST_RES, ...
    'GW_rdos_for1_w100_all_quantiles.mat');

% ------------------------------------------------------------------------
% Checks
% ------------------------------------------------------------------------
assert(exist(GLOBE_DIR, 'dir') == 7, ...
    'Globe folder not found: %s', GLOBE_DIR);

assert(exist(GLOBE_INTRO_RES, 'dir') == 7, ...
    'Globe introduction results folder not found: %s', ...
    GLOBE_INTRO_RES);

assert(exist(GLOBE_FORECAST_RES, 'dir') == 7, ...
    'Globe forecast results folder not found: %s', ...
    GLOBE_FORECAST_RES);

assert(exist(QUANTILES_FILE, 'file') == 2, ...
    'Missing file: %s', QUANTILES_FILE);

assert(exist(WEIGHTS_FILE, 'file') == 2, ...
    'Missing file: %s', WEIGHTS_FILE);

assert(exist(BIC_FILE, 'file') == 2, ...
    'Missing file: %s', BIC_FILE);

assert(exist(GW_RDOS_FILE, 'file') == 2, ...
    'Missing file: %s', GW_RDOS_FILE);

% ------------------------------------------------------------------------
% 2. Load data and settings
% ------------------------------------------------------------------------
load(QUANTILES_FILE, 'QUANTILES_monthly');

Y     = QUANTILES_monthly.Globe(:,1);   % Global mean temperature
years = QUANTILES_monthly.years(:);
[t,~] = size(Y);

Weights_BIC = load(WEIGHTS_FILE);
BIC         = load(BIC_FILE);

w   = 120;     % Initial estimation window: 1880--1999
f   = 1;       % One-step-ahead forecast
mod = 14;      % Number of forecasting models

% ------------------------------------------------------------------------
% 3. Real-time forecasts with all models
% ------------------------------------------------------------------------
historical = [years, Y];
training   = [years(1:w), Y(1:w)]; %#ok<NASGU>

FOR       = NaN(t-w+1, mod);
FOR_years = NaN(t-w+1, 1);
years_ext = [years; years(end)+1];  % Forecast year after the last observation

for i = 1:t-w+1
    Z = Y(1:w+i-1, 1);
    FOR_years(i) = years_ext(w+i);
    [FORECAST, ~, ~] = compute_forecast_models(Z, f);
    FOR(i,:) = FORECAST(:)';
end

forecast = [FOR_years, FOR];

% Fan plot
fig = figure('Visible','off');
fanplot(historical, forecast, 'NumQuantiles', 10, ...
    'FanLineColor', 'blue', 'HistoricalLineWidth', 1.8, ...
    'ForecastLineColor', 'red');
safe_fontsize(14);
ylabel('Mean temperature');
title({'Forecast mean temperature in real time'});
export_figure(fig, fullfile(FIG_DIR, 'Figure_fan_mean_real_time_Globe'), [20 16], 14);
close(fig);

% RMSE for individual models
RMSE = NaN(1, mod);
for j = 1:mod
    RMSE(j) = sqrt(sum((Y(w+1:t) - forecast(1:end-1,j+1)).^2) / (t-w));
end

% ------------------------------------------------------------------------
% 4. Benchmark, selected model, and forecast combinations
% ------------------------------------------------------------------------
FOR2 = NaN(size(FOR,1), 9);

% Benchmark model
FOR2(:,1) = FOR(:,2);

% Selected model with SBIC: polynomial trend with average slope
FOR2(:,2) = FOR(:,4);

% Combination 0: all models
BIC0        = BIC(1,:);
Weights_BIC = Weights_BIC(1,:);
comb0       = FOR * Weights_BIC';
Weig0       = Weights_BIC; %#ok<NASGU>

% Combination 1: remove two extreme models
A = FOR;
[~, maxpos] = max(mean(A));
[~, minpos] = min(mean(A));
A(:,[minpos,maxpos]) = [];
pos_comb1 = [minpos,maxpos]; %#ok<NASGU>
BIC1 = BIC0;
BIC1([minpos,maxpos]) = [];
Weig1 = exp(-0.5*BIC1) ./ sum(exp(-0.5*BIC1));
comb1 = A * Weig1';

% Combination 2: remove four extreme models
[~, maxpos] = max(mean(A));
[~, minpos] = min(mean(A));
A(:,[minpos,maxpos]) = [];
pos_comb2 = [minpos,maxpos]; %#ok<NASGU>
BIC2 = BIC1;
BIC2([minpos,maxpos]) = [];
Weig2 = exp(-0.5*BIC2) ./ sum(exp(-0.5*BIC2));
comb2 = A * Weig2';

% Combination 3: remove six extreme models
[~, maxpos] = max(mean(A));
[~, minpos] = min(mean(A));
A(:,[minpos,maxpos]) = [];
pos_comb3 = [minpos,maxpos]; %#ok<NASGU>
BIC3 = BIC2;
BIC3([minpos,maxpos]) = [];
Weig3 = exp(-0.5*BIC3) ./ sum(exp(-0.5*BIC3));
comb3 = A * Weig3';

FOR2(:,3) = comb0;
FOR2(:,4) = comb1;
FOR2(:,5) = comb2;
FOR2(:,6) = comb3;

% ------------------------------------------------------------------------
% 5. Pareto-superior forecast combinations
% ------------------------------------------------------------------------
q = 1;   % mean temperature is the first series
S = load(GW_RDOS_FILE, 'RDOS');
selected_models = S.RDOS.GW.models(q,1:mod);
FOR_pareto = FOR(:, selected_models == 1);

comb2_pareto = mean(FOR_pareto, 2);

m = sum(selected_models == 1);
BIC_pareto = BIC(q, selected_models == 1);
Weig_sbic_pareto = exp(-0.5*BIC_pareto) ./ sum(exp(-0.5*BIC_pareto));
comb1_pareto = FOR_pareto * Weig_sbic_pareto';

beta = compute_betas_sum1_all(Y(w+f:t), FOR_pareto(1:end-1,:));
beta = beta(2:end);
comb3_pareto = FOR_pareto * beta;

FOR2(:,7) = comb1_pareto;
FOR2(:,8) = comb2_pareto;
FOR2(:,9) = comb3_pareto;

% ------------------------------------------------------------------------
% 6. Figures
% ------------------------------------------------------------------------
models = {'mean','linear-trend','pol-trend','pol-trend-av-sl', ...
          'pol-trend-log','struct-breaks','pol-trend-arp', ...
          'pol-trend-arp-av-sl','arp','rw','rwd','ima','arfima','arp20'};

options_forecast = {'original','benchmark model','selected sbic model (pol-trend-av-sl)', ...
                    'comb0','comb1','comb2','comb3', ...
                    'comb1-pareto','comb2-pareto','comb3-pareto'};

k = length(options_forecast);

% RMSE figure
fig = figure('Visible','off');
subplot(2,1,1);
bar(RMSE);
safe_fontsize(14);
set(gca, 'XTickLabel', models, 'XTick', 1:mod);
xtickangle(45);
ylim([0 1.2]);
title('RMSE: individual models');

RMSE2 = NaN(1, k-1);
for j = 1:k-1
    RMSE2(j) = sqrt(sum((Y(w+1:t) - FOR2(1:end-1,j)).^2) / (t-w));
end

subplot(2,1,2);
bar(RMSE2);
safe_fontsize(14);
ylim([0 1.2]);
set(gca, 'XTickLabel', options_forecast(2:end), 'XTick', 1:k-1);
xtickangle(45);
title('RMSE: benchmark and forecast combinations');
export_figure(fig, fullfile(FIG_DIR, 'Figure_rmse_mean_real_time_Globe'), [20 16], 14);
close(fig);

% Full real-time forecast comparison
colors = get_forecast_colors();
FOR3 = [NaN(w,k-1); FOR2(1:end-1,:)];

fig = figure('Visible','off');
plot(Y, 'LineWidth', 2, 'Color', [0 0 1]);
safe_fontsize(14);
step = 5;
idx = 1:step:t;
set(gca, 'XTick', idx, 'XTickLabel', years(idx));
xtickangle(45);
hold on;
for i = 1:k-1
    plot(FOR3(:,i), 'LineWidth', 1, 'LineStyle', '--', 'Color', colors{i});
end
legend(options_forecast, 'Location', 'best');
hold off;
export_figure(fig, fullfile(FIG_DIR, 'Figure_for_mean_real_time_Globe'), [20 16], 14);
close(fig);

% Paper figure: selected subset of forecasts
fig = figure('Visible','off');
index_for = [1,2,6,7];
FOR3_bis = FOR3(:, index_for);
k_bis = length(index_for);
options_forecast_bis = options_forecast([1,2,3,7,8]);

plot(Y, 'LineWidth', 2, 'Color', [0 0 1]);
safe_fontsize(14);
set(gca, 'XTick', idx, 'XTickLabel', years(idx));
xtickangle(45);
hold on;
for i = 1:k_bis
    plot(FOR3_bis(:,i), 'LineWidth', 1, 'LineStyle', '--', 'Color', colors{i});
end
legend(options_forecast_bis, 'Location', 'best');
hold off;
export_figure(fig, fullfile(FIG_DIR, 'Figure_for_mean_real_time_Globe_bis'), [20 16], 14);
close(fig);

% ------------------------------------------------------------------------
% 7. Save numerical results
% ------------------------------------------------------------------------
save(fullfile(RES_DIR, 'real_time_mean_results_Globe.mat'), ...
    'Y', 'years', 'w', 'f', 'mod', 'historical', 'forecast', 'FOR', 'FOR2', ...
    'RMSE', 'RMSE2', 'models', 'options_forecast', 'selected_models', ...
    'Weig_sbic_pareto', 'm');

fprintf('Real-time mean forecast exercise completed.\n');
fprintf('Figures saved in: %s\n', FIG_DIR);
fprintf('Results saved in: %s\n', RES_DIR);

% ========================================================================
% Local functions
% ========================================================================

function add_if_exists(dir_path)
    if exist(dir_path, 'dir')
        addpath(genpath(dir_path));
    end
end


function colors = get_forecast_colors()
    colors = { ...
        [1 0 0], ...       % red
        [0 0.6 0], ...     % green
        [0.9 0.7 0], ...   % yellow/orange
        [0.5 0 0.5], ...   % purple
        [0.5 0.5 0.5], ... % grey
        [1 0 1], ...       % magenta
        [0.8 0 0], ...     % dark red
        [0.5 0 0], ...     % maroon
        [0 0 0.5]};        % navy
end

function safe_fontsize(sz)
    try
        fontsize(sz, 'points');
    catch
        set(gca, 'FontSize', sz);
    end
end

function export_figure(fig, file_base, paper_size, font_size)
    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperSize', paper_size);
    set(fig, 'PaperPosition', [0 0 paper_size]);
    set(findall(fig, 'Type', 'axes'), 'FontSize', font_size);
    print(fig, '-dpdf', file_base);
    print(fig, '-dpng', file_base);
    print(fig, '-deps', file_base);
end
