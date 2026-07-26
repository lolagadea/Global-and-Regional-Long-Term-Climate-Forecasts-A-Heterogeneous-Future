% ========================================================================
% main_real_time_density_portable.m
% ------------------------------------------------------------------------
% Real-time density forecast exercise for Globe temperature quantiles.
%
% Portable version adapted to the FINAL_IJF/CODES_IJF structure.
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

if ~exist(FIG_DIR, 'dir'), mkdir(FIG_DIR); end
if ~exist(RES_DIR, 'dir'), mkdir(RES_DIR); end

% ------------------------------------------------------------------------
% 2. Add paths
% ------------------------------------------------------------------------
if exist(FUN_FORECAST_DIR, 'dir') == 7
    addpath(genpath(FUN_FORECAST_DIR));
end

if exist(FUN_GENERAL_DIR, 'dir') == 7
    addpath(genpath(FUN_GENERAL_DIR));
end

if exist(ARFIMA_DIR, 'dir') == 7
    addpath(genpath(ARFIMA_DIR));
end

addpath(genpath(GLOBE_FORECAST_RES));

% ------------------------------------------------------------------------
% 3. Input files
% ------------------------------------------------------------------------
QUANTILES_FILE = fullfile( ...
    GLOBE_INTRO_RES, ...
    'QUANTILES_monthly_Globe_1880_2023.mat');

BIC_FILE = fullfile( ...
    GLOBE_FORECAST_RES, ...
    'BIC.dat');

GW_RDOS_FILE = fullfile( ...
    GLOBE_FORECAST_RES, ...
    'GW_rdos_for1_w100_all_quantiles.mat');

SELECTED_MODEL_FILE = fullfile( ...
    GLOBE_FORECAST_RES, ...
    'selected_the_model_allQ.mat');

SELECTED_MODEL2_FILE = fullfile( ...
    GLOBE_FORECAST_RES, ...
    'selected_the_model2_allQ.mat');

RDOS_ONE_MODEL_FILE = fullfile( ...
    GLOBE_FORECAST_RES, ...
    'RDOS_one_model.mat');

assert(exist(GLOBE_DIR, 'dir') == 7, ...
    'Globe folder not found: %s', GLOBE_DIR);

assert(exist(GLOBE_INTRO_RES, 'dir') == 7, ...
    'Globe introduction results folder not found: %s', GLOBE_INTRO_RES);

assert(exist(GLOBE_FORECAST_RES, 'dir') == 7, ...
    'Globe forecast results folder not found: %s', GLOBE_FORECAST_RES);

assert(exist(QUANTILES_FILE, 'file') == 2, ...
    'Missing file: %s', QUANTILES_FILE);

assert(exist(BIC_FILE, 'file') == 2, ...
    'Missing file: %s', BIC_FILE);

assert(exist(GW_RDOS_FILE, 'file') == 2, ...
    'Missing file: %s', GW_RDOS_FILE);

assert(exist(SELECTED_MODEL_FILE, 'file') == 2, ...
    'Missing file: %s', SELECTED_MODEL_FILE);

assert(exist(SELECTED_MODEL2_FILE, 'file') == 2, ...
    'Missing file: %s', SELECTED_MODEL2_FILE);

assert(exist(RDOS_ONE_MODEL_FILE, 'file') == 2, ...
    'Missing file: %s', RDOS_ONE_MODEL_FILE);

% ------------------------------------------------------------------------
% 4. Load data
% ------------------------------------------------------------------------
load(QUANTILES_FILE, 'QUANTILES_monthly');

QUANTILES = QUANTILES_monthly.Globe;
years     = QUANTILES_monthly.years;
years     = years';

models = {'mean','linear-trend','pol-trend','pol-trend-av-sl', ...
          'pol-trend-log','struct-breaks','pol-trend-arp', ...
          'pol-trend-arp-av-sl','arp','rw','rwd','ima','arfima','arp20'};

BIC = load(BIC_FILE);

[t,n] = size(QUANTILES);
select_w = 3;
select_f = 1;
mod = 14;
Q = 9:19;

f  = 1;
ww = 100;
f2 = strcat('for',num2str(f));
w2 = strcat('w',num2str(ww));

w = 120;
years = [years; years(end)+1]; % We also forecast 2024.

% ========================================================================
% Method 0
% ========================================================================
S = load(GW_RDOS_FILE, 'RDOS');
RDOS = S.RDOS;

FOR = NaN(t-w, mod, length(Q));
for q = 1:length(Q)
    selected_mod0 = RDOS.GW.models(Q(q), 1:mod);
    m = sum(selected_mod0 == 1);
    Y = QUANTILES(:, Q(q));

    for i = 1:t-w+1
        Z = Y(1:w+i-1, 1);
        FOR_years(i) = years(w+i);
        [FORECAST, CI_fit, RES] = compute_forecast_models(Z, f);
        FOR(i,:,q) = FORECAST';
    end

    % Combined model with BIC weights.
    BIC2 = BIC(Q(q), selected_mod0 == 1);
    Weig_sbic = NaN(1, m);
    for k = 1:m
        Weig_sbic(k) = exp(-1/2*BIC2(k)) / sum(exp(-1/2*BIC2(1,:)));
    end
    FORECAST0(:,q) = FOR(:, selected_mod0 == 1, q) * Weig_sbic';
end

figure(1)
[tf] = size(FORECAST0, 1);
for i = 1:tf
    a = my_round(min(FORECAST0(i,:)));
    b = my_round(max(FORECAST0(i,:)));
    pst = a:(b-a)/100:b;
    [fi,xi] = ksdensity(FORECAST0(i,:), pst);
    plot(xi, fi)
    hold on
    A(i) = max(FORECAST0(i,:));
    B(i) = min(FORECAST0(i,:));
end
hold off

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[20 16]);
set(gcf,'PaperPosition',[0 0 20 16]);
set(findall(gcf,'type','axes'),'FontSize',14);

print(gcf, '-dpdf', fullfile(FIG_DIR, 'Figure_density_real_time_mod0_Globe'));
print(gcf, '-dpng', fullfile(FIG_DIR, 'Figure_density_real_time_mod0_Globe'));
print(gcf, '-deps', fullfile(FIG_DIR, 'Figure_density_real_time_mod0_Globe'));

% ========================================================================
% Method 1
% ========================================================================
S1 = load(SELECTED_MODEL_FILE, 'MODELS');
S2 = load(SELECTED_MODEL2_FILE, 'MODELS2');
MODELS  = S1.MODELS;
MODELS2 = S2.MODELS2;

selected_models_loss  = (MODELS(:,select_w,select_f))';  % method 1
selected_models_loss2 = (MODELS2(:,select_w,select_f))';

if length(find(selected_models_loss == 1)) == 0
    disp('Warning: there are no selected models; using selected_models_loss2.')
    selected_mod1 = selected_models_loss2;
else
    selected_mod1 = selected_models_loss;
end

FOR = NaN(t-w, mod, length(Q));
for q = 1:length(Q)
    Y = QUANTILES(:, Q(q));
    for i = 1:t-w+1
        Z = Y(1:w+i-1, 1);
        FOR_years(i) = years(w+i);
        [FORECAST, CI_fit, RES] = compute_forecast_models(Z, f);
        FOR(i,:,q) = FORECAST';
    end
end

FORECAST1 = FOR(:, selected_mod1 == 1, :);
if sum(selected_mod1 == 1)
    FORECAST1 = squeeze(FORECAST1);
else
    % Combined model with BIC weights.
    m = sum(selected_mod1 == 1);
    BIC2 = BIC(Q(q), selected_mod1 == 1);
    Weig_sbic = NaN(1, m);
    for k = 1:m
        Weig_sbic(k) = exp(-1/2*BIC2(k)) / sum(exp(-1/2*BIC2(1,:)));
    end
    FORECAST1 = FOR(:, selected_mod1 == 1, :) * Weig_sbic';
end

figure(2)
[tf] = size(FORECAST1, 1);
for i = 1:tf
    a = my_round(min(FORECAST1(i,:)));
    b = my_round(max(FORECAST1(i,:)));
    pst = a:(b-a)/100:b;
    [fi,xi] = ksdensity(FORECAST1(i,:), pst);
    plot(xi, fi)
    hold on
    A(i) = max(FORECAST1(i,:));
    B(i) = min(FORECAST1(i,:));
end
hold off

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[20 16]);
set(gcf,'PaperPosition',[0 0 20 16]);
set(findall(gcf,'type','axes'),'FontSize',14);

print(gcf, '-dpdf', fullfile(FIG_DIR, 'Figure_density_real_time_mod1_Globe'));
print(gcf, '-dpng', fullfile(FIG_DIR, 'Figure_density_real_time_mod1_Globe'));
print(gcf, '-deps', fullfile(FIG_DIR, 'Figure_density_real_time_mod1_Globe'));

% ========================================================================
% Method 2
% ========================================================================
S = load(RDOS_ONE_MODEL_FILE, 'RDOS_one_model');
RDOS_one_model = S.RDOS_one_model;

selected_mod2 = RDOS_one_model.select(:, select_w, select_f); % method 2

FOR = NaN(t-w, mod, length(Q));
for q = 1:length(Q)
    Y = QUANTILES(:, Q(q));
    for i = 1:t-w+1
        Z = Y(1:w+i-1, 1);
        FOR_years(i) = years(w+i);
        [FORECAST, CI_fit, RES] = compute_forecast_models(Z, f);
        FOR(i,:,q) = FORECAST';
    end
end

FORECAST2 = FOR(:, selected_mod2 == 1, :);
if sum(selected_mod2 == 1)
    FORECAST2 = squeeze(FORECAST2);
else
    % Combined model with BIC weights.
    m = sum(selected_mod2 == 1);
    BIC2 = BIC(Q(q), selected_mod2 == 1);
    Weig_sbic = NaN(1, m);
    for k = 1:m
        Weig_sbic(k) = exp(-1/2*BIC2(k)) / sum(exp(-1/2*BIC2(1,:)));
    end
    FORECAST2 = FOR(:, selected_mod2 == 1, :) * Weig_sbic';
end

figure(3)
[tf] = size(FORECAST2, 1);
for i = 1:tf
    a = my_round(min(FORECAST2(i,:)));
    b = my_round(max(FORECAST2(i,:)));
    pst = a:(b-a)/100:b;
    [fi,xi] = ksdensity(FORECAST2(i,:), pst);
    plot(xi, fi)
    hold on
    A(i) = max(FORECAST2(i,:));
    B(i) = min(FORECAST2(i,:));
end
hold off

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[20 16]);
set(gcf,'PaperPosition',[0 0 20 16]);
set(findall(gcf,'type','axes'),'FontSize',14);

print(gcf, '-dpdf', fullfile(FIG_DIR, 'Figure_density_real_time_mod2_Globe'));
print(gcf, '-dpng', fullfile(FIG_DIR, 'Figure_density_real_time_mod2_Globe'));
print(gcf, '-deps', fullfile(FIG_DIR, 'Figure_density_real_time_mod2_Globe'));

% ========================================================================
% True quantile densities, 2000--2023
% ========================================================================
Z = QUANTILES(w:end-1, 9:19);

figure(4)
for i = 1:tf-1
    a = my_round(min(Z(i,:)));
    b = my_round(max(Z(i,:)));
    pst = a:(b-a)/100:b;
    [fi,xi] = ksdensity(Z(i,:), pst);
    plot(xi, fi)
    hold on
    A(i) = max(Z(i,:));
    B(i) = min(Z(i,:));
end
hold off

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[20 16]);
set(gcf,'PaperPosition',[0 0 20 16]);
set(findall(gcf,'type','axes'),'FontSize',14);

print(gcf, '-dpdf', fullfile(FIG_DIR, 'Figure_density_real_quantiles_2000-2023_Globe'));
print(gcf, '-dpng', fullfile(FIG_DIR, 'Figure_density_real_quantiles_2000-2023_Globe'));
print(gcf, '-deps', fullfile(FIG_DIR, 'Figure_density_real_quantiles_2000-2023_Globe'));

% ========================================================================
% Save results
% ========================================================================
RDOS_density_real_time_Globe.mod0      = FORECAST0;
RDOS_density_real_time_Globe.mod1      = FORECAST1;
RDOS_density_real_time_Globe.mod2      = FORECAST2;
RDOS_density_real_time_Globe.true      = Z;
RDOS_density_real_time_Globe.quantiles = QUANTILES;
RDOS_density_real_time_Globe.years     = years;
RDOS_density_real_time_Globe.years_for = FOR_years;

save(fullfile(RES_DIR, 'RDOS_density_real_time_Globe.mat'), ...
     'RDOS_density_real_time_Globe');

fprintf('Real-time density forecast completed.\n');
fprintf('Figures saved in: %s\n', FIG_DIR);
fprintf('Results saved in: %s\n', RES_DIR);
