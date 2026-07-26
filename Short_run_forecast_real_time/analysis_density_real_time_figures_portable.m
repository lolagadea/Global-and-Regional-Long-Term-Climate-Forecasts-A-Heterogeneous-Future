% ========================================================================
% analysis_density_real_time_figures_portable.m
% ------------------------------------------------------------------------
% Presentation figures for the real-time density forecast exercise.
%
% This script reads the output produced by main_real_time_density_portable:
%
%   Results/RDOS_density_real_time_Globe.mat
%
% and generates two figures:
%   1. Method 0, Method 1, Method 2 and Observed densities (2x2)
%   2. Method 2 and Observed densities only (2x1, paper version)
%
% The script is portable: all paths are defined relative to the location of
% this file.
% ========================================================================

clear; clc;
warning('off');

% ------------------------------------------------------------------------
% 1. Portable paths
% ------------------------------------------------------------------------
THIS_DIR = fileparts(mfilename('fullpath'));

FIG_DIR = fullfile(THIS_DIR, 'Figures');
RES_DIR = fullfile(THIS_DIR, 'Results');

if ~exist(FIG_DIR, 'dir'), mkdir(FIG_DIR); end
if ~exist(RES_DIR, 'dir'), mkdir(RES_DIR); end

DENSITY_FILE = fullfile(RES_DIR, 'RDOS_density_real_time_Globe.mat');

assert(exist(DENSITY_FILE, 'file') == 2, ...
    'Missing file: %s\nRun main_real_time_density_portable first.', ...
    DENSITY_FILE);

% ------------------------------------------------------------------------
% 2. Load real-time density forecasts
% ------------------------------------------------------------------------
load(DENSITY_FILE, 'RDOS_density_real_time_Globe');

FORECAST0 = RDOS_density_real_time_Globe.mod0;
FORECAST1 = RDOS_density_real_time_Globe.mod1;
FORECAST2 = RDOS_density_real_time_Globe.mod2;
OBS       = RDOS_density_real_time_Globe.true;
years_for = RDOS_density_real_time_Globe.years_for;

% Years to display. These indices typically correspond to:
% 2000, 2005, 2010, 2015, 2020, 2024.
D = [1, 6, 11, 16, 21, 24];
year_labels = string(years_for(D));

% ------------------------------------------------------------------------
% 3. Presentation options
% ------------------------------------------------------------------------
USE_FIXED_BW  = true;    % fixed bandwidth improves comparability
SHARE_X_RANGE = false;   % true forces the same x-axis range in all panels

cmap = parula(numel(D));

% Fixed bandwidth using Silverman's rule over all relevant distributions.
BW = [];
if USE_FIXED_BW
    allvals = [FORECAST0(D,:); FORECAST1(D,:); FORECAST2(D,:); OBS(D,:)];
    allvals = allvals(:);
    sig = std(allvals, 'omitnan');
    n   = sum(~isnan(allvals));
    BW  = 1.06 * sig * n^(-1/5);
end

% Optional common x-axis range across all panels.
x_min = [];
x_max = [];
if SHARE_X_RANGE
    global_vals = [FORECAST0(D,:), FORECAST1(D,:), FORECAST2(D,:), OBS(D,:)];
    x_min = floor(min(global_vals, [], 'all')) - 0.5;
    x_max =  ceil(max(global_vals, [], 'all')) + 0.5;
end

% ========================================================================
% Figure 1: Method 0, Method 1, Method 2 and Observed densities
% ========================================================================
FIG_NAME = 'Figure_density_real_time_Globe_2000_2024';

fig = figure('Visible', 'off'); clf
fig.Position = [100 100 1100 850];

tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

plot_density_panel(FORECAST0, 'Method 0', D, cmap, year_labels, BW, ...
    SHARE_X_RANGE, x_min, x_max);
plot_density_panel(FORECAST1, 'Method 1', D, cmap, year_labels, BW, ...
    SHARE_X_RANGE, x_min, x_max);
plot_density_panel(FORECAST2, 'Method 2', D, cmap, year_labels, BW, ...
    SHARE_X_RANGE, x_min, x_max);
plot_density_panel(OBS, 'Observed', D, cmap, year_labels, BW, ...
    SHARE_X_RANGE, x_min, x_max);

lg = legend('NumColumns', 3);
lg.Layout.Tile = 'south';

set(findall(fig, 'Type', 'axes'), 'FontSize', 12);
set(fig, 'PaperUnits', 'centimeters', ...
         'PaperSize', [22 18], ...
         'PaperPosition', [0 0 22 18]);

print(fig, fullfile(FIG_DIR, FIG_NAME), '-dpdf', '-painters');
print(fig, fullfile(FIG_DIR, FIG_NAME), '-dpng', '-r300');
print(fig, fullfile(FIG_DIR, FIG_NAME), '-deps');
close(fig);

% ========================================================================
% Figure 2: Paper version, Method 2 vs Observed
% ========================================================================
FIG_NAME_BIS = 'Figure_density_real_time_Globe_2000_2024_bis';

% For the paper version, use a common x-axis range for direct comparison.
SHARE_X_RANGE_BIS = true;

allvals_bis = [FORECAST2(D,:); OBS(D,:)];
allvals_bis = allvals_bis(:);
sig_bis = std(allvals_bis, 'omitnan');
n_bis   = sum(~isnan(allvals_bis));
BW_BIS  = 1.06 * sig_bis * n_bis^(-1/5);

x_min_bis = floor(min([FORECAST2(D,:), OBS(D,:)], [], 'all')) - 0.5;
x_max_bis =  ceil(max([FORECAST2(D,:), OBS(D,:)], [], 'all')) + 0.5;

fig = figure('Visible', 'off'); clf
fig.Position = [100 100 900 800];

tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

plot_density_panel(FORECAST2, 'Method 2', D, cmap, year_labels, BW_BIS, ...
    SHARE_X_RANGE_BIS, x_min_bis, x_max_bis);
plot_density_panel(OBS, 'Observed', D, cmap, year_labels, BW_BIS, ...
    SHARE_X_RANGE_BIS, x_min_bis, x_max_bis);

lg = legend('NumColumns', 3);
lg.Layout.Tile = 'south';

set(findall(fig, 'Type', 'axes'), 'FontSize', 12);
set(fig, 'PaperUnits', 'centimeters', ...
         'PaperSize', [20 16], ...
         'PaperPosition', [0 0 20 16]);

print(fig, fullfile(FIG_DIR, FIG_NAME_BIS), '-dpdf', '-painters');
print(fig, fullfile(FIG_DIR, FIG_NAME_BIS), '-dpng', '-r300');
print(fig, fullfile(FIG_DIR, FIG_NAME_BIS), '-deps');
close(fig);

fprintf('Real-time density figures completed.\n');
fprintf('Figures saved in: %s\n', FIG_DIR);

% ========================================================================
% Local functions
% ========================================================================
function plot_density_panel(DATA, title_txt, D, cmap, year_labels, BW, ...
                            share_x_range, x_min, x_max)
% plot_density_panel
% ------------------------------------------------------------------------
% Plots kernel density estimates for selected forecast/observed years.
% ------------------------------------------------------------------------

    nexttile; hold on; box on; grid on
    LW = 1.9;

    if share_x_range
        xg = linspace(x_min, x_max, 400);
    else
        local_min = floor(min(DATA(D,:), [], 'all')) - 0.5;
        local_max =  ceil(max(DATA(D,:), [], 'all')) + 0.5;
        xg = linspace(local_min, local_max, 400);
    end

    ymax = 0;
    for ii = 1:numel(D)
        if isempty(BW)
            [fi, xi] = ksdensity(DATA(D(ii),:), xg);
        else
            [fi, xi] = ksdensity(DATA(D(ii),:), xg, 'Bandwidth', BW);
        end
        plot(xi, fi, ...
             'LineWidth', LW, ...
             'Color', cmap(ii,:), ...
             'DisplayName', char(year_labels(ii)));
        ymax = max(ymax, max(fi));
    end

    ylim([0, ymax * 1.08]);
    xlabel('Temperature level (°C)');
    ylabel('Density');
    title(title_txt);
    hold off
end
