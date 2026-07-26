% ========================================================================
% compare_results_regional_portable.m
% ------------------------------------------------------------------------
% Regional forecast comparison.
%
% This script compares forecasted temperature increases across regions,
% quantiles, forecast horizons, and forecasting methods. The increases are
% computed relative to two historical reference periods:
%   - 1986--2005
%   - 1995--2014
%
% The script is portable: all paths are defined relative to the location of
% this file. It assumes the following folder structure:
%
%   Regional/
%       Comparing/
%           compare_results_regional_portable.m
%           Figures/
%           Tables/
%           Results/
%       Introduction/
%           results/
%               QUANTILES_monthly_1960_2023.mat
%       Forecast/
%           <RegionFolder>/
%               Results_for/
%                   RDOS_FOR_allQ_method*_w25_<Region>_1960_2023.mat
%
% Important naming convention:
%   - APC = Arctic
%   - ANC = Antarctic
% ========================================================================

clear; clc;
warning('off');

% ------------------------------------------------------------------------
% 1. Portable paths
% ------------------------------------------------------------------------
THIS_DIR     = fileparts(mfilename('fullpath'));
REGIONAL_DIR = fileparts(THIS_DIR);

FORECAST_DIR = fullfile(REGIONAL_DIR, 'Forecast');
INTRO_RESULTS_DIR = fullfile(REGIONAL_DIR, 'Introduction', 'results');

FIG_DIR = fullfile(THIS_DIR, 'Figures');
TAB_DIR = fullfile(THIS_DIR, 'Tables');
RES_DIR = fullfile(THIS_DIR, 'Results');

DATA_FILE = fullfile(INTRO_RESULTS_DIR, 'QUANTILES_monthly_1960_2023.mat');

if ~exist(FORECAST_DIR, 'dir')
    error('Forecast folder not found: %s', FORECAST_DIR);
end
if ~exist(INTRO_RESULTS_DIR, 'dir')
    error('Introduction results folder not found: %s', INTRO_RESULTS_DIR);
end
if ~exist(DATA_FILE, 'file')
    error('Quantiles file not found: %s', DATA_FILE);
end
if ~exist(FIG_DIR, 'dir'), mkdir(FIG_DIR); end
if ~exist(TAB_DIR, 'dir'), mkdir(TAB_DIR); end
if ~exist(RES_DIR, 'dir'), mkdir(RES_DIR); end

% ------------------------------------------------------------------------
% 2. Settings
% ------------------------------------------------------------------------
regions_print = {'Globe','Arctic','Europe','NorthAmerica','SouthAmerica', ...
                 'Asia','Africa','Australia','Antarctic'};

regions_field = {'Globe','APC','Europe','NorthAm','SouthAm', ...
                 'Asia','Africa','Australia','ANC'};

regions_folder = {'Globe','Arctic','Europe','NorthAm','SouthAm', ...
                  'Asia','Africa','Australia','Antarctic'};

regions_file   = {'Globe','Arctic','Europe','NorthAm','SouthAm', ...
                  'Asia','Africa','Australia','Antarctic'};

quantile_names = {'q05','q10','q20','q30','q40','q50','q60','q70','q80','q90','q95'};
quantile_cols  = 9:19;

reference_names = {'mean 1986-2005', 'mean 1995-2014'};
reference_idx   = {27:46, 36:55};

methods  = 0:2;
horizons = [1 10 25];

Q = numel(quantile_names);
R = numel(regions_print);
P = numel(reference_names);

% ------------------------------------------------------------------------
% 3. Load regional quantiles and construct reference means
% ------------------------------------------------------------------------
load(DATA_FILE, 'QUANTILES_monthly');

QDATA = NaN(size(QUANTILES_monthly.Globe,1), Q, R);
for r = 1:R
    field_name = regions_field{r};
    if ~isfield(QUANTILES_monthly, field_name)
        error('Field %s not found in QUANTILES_monthly.', field_name);
    end
    QDATA(:,:,r) = QUANTILES_monthly.(field_name)(:, quantile_cols);
end

% MATRIX(q, r, p): mean of quantile q in region r for reference period p.
MATRIX = NaN(Q, R, P);
for p = 1:P
    idx = reference_idx{p};
    MATRIX(:,:,p) = squeeze(mean(QDATA(idx,:,:), 1));
end

% ------------------------------------------------------------------------
% 4. Descriptive table
% ------------------------------------------------------------------------
write_latex_regional_table( ...
    MATRIX, quantile_names, regions_print, reference_names, ...
    'Descriptive values of quantiles by regions (CRU data 1960--2023)', ...
    'tab-descriptive-regional-1960-2023', ...
    fullfile(TAB_DIR, 'Table_descriptive_regional_1960_2023.tex'));

% ------------------------------------------------------------------------
% 5. Forecast comparisons: methods 0, 1 and 2
% ------------------------------------------------------------------------
RESULTS = struct();

for m = methods
    FOR = load_regional_forecasts(m, FORECAST_DIR, regions_print, regions_folder, regions_file, MATRIX);

    for h_id = 1:numel(horizons)
        h = horizons(h_id);

        Matrix_mean1 = squeeze(FOR(:,h_id,:,1));  % Q x R
        Matrix_mean2 = squeeze(FOR(:,h_id,:,2));  % Q x R

        RESULTS.(sprintf('method%d',m)).(sprintf('h%d',h)).mean1986_2005 = Matrix_mean1;
        RESULTS.(sprintf('method%d',m)).(sprintf('h%d',h)).mean1995_2014 = Matrix_mean2;

        fig_name = fullfile(FIG_DIR, sprintf('Figure_compare_regional_method%d_h%d', m, h));
        plot_regional_comparison(Matrix_mean1, Matrix_mean2, quantile_names, regions_print, h, fig_name);
    end
end

% ------------------------------------------------------------------------
% 6. Tables for method 2, h = 10 and h = 25
% ------------------------------------------------------------------------
for h = [10 25]
    table_array = NaN(Q, R, P);
    table_array(:,:,1) = RESULTS.method2.(sprintf('h%d',h)).mean1986_2005;
    table_array(:,:,2) = RESULTS.method2.(sprintf('h%d',h)).mean1995_2014;

    write_latex_regional_table( ...
        table_array, quantile_names, regions_print, reference_names, ...
        sprintf('Temperature increase h=%d (method 2)', h), ...
        sprintf('tab-for%d-regional-method2', h), ...
        fullfile(TAB_DIR, sprintf('Table_for%d_regional_method2.tex', h)));
end

% ------------------------------------------------------------------------
% 7. Paper figures: method 2, only relative to 1986--2005
% ------------------------------------------------------------------------
for h = horizons
    Matrix_mean1 = RESULTS.method2.(sprintf('h%d',h)).mean1986_2005;
    fig_name = fullfile(FIG_DIR, sprintf('Figure_compare_regional_method2_h%d_bis', h));
    plot_regional_comparison_single(Matrix_mean1, quantile_names, regions_print, h, fig_name);
end

save(fullfile(RES_DIR, 'regional_comparison_results.mat'), 'RESULTS', 'MATRIX', ...
     'regions_print', 'regions_field', 'regions_folder', 'regions_file', ...
     'quantile_names', 'horizons', 'methods');

fprintf('Regional comparison completed.\n');
fprintf('Figures saved in: %s\n', FIG_DIR);
fprintf('Tables saved in:  %s\n', TAB_DIR);

% ========================================================================
% Local functions
% ========================================================================

function FOR = load_regional_forecasts(method_id, forecast_dir, regions_print, regions_folder, regions_file, reference_means)
% load_regional_forecasts
% ------------------------------------------------------------------------
% Loads RDOS_FOR files for all regions and computes forecasted temperature
% increases relative to the two reference-period means.
%
% Output:
%   FOR(q, h, r, p)
%       q = quantile
%       h = horizon index: 1, 10, 25
%       r = region
%       p = reference period
% ------------------------------------------------------------------------

    R = numel(regions_file);
    Q = size(reference_means, 1);
    P = size(reference_means, 3);
    H = 3;

    FOR = NaN(Q, H, R, P);

    for r = 1:R
        region_dir = fullfile(forecast_dir, regions_folder{r}, 'Results_for');
        file_name  = sprintf('RDOS_FOR_allQ_method%d_w25_%s_1960_2023.mat', method_id, regions_file{r});
        file_path  = fullfile(region_dir, file_name);

        if ~exist(region_dir, 'dir')
            error('Results_for folder not found for %s: %s', regions_print{r}, region_dir);
        end
        if ~exist(file_path, 'file')
            error('Required RDOS file not found: %s', file_path);
        end

        S = load(file_path, 'RDOS_FOR');
        if ~isfield(S, 'RDOS_FOR')
            error('Variable RDOS_FOR not found in file: %s', file_path);
        end

        forecast_values = squeeze(S.RDOS_FOR(:,1,:,:));  % Q x H

        for p = 1:P
            FOR(:,:,r,p) = forecast_values - reference_means(:,r,p);
        end
    end
end

function plot_regional_comparison(Matrix_mean1, Matrix_mean2, quantile_names, regions_print, h, fig_name)
% plot_regional_comparison
% ------------------------------------------------------------------------
% Creates two-panel bar charts comparing regional forecasted temperature
% increases relative to the 1986--2005 and 1995--2014 reference periods.
% ------------------------------------------------------------------------

    Q = numel(quantile_names);
    R = numel(regions_print);

    fig = figure('Visible','off');

    subplot(2,1,1);
    b = bar(Matrix_mean1', 'FaceColor', 'flat');
    for k = 1:Q
        b(k).CData = k;
    end
    b(Q).FaceColor = [1 0 0];
    ylim([-2,6]);
    ylabel('temperature in celsius degrees', 'FontSize', 10);
    set(gca, 'XTick', 1:R, 'XTickLabel', regions_print, 'FontSize', 12);
    title(sprintf('Increase in temperature by quantiles with respect to the mean 1986-2005 (h=%d)', h));

    subplot(2,1,2);
    b = bar(Matrix_mean2', 'FaceColor', 'flat');
    for k = 1:Q
        b(k).CData = k;
    end
    b(Q).FaceColor = [1 0 0];
    ylim([-2,6]);
    leg = legend(quantile_names, 'Location', 'best');
    set(leg, 'Position', [0.919843750139696 0.378577543013232 0.0359374997206032 0.217914432366901], ...
             'FontSize', 7);
    ylabel('temperature in celsius degrees', 'FontSize', 10);
    set(gca, 'XTick', 1:R, 'XTickLabel', regions_print, 'FontSize', 12);
    title(sprintf('Increase in temperature by quantiles with respect to the mean 1995-2014 (h=%d)', h));

    export_regional_figure(fig, fig_name);
    close(fig);
end

function plot_regional_comparison_single(Matrix_mean1, quantile_names, regions_print, h, fig_name)
% plot_regional_comparison_single
% ------------------------------------------------------------------------
% Creates the paper version of the regional comparison figure, using only
% the 1986--2005 reference period.
% ------------------------------------------------------------------------

    Q = numel(quantile_names);
    R = numel(regions_print);

    fig = figure('Visible','off');
    b = bar(Matrix_mean1', 'FaceColor', 'flat');
    for k = 1:Q
        b(k).CData = k;
    end
    b(Q).FaceColor = [1 0 0];

    ylim([-2,6]);
    legend(quantile_names, 'Location', 'best');
    ylabel('temperature in celsius degrees', 'FontSize', 10);
    set(gca, 'XTick', 1:R, 'XTickLabel', regions_print, 'FontSize', 12);
    title(sprintf('Increase in temperature by quantiles with respect to the mean 1986-2005 (h=%d)', h));

    export_regional_figure(fig, fig_name);
    close(fig);
end

function export_regional_figure(fig, fig_name)
% export_regional_figure
% ------------------------------------------------------------------------
% Exports a figure in PNG, PDF, and EPS formats using the same page setup as
% the original script.
% ------------------------------------------------------------------------

    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperSize', [20 16]);
    set(fig, 'PaperPosition', [0 0 20 16]);

    axes_handles = findall(fig, 'Type', 'axes');
    set(axes_handles, 'FontSize', 12);

    for i = 1:numel(axes_handles)
        set(get(axes_handles(i), 'YLabel'), 'FontSize', 10);
    end

    print(fig, '-dpng', fig_name);
    print(fig, '-dpdf', fig_name);
    print(fig, '-deps', fig_name);
end

function write_latex_regional_table(A, quantile_names, regions_print, reference_names, table_caption, table_label, file_path)
% write_latex_regional_table
% ------------------------------------------------------------------------
% Writes a LaTeX table with quantiles by region and reference period.
%
% Input:
%   A(q, r, p)
%       q = quantile
%       r = region
%       p = reference period/statistic
% ------------------------------------------------------------------------

    Q = numel(quantile_names);
    R = numel(regions_print);
    P = numel(reference_names);

    fid = fopen(file_path, 'w');
    if fid == -1
        error('Cannot open output table: %s', file_path);
    end

    fprintf(fid, ' \\begin{table}[h!]\\caption{%s}\\label{%s}\\begin{center}\\scalebox{0.7}{\\begin{tabular}{ll|', table_caption, table_label);
    fprintf(fid, repmat('c', 1, R));
    fprintf(fid, '} \\hline \\hline \n');

    fprintf(fid, 'quantile & statistics');
    for r = 1:R
        fprintf(fid, ' & \\rotatebox{90}{%s}', regions_print{r});
    end
    fprintf(fid, ' \\\\ \\hline \n');

    for q = 1:Q
        for p = 1:P
            if p == 1
                q_label = quantile_names{q};
            else
                q_label = '';
            end

            fprintf(fid, '%s & %s', q_label, reference_names{p});
            for r = 1:R
                fprintf(fid, ' & %5.2f', A(q,r,p));
            end
            fprintf(fid, ' \\\\ \n');
        end
        fprintf(fid, '\\cline{2-%d} \n', R+2);
    end

    fprintf(fid, '\\hline \\hline \\end{tabular}}\\end{center}\\end{table}\n');
    fclose(fid);
end
