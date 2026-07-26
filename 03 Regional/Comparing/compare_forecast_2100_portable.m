% ========================================================================
% compare_forecast_2100_portable.m
% ------------------------------------------------------------------------
% Regional comparison of forecasts for year 2100.
%
% This script compares forecasted temperature increases across regions,
% quantiles, and forecasting methods. The increases are computed relative
% to two historical reference periods:
%   - 1986--2005
%   - 1995--2014
%
% The script is portable: all paths are defined relative to the location of
% this file. It assumes the following folder structure:
%
%   03 Regional/
%       Comparing/
%           compare_forecast_2100_portable.m
%           Figures/
%           Tables/
%           Results/
%       Introduction/
%           results/
%               QUANTILES_monthly_1960_2023.mat
%       Forecast/
%           Africa/Results_for/
%           Antarctic/Results_for/
%           Arctic/Results_for/
%           Asia/Results_for/
%           Australia/Results_for/
%           Europe/Results_for/
%           Globe/Results_for/
%           NorthAm/Results_for/
%           SouthAm/Results_for/
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
THIS_DIR      = fileparts(mfilename('fullpath'));
REGIONAL_DIR  = fileparts(THIS_DIR);

FORECAST_DIR      = fullfile(REGIONAL_DIR, 'Forecast');
INTRO_RESULTS_DIR = fullfile(REGIONAL_DIR, 'Introduction', 'results');

FIG_DIR = fullfile(THIS_DIR, 'Figures');
TAB_DIR = fullfile(THIS_DIR, 'Tables');
RES_DIR = fullfile(THIS_DIR, 'Results');

if ~exist(FIG_DIR, 'dir'), mkdir(FIG_DIR); end
if ~exist(TAB_DIR, 'dir'), mkdir(TAB_DIR); end
if ~exist(RES_DIR, 'dir'), mkdir(RES_DIR); end

DATA_FILE = fullfile(INTRO_RESULTS_DIR, 'QUANTILES_monthly_1960_2023.mat');

if ~exist(DATA_FILE, 'file')
    error('Quantiles file not found: %s', DATA_FILE);
end

% ------------------------------------------------------------------------
% 2. Settings
% ------------------------------------------------------------------------
regions_print  = {'Globe','Arctic','Europe','NorthAmerica','SouthAmerica', ...
                  'Asia','Africa','Australia','Antarctic'};

regions_folder = {'Globe','Arctic','Europe','NorthAm','SouthAm', ...
                  'Asia','Africa','Australia','Antarctic'};

regions_file   = {'Globe','Arctic','Europe','NorthAm','SouthAm', ...
                  'Asia','Africa','Australia','Antarctic'};

regions_field  = {'Globe','APC','Europe','NorthAm','SouthAm', ...
                  'Asia','Africa','Australia','ANC'};

quantile_names = {'q05','q10','q20','q30','q40','q50','q60','q70','q80','q90','q95'};
quantile_cols  = 9:19;

reference_names = {'mean 1986-2005', 'mean 1995-2014'};
reference_idx   = {27:46, 36:55};

methods = 0:2;

Q = numel(quantile_names);
R = numel(regions_print);
P = numel(reference_names);

% ------------------------------------------------------------------------
% 3. Load regional quantiles and construct reference means
% ------------------------------------------------------------------------
load(DATA_FILE, 'QUANTILES_monthly');

QDATA = NaN(size(QUANTILES_monthly.Globe, 1), Q, R);

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
% 4. Load 2100 forecasts and compute increases
% ------------------------------------------------------------------------
RESULTS = struct();

for m = methods
    [forecast_2100, ci05_2100, ci95_2100] = load_2100_forecasts( ...
        m, FORECAST_DIR, regions_folder, regions_file);

    aum_temp      = NaN(Q, R, P);
    aum_temp_ci05 = NaN(Q, R, P);
    aum_temp_ci95 = NaN(Q, R, P);

    for p = 1:P
        aum_temp(:,:,p)      = forecast_2100 - MATRIX(:,:,p);
        aum_temp_ci05(:,:,p) = ci05_2100    - MATRIX(:,:,p);
        aum_temp_ci95(:,:,p) = ci95_2100    - MATRIX(:,:,p);
    end

    method_name = sprintf('method%d', m);
    RESULTS.(method_name).forecast_2100 = forecast_2100;
    RESULTS.(method_name).ci05_2100     = ci05_2100;
    RESULTS.(method_name).ci95_2100     = ci95_2100;
    RESULTS.(method_name).increase      = aum_temp;
    RESULTS.(method_name).increase_ci05 = aum_temp_ci05;
    RESULTS.(method_name).increase_ci95 = aum_temp_ci95;

    % Two-panel figure for the two reference periods.
    fig_name = fullfile(FIG_DIR, sprintf('Figure_compare_regional_method%d_2100', m));
    plot_2100_bar_two_periods(aum_temp, quantile_names, regions_print, m, fig_name);

    % LaTeX table for this method.
    write_latex_regional_table( ...
        aum_temp, quantile_names, regions_print, reference_names, ...
        sprintf('Temperature increase in 2100 (method %d)', m), ...
        sprintf('tab-for2100-regional-method%d', m), ...
        fullfile(TAB_DIR, sprintf('Table_for2100_regional_method%d.tex', m)));
end

% ------------------------------------------------------------------------
% 5. Additional CI figures for methods 1 and 2
% ------------------------------------------------------------------------
for m = [1 2]
    method_name = sprintf('method%d', m);

    for p = 1:P
        fig_name = fullfile(FIG_DIR, sprintf( ...
            'Figure_regional_quantiles_ci_period%d_method%d_2100', p, m));

        plot_2100_quantile_ci_panels( ...
            RESULTS.(method_name).increase(:,:,p), ...
            RESULTS.(method_name).increase_ci05(:,:,p), ...
            RESULTS.(method_name).increase_ci95(:,:,p), ...
            quantile_names, regions_print, reference_names{p}, m, fig_name);
    end
end

% ------------------------------------------------------------------------
% 6. Paper figure: method 2, relative to 1986--2005
% ------------------------------------------------------------------------
fig_name = fullfile(FIG_DIR, 'Figure_compare_regional_method2_2100_bis');
plot_2100_bar_single_period( ...
    RESULTS.method2.increase(:,:,1), ...
    quantile_names, regions_print, fig_name);

% ------------------------------------------------------------------------
% 7. Save results
% ------------------------------------------------------------------------
save(fullfile(RES_DIR, 'regional_forecast_2100_results.mat'), ...
     'RESULTS', 'MATRIX', 'regions_print', 'regions_folder', ...
     'regions_file', 'regions_field', 'quantile_names', 'methods');

fprintf('Regional 2100 comparison completed.\n');
fprintf('Figures saved in: %s\n', FIG_DIR);
fprintf('Tables saved in:  %s\n', TAB_DIR);
fprintf('Results saved in: %s\n', RES_DIR);

% ========================================================================
% Local functions
% ========================================================================

function [forecast_2100, ci05_2100, ci95_2100] = load_2100_forecasts( ...
    method_id, forecast_dir, regions_folder, regions_file)
% load_2100_forecasts
% ------------------------------------------------------------------------
% Loads for_2100 files for all regions from:
%   Forecast/<RegionFolder>/Results_for/for_2100_<RegionFile>_1960_2023.mat
%
% Output matrices have dimension Q x R.
% ------------------------------------------------------------------------

    R = numel(regions_folder);

    forecast_2100 = [];
    ci05_2100     = [];
    ci95_2100     = [];

    method_field = sprintf('m%d', method_id);

    for r = 1:R
        results_dir = fullfile(forecast_dir, regions_folder{r}, 'Results_for');

        if ~exist(results_dir, 'dir')
            error('Results_for folder not found for %s: %s', regions_folder{r}, results_dir);
        end

        file_name = sprintf('for_2100_%s_1960_2023.mat', regions_file{r});
        file_path = fullfile(results_dir, file_name);

        if ~exist(file_path, 'file')
            error('Required 2100 forecast file not found: %s', file_path);
        end

        S = load(file_path, 'FOR_comb');

        if ~isfield(S, 'FOR_comb')
            error('Variable FOR_comb not found in file: %s', file_path);
        end

        if ~isfield(S.FOR_comb, method_field)
            error('Field FOR_comb.%s not found in file: %s', method_field, file_path);
        end

        F = S.FOR_comb.(method_field);

        if r == 1
            Q = numel(F.for);
            forecast_2100 = NaN(Q, R);
            ci05_2100     = NaN(Q, R);
            ci95_2100     = NaN(Q, R);
        end

        forecast_2100(:,r) = F.for(:);
        ci05_2100(:,r)     = F.ci.fit(:,1);
        ci95_2100(:,r)     = F.ci.fit(:,2);
    end
end

function plot_2100_bar_two_periods(aum_temp, quantile_names, regions_print, method_id, fig_name)
% plot_2100_bar_two_periods
% ------------------------------------------------------------------------
% Creates the two-panel bar chart for the two reference periods.
% ------------------------------------------------------------------------

    Q = numel(quantile_names);
    R = numel(regions_print);

    fig = figure('Visible','off');

    subplot(2,1,1);
    b = bar(aum_temp(:,:,1)', 'FaceColor', 'flat');
    for k = 1:Q
        b(k).CData = k;
    end
    b(Q).FaceColor = [1 0 0];
    ylim([-9,9]);
    ylabel('temperature in celsius degrees');
    set(gca, 'XTick', 1:R, 'XTickLabel', regions_print, 'FontSize', 12);
    title({'Increase in temperature by quantiles with respect to the mean,', ...
           sprintf('1986--2005 (method %d, year 2100)', method_id)});

    subplot(2,1,2);
    b = bar(aum_temp(:,:,2)', 'FaceColor', 'flat');
    for k = 1:Q
        b(k).CData = k;
    end
    b(Q).FaceColor = [1 0 0];
    ylim([-9,9]);
    leg = legend(quantile_names, 'Location', 'best');
    set(leg, ...
        'Position', [0.919843750139696 0.378577543013232 0.0359374997206032 0.217914432366901], ...
        'FontSize', 7);
    ylabel('temperature in celsius degrees');
    set(gca, 'XTick', 1:R, 'XTickLabel', regions_print, 'FontSize', 12);
    title({'Increase in temperature by quantiles with respect to the mean,', ...
           sprintf('1995--2014 (method %d, year 2100)', method_id)});

    export_regional_figure(fig, fig_name);
    close(fig);
end

function plot_2100_bar_single_period(aum_temp_period1, quantile_names, regions_print, fig_name)
% plot_2100_bar_single_period
% ------------------------------------------------------------------------
% Creates the paper version of the 2100 regional comparison figure, using
% only the 1986--2005 reference period.
% ------------------------------------------------------------------------

    Q = numel(quantile_names);
    R = numel(regions_print);

    fig = figure('Visible','off');

    b = bar(aum_temp_period1', 'FaceColor', 'flat');
    for k = 1:Q
        b(k).CData = k;
    end
    b(Q).FaceColor = [1 0 0];

    ylim([-1,9]);
    ylabel('temperature in celsius degrees');
    set(gca, 'XTick', 1:R, 'XTickLabel', regions_print, 'FontSize', 12);
    title('Increase in temperature by quantiles with respect to the mean 1986--2005, method 2 (2100)');

    export_regional_figure(fig, fig_name);
    close(fig);
end

function plot_2100_quantile_ci_panels(aum_temp, ci05, ci95, quantile_names, regions_print, reference_name, method_id, fig_name)
% plot_2100_quantile_ci_panels
% ------------------------------------------------------------------------
% Creates 11 panels, one for each quantile. Each panel reports the point
% forecast and the lower/upper CI across regions.
% ------------------------------------------------------------------------

    Q = numel(quantile_names);
    R = numel(regions_print);

    fig = figure('Visible','off');

    for q = 1:Q
        subplot(4,3,q);

        A = [aum_temp(q,:); ci05(q,:); ci95(q,:)];
        plot1 = plot(A', ...
            'Marker', '_', ...
            'LineWidth', 2, ...
            'LineStyle', 'none', ...
            'MarkerFaceColor', [1 0 0], ...
            'Color', [1 0 0]);

        set(plot1(1), 'Marker', 'o', 'MarkerFaceColor', 'none', 'Color', [0 0.447 0.741]);
        set(plot1(3), 'MarkerEdgeColor', [1 0 0]);

        xlim([0 R+1]);
        ylim([-2 10]);
        set(gca, 'XTick', 0:R+1, ...
                 'XTickLabel', [{''}, regions_print, {''}], ...
                 'FontSize', 12);
        title(quantile_names{q});
    end

    sgtitle(sprintf('Temperature increase in 2100 with respect to %s (method %d)', reference_name, method_id));

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
