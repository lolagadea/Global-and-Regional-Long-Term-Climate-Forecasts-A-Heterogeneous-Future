% ========================================================================
% make_figure_bar_paper_portable.m
% ------------------------------------------------------------------------
% Atlas comparison: CMIP6 projections vs model-based forecasts.
% Portable version.
% ========================================================================

clear; clc;
warning('off');

% ------------------------------------------------------------------------
% 1. Portable paths
% ------------------------------------------------------------------------
THIS_DIR = fileparts(mfilename('fullpath'));

RESULTS_DIR = fullfile(THIS_DIR,'results');
FIG_DIR     = fullfile(THIS_DIR,'Figures');

if ~exist(FIG_DIR,'dir')
    mkdir(FIG_DIR);
end

addpath(RESULTS_DIR);

% ------------------------------------------------------------------------
% 2. Historical reference periods
% ------------------------------------------------------------------------
load TEMP_CMIP6_historical_mon_185001_201412

pa = mean(TEMP.data.yearly(1:51));      % 1850-1900
p1 = mean(TEMP.data.yearly(137:156));   % 1986-2005
p2 = mean(TEMP.data.yearly(147:165));   % 1995-2014

% ------------------------------------------------------------------------
% 3. CMIP6 scenario projections
% ------------------------------------------------------------------------
load TEMP_2015_2100_CMIP6_Globe_SSP126
t126_2048 = TEMP.data.yearly(34);
t126_2100 = TEMP.data.yearly(86);

load TEMP_2015_2100_CMIP6_Globe_SSP245
t245_2048 = TEMP.data.yearly(34);
t245_2100 = TEMP.data.yearly(86);

load TEMP_2015_2100_CMIP6_Globe_SSP370
t370_2048 = TEMP.data.yearly(34);
t370_2100 = TEMP.data.yearly(86);

load TEMP_2015_2100_CMIP6_Globe_SSP585
t585_2048 = TEMP.data.yearly(34);
t585_2100 = TEMP.data.yearly(86);

matrix = NaN(2,4);

matrix(1,:) = [t126_2048, t245_2048, t370_2048, t585_2048];
matrix(2,:) = [t126_2100, t245_2100, t370_2100, t585_2100];

proy1 = matrix - pa;
proy2 = matrix - p1;
proy3 = matrix - p2;

% ------------------------------------------------------------------------
% 4. Model-based forecasts from the paper
% ------------------------------------------------------------------------
for_w25  = 13.46;
for_2100 = 14.175;

periods = [11.12, 12.17, 12.39];

FOR_w25  = for_w25  - periods;
FOR_2100 = for_2100 - periods;

W25_2048 = [proy1(1,:); proy2(1,:); proy3(1,:)];
W2100    = [proy1(2,:); proy2(2,:); proy3(2,:)];

% ------------------------------------------------------------------------
% 5. Figure: 2048
% ------------------------------------------------------------------------
scenario_titles = {'SSP1-2.6','SSP2-4.5','SSP3-7.0','SSP5-8.5'};
xlabels_ref     = {'1850-1900','1986-2005','1995-2014'};

make_atlas_bar_figure( ...
    W25_2048, FOR_w25, scenario_titles, xlabels_ref, ...
    fullfile(FIG_DIR, 'Figure_atlas_Globe_mean_2048'));

% ------------------------------------------------------------------------
% 6. Figure: 2100
% ------------------------------------------------------------------------
make_atlas_bar_figure( ...
    W2100, FOR_2100, scenario_titles, xlabels_ref, ...
    fullfile(FIG_DIR, 'Figure_atlas_Globe_mean_2100'));

fprintf('Atlas bar figures saved in:\n%s\n', FIG_DIR);

% ========================================================================
% Local function
% ========================================================================
function make_atlas_bar_figure(projections, forecast, scenario_titles, xlabels_ref, out_name)

    figure;

    for s = 1:4
        subplot(2,2,s)

        A = [projections(:,s), forecast(:)];
        bar(A);

        xticks(1:3);
        xticklabels(xlabels_ref);

        title(scenario_titles{s}, 'FontSize', 14);
        legend('CMIP6 projections','model-based forecasts', ...
            'FontSize', 14, ...
            'Location', 'best');

        set(gca, 'FontSize', 14);
    end

    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperSize',[20 16]);
    set(gcf,'PaperPosition',[0 0 20 16]);
    set(findall(gcf,'type','axes'),'FontSize',14);

    print(gcf, out_name, '-dpdf');
    print(gcf, out_name, '-dpng');
    print(gcf, out_name, '-deps');

end