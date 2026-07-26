% ========================================================================
% analysis_atlas2_portable.m
% ------------------------------------------------------------------------
% Density forecasts vs Atlas projections: Australia.
% Portable version.
% ========================================================================

clear; clc;
warning('off');

THIS_DIR = fileparts(mfilename('fullpath'));

RESULTS_DIR = fullfile(THIS_DIR,'results');
FIG_DIR     = fullfile(THIS_DIR,'Figures');

if ~exist(FIG_DIR,'dir')
    mkdir(FIG_DIR);
end

addpath(RESULTS_DIR);

CODES_DIR = fileparts(fileparts(fileparts(THIS_DIR)));

FORECAST_RES_DIR = fullfile(CODES_DIR, ...
    '03 Regional','Forecast','Australia','Results_for');

load(fullfile(FORECAST_RES_DIR, ...
    'for_2100_Australia_1960_2023.mat'))
for_2100 = FOR_comb.m2.for;

load(fullfile(FORECAST_RES_DIR, ...
    'RDOS_FIT_allQ_method2_w25_Australia_1960_2023.mat'))
for_2050 = RDOS_FIT(:,1,1,3);

load TEMP_2006_2100_Australia_rcp26
atlas_2048_r26 = TEMP.data.yearly(43);
atlas_2100_r26 = TEMP.data.yearly(95);

load TEMP_2006_2100_Australia_rcp45
atlas_2048_r45 = TEMP.data.yearly(43);
atlas_2100_r45 = TEMP.data.yearly(95);

load TEMP_2006_2100_Australia_rcp85
atlas_2048_r85 = TEMP.data.yearly(43);
atlas_2100_r85 = TEMP.data.yearly(95);

c_r26 = [0 0.6 0];
c_r45 = [1 0.5 0];
c_r85 = [0.8 0 0];
c_med = [0 0 1];

grid_for = @(x) linspace(floor(min(x))-1, ceil(max(x))+1, 500);

figure(1); clf
tiledlayout(2,1,'TileSpacing','compact','Padding','compact')

nexttile
xg = grid_for(for_2050);
[fi,xi] = ksdensity(for_2050, xg);

area(xi, fi, 'FaceAlpha', 0.15, 'FaceColor', [0.7 0.7 0.7], ...
    'LineStyle', 'none');
hold on
plot(xi, fi, 'k', 'LineWidth', 1.5)

xline(atlas_2048_r26, '-', 'RCP2.6', 'LineWidth', 1.6, 'Color', c_r26);
xline(atlas_2048_r45, '-', 'RCP4.5', 'LineWidth', 1.6, 'Color', c_r45);
xline(atlas_2048_r85, '-', 'RCP8.5', 'LineWidth', 1.6, 'Color', c_r85);
xline(median(for_2050,'omitnan'), '-', 'Median', 'LineWidth', 1.6, 'Color', c_med);

ylim([0, max(fi)*1.05]); box on
xlabel('Temperature level (°C)')
ylabel('Density')
title('2048 (h=25)')
legend({'Density forecast','RCP2.6','RCP4.5','RCP8.5','Median'}, ...
    'Location','best')

nexttile
xg = grid_for(for_2100);
[fi,xi] = ksdensity(for_2100, xg);

area(xi, fi, 'FaceAlpha', 0.15, 'FaceColor', [0.7 0.7 0.7], ...
    'LineStyle', 'none');
hold on
plot(xi, fi, 'k', 'LineWidth', 1.5)

xline(atlas_2100_r26, '-', 'RCP2.6', 'LineWidth', 1.6, 'Color', c_r26);
xline(atlas_2100_r45, '-', 'RCP4.5', 'LineWidth', 1.6, 'Color', c_r45);
xline(atlas_2100_r85, '-', 'RCP8.5', 'LineWidth', 1.6, 'Color', c_r85);
xline(median(for_2100,'omitnan'), '-', 'Median', 'LineWidth', 1.6, 'Color', c_med);

ylim([0, max(fi)*1.05]); box on
xlabel('Temperature level (°C)')
ylabel('Density')
title('2100 (h=77)')
legend({'Density forecast','RCP2.6','RCP4.5','RCP8.5','Median'}, ...
    'Location','best')

OUT_NAME = fullfile(FIG_DIR,'Figure_atlas_Australia');

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[20 16]);
set(gcf,'PaperPosition',[0 0 20 16]);

print(gcf, OUT_NAME, '-dpdf','-painters');
print(gcf, OUT_NAME, '-dpng','-r300');
print(gcf, OUT_NAME, '-deps');

fprintf('Saved:\n%s\n', OUT_NAME);