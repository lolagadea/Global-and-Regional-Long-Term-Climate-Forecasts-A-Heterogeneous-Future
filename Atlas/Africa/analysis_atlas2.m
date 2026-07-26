% ========================================================================
% analysis_atlas2_portable.m
% ------------------------------------------------------------------------
% Density forecasts vs Atlas projections:
% Africa example.
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

CODES_DIR = fileparts(fileparts(fileparts(THIS_DIR)));

FORECAST_RES_DIR = fullfile(CODES_DIR, ...
    '03 Regional','Forecast','Africa','Results_for');

if ~exist(FIG_DIR,'dir')
    mkdir(FIG_DIR);
end

addpath(RESULTS_DIR);

% ------------------------------------------------------------------------
% 2. Load forecast distributions
% ------------------------------------------------------------------------
load(fullfile(FORECAST_RES_DIR, 'for_2100_Africa_1960_2023.mat'))
for_2100 = FOR_comb.m2.for;

load(fullfile(FORECAST_RES_DIR, ...
    'RDOS_FIT_allQ_method2_w25_Africa_1960_2023.mat'))
for_2050 = RDOS_FIT(:,1,1,3);

% ------------------------------------------------------------------------
% 3. Load Atlas temperature levels
% ------------------------------------------------------------------------
load TEMP_2006_2100_Africa_rcp26
atlas_2048_r26 = TEMP.data.yearly(43);
atlas_2100_r26 = TEMP.data.yearly(95);

load TEMP_2006_2100_Africa_rcp45
atlas_2048_r45 = TEMP.data.yearly(43);
atlas_2100_r45 = TEMP.data.yearly(95);

load TEMP_2006_2100_Africa_rcp85
atlas_2048_r85 = TEMP.data.yearly(43);
atlas_2100_r85 = TEMP.data.yearly(95);


% ------------------------------------------------------------------------
% 3. Make the figure
% ------------------------------------------------------------------------

% --- Colores para RCPs y mediana ---
c_r26 = [0 0.6 0];   % verde
c_r45 = [1 0.5 0];   % naranja
c_r85 = [0.8 0 0];   % rojo
c_med = [0 0 1];     % azul

grid_for = @(x) linspace(floor(min(x))-1, ceil(max(x))+1, 500);

figure(1); clf
tiledlayout(2,1,'TileSpacing','compact','Padding','compact')

% ===== Panel 1: ~2048 (h=25) =====
nexttile
xg = grid_for(for_2050);
[fi,xi] = ksdensity(for_2050, xg);
area(xi, fi, 'FaceAlpha', 0.15, 'FaceColor', [0.7 0.7 0.7], 'LineStyle', 'none'); hold on
plot(xi, fi, 'k', 'LineWidth', 1.5)

xline(atlas_2048_r26, '-', 'RCP2.6', 'LineWidth', 1.6, 'Color', c_r26);
xline(atlas_2048_r45, '-', 'RCP4.5', 'LineWidth', 1.6, 'Color', c_r45);
xline(atlas_2048_r85, '-', 'RCP8.5', 'LineWidth', 1.6, 'Color', c_r85);
xline(median(for_2050,'omitnan'), '-', 'Median', 'LineWidth', 1.6, 'Color', c_med);

ylim([0, max(fi)*1.05]); box on
xlabel('Temperature level (°C)')
ylabel('Density')
%title('Africa — 2048 (h=25)')
title('2048 (h=25)')
legend({'Density forecast','RCP2.6','RCP4.5','RCP8.5','Median'}, ...
       'Location','best')

% ===== Panel 2: 2100 (h=77) =====
nexttile
xg = grid_for(for_2100);
[fi,xi] = ksdensity(for_2100, xg);
area(xi, fi, 'FaceAlpha', 0.15, 'FaceColor', [0.7 0.7 0.7], 'LineStyle', 'none'); hold on
plot(xi, fi, 'k', 'LineWidth', 1.5)

xline(atlas_2100_r26, '-', 'RCP2.6', 'LineWidth', 1.6, 'Color', c_r26);
xline(atlas_2100_r45, '-', 'RCP4.5', 'LineWidth', 1.6, 'Color', c_r45);
xline(atlas_2100_r85, '-', 'RCP8.5', 'LineWidth', 1.6, 'Color', c_r85);
xline(median(for_2100,'omitnan'), '-', 'Median', 'LineWidth', 1.6, 'Color', c_med);

ylim([0, max(fi)*1.05]); box on
xlabel('Temperature level (°C)')
ylabel('Density')
%title('Africa — 2100 (h=77)')
title('2100 (h=77)')
legend({'Density forecast','RCP2.6','RCP4.5','RCP8.5','Median'}, ...
       'Location','best')

% --- Export ---
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[20 16]);
set(gcf,'PaperPosition',[0 0 20 16]);
print(gcf,'Figure_atlas_Africa','-dpdf','-painters');
print(gcf,'Figure_atlas_Africa','-dpng','-r300');
print(gcf,'Figure_atlas_Africa','-deps');