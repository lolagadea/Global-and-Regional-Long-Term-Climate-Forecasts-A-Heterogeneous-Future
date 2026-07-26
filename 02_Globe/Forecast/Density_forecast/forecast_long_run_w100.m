% ========================================================================
% forecast_long_run_w100.m
% ------------------------------------------------------------------------
% Long-run density forecasts for selected temperature quantiles.
%
% This script extends the Method 2 density forecast up to 2100 and produces
% the long-run figures for q05, q50 and q95 relative to alternative reference
% periods.
%
% Inputs
%   - Introduction/Results/QUANTILES_monthly_Globe_1880_2023.mat
%   - Density_forecast/Results/RDOS_FOR_allQ_method2_w100_Globe_1880_2023.mat
%   - run_forecast_competition/Results_for/RDOS_one_model.mat
%   - run_forecast_competition/Results_for/BIC.dat
%
% Outputs
%   - Figures/Figure_Globe_1880_2023_f2024_2100_method2_w100.png/.pdf
%   - Figures/Figure_Globe_1880_2023_f2024_2100_method2_w100_bis.png/.pdf
%
% Notes
%   This version is portable and designed for the IJF replication package.
%   Paths are defined relative to the location of this file.
% ========================================================================

clear; clc;
warning('off');

% ------------------------------------------------------------------------
% 1. Portable paths
% ------------------------------------------------------------------------

THIS_DIR     = fileparts(mfilename('fullpath'));          % Density_forecast
FORECAST_DIR = fileparts(THIS_DIR);                       % Forecast
GLOBE_DIR    = fileparts(FORECAST_DIR);                   % 02_Globe
ROOT_DIR     = fileparts(GLOBE_DIR);                      % CODES_IJF

FUN_GENERAL_DIR  = fullfile(ROOT_DIR,'functions');
FUN_FORECAST_DIR = fullfile(ROOT_DIR,'functions_for');
INTRO_RES_DIR    = fullfile(GLOBE_DIR,'Introduction','Results');
COMP_RES_DIR     = fullfile(FORECAST_DIR,'run_forecast_competition','Results_for');

FIG_DIR = fullfile(THIS_DIR,'Figures');
TAB_DIR = fullfile(THIS_DIR,'Tables'); %#ok<NASGU>
RES_DIR = fullfile(THIS_DIR,'Results');

if ~exist(FIG_DIR,'dir'), mkdir(FIG_DIR); end
if ~exist(TAB_DIR,'dir'), mkdir(TAB_DIR); end
if ~exist(RES_DIR,'dir'), mkdir(RES_DIR); end

if exist(FUN_GENERAL_DIR,'dir'),  addpath(genpath(FUN_GENERAL_DIR));  end
if exist(FUN_FORECAST_DIR,'dir'), addpath(genpath(FUN_FORECAST_DIR)); end

% ------------------------------------------------------------------------
% 2. Inputs and settings
% ------------------------------------------------------------------------

FILE_QUANTILES = fullfile(INTRO_RES_DIR,'QUANTILES_monthly_Globe_1880_2023.mat');
FILE_RDOS_FOR  = fullfile(RES_DIR,'RDOS_FOR_allQ_method2_w100_Globe_1880_2023.mat');
FILE_RDOS_ONE  = fullfile(COMP_RES_DIR,'RDOS_one_model.mat');
FILE_BIC       = fullfile(COMP_RES_DIR,'BIC.dat');

assert(exist(FILE_QUANTILES,'file')==2, 'Required file not found:\n%s', FILE_QUANTILES);
assert(exist(FILE_RDOS_FOR,'file')==2,  'Required file not found:\n%s', FILE_RDOS_FOR);
assert(exist(FILE_RDOS_ONE,'file')==2,  'Required file not found:\n%s', FILE_RDOS_ONE);
assert(exist(FILE_BIC,'file')==2,       'Required file not found:\n%s', FILE_BIC);

load(FILE_QUANTILES,'QUANTILES_monthly');
load(FILE_RDOS_FOR,'RDOS_FOR');
load(FILE_RDOS_ONE,'RDOS_one_model');
BIC = load(FILE_BIC);

QUANTILES = QUANTILES_monthly.Globe;

F        = [1, 10, 25];
h        = length(F);
W        = 100; %#ok<NASGU>
select_w = 3;
Q        = 9:19;

% Mean reference periods:
%   1880-1900, 1986-2005, 1995-2014
MATRIX_mean = [mean(QUANTILES(1:21,  Q))', ...
               mean(QUANTILES(107:126,Q))', ...
               mean(QUANTILES(116:135,Q))'];

% Z stores forecasts for 2024--2100:
% columns 1--3: horizons h = 1, 10, 25
% column 4    : recursive/long-run forecast up to 2100
Z = NaN(length(Q),77,h+1);   % 77 = 2100 - 2023
RDOS_m2 = squeeze(RDOS_FOR(:,1,1,:));

for i = 1:h
    f = F(i);
    Z(:,f,i) = RDOS_m2(:,i);
end

Z_q05 = squeeze(Z(1,:,:));
Z_q50 = squeeze(Z(6,:,:));
Z_q95 = squeeze(Z(11,:,:));

% ------------------------------------------------------------------------
% 3. Long-run forecasts, 2049--2100
% ------------------------------------------------------------------------

selected_models = RDOS_one_model.select(:,select_w,3);   % Method 2
idx_models = (selected_models == 1);

f_long = 26:77;
FOR_comb = NaN(length(f_long),length(Q));

for i = 1:length(f_long)
    for q = 1:length(Q)
        bic_q = BIC(Q(q),idx_models);
        weights_bic = exp(-0.5*bic_q) ./ sum(exp(-0.5*bic_q));

        y_q = QUANTILES(:,Q(q));
        [FORECAST,~,~] = compute_forecast_models(y_q,f_long(i));
        forecast_selected = FORECAST(idx_models);

        FOR_comb(i,q) = weights_bic * forecast_selected;
    end
end

Z_q05(26:end,4) = FOR_comb(:,1);
Z_q50(26:end,4) = FOR_comb(:,6);
Z_q95(26:end,4) = FOR_comb(:,11);

% ------------------------------------------------------------------------
% 4. Temperature increases relative to reference periods
% ------------------------------------------------------------------------

m_q05 = MATRIX_mean(1,:);
m_q50 = MATRIX_mean(6,:);
m_q95 = MATRIX_mean(11,:);

aum_temp_q05_p1 = Z_q05 - m_q05(1);
aum_temp_q05_p2 = Z_q05 - m_q05(2);
aum_temp_q05_p3 = Z_q05 - m_q05(3);

aum_temp_q50_p1 = Z_q50 - m_q50(1);
aum_temp_q50_p2 = Z_q50 - m_q50(2);
aum_temp_q50_p3 = Z_q50 - m_q50(3);

aum_temp_q95_p1 = Z_q95 - m_q95(1);
aum_temp_q95_p2 = Z_q95 - m_q95(2);
aum_temp_q95_p3 = Z_q95 - m_q95(3);

% ------------------------------------------------------------------------
% 5. Figure with three reference periods
% ------------------------------------------------------------------------

time = 2024:2100;
T = length(time);
C = {'b','r','g'};

fig1 = figure(1); clf;
set(fig1,'Color','w');

% q05
subplot(3,1,1);
plot_reference_set(aum_temp_q05_p1,aum_temp_q05_p2,aum_temp_q05_p3,h,C);
set(gca,'XTick',1:T,'XTickLabel',time,'FontSize',12);
ylim([0 4]);
title('q05');

% q50
subplot(3,1,2);
plot_reference_set(aum_temp_q50_p1,aum_temp_q50_p2,aum_temp_q50_p3,h,C);
set(gca,'XTick',1:T,'XTickLabel',time,'FontSize',12);
ylim([0 4]);
title('q50');

% q95
subplot(3,1,3);
plot_reference_set(aum_temp_q95_p1,aum_temp_q95_p2,aum_temp_q95_p3,h,C);
set(gca,'XTick',1:T,'XTickLabel',time,'FontSize',12);
ylim([0 4]);
title('q95');

annotation('textbox',[0.554166666666666 0.00229170305676858 0.04 0.04], ...
    'String',{'1986-2005'},'FitBoxToText','off','EdgeColor',[1 0 0]);
annotation('textbox',[0.452604166666666 0.00398982511923703 0.04 0.04], ...
    'String',{'1880-1900'},'FitBoxToText','off','EdgeColor',[0 0 1]);
annotation('textbox',[0.503124999999999 0.00491703056768563 0.04 0.04], ...
    'String',{'1995-2014'},'FitBoxToText','off','EdgeColor',[0 1 0]);

set(gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0]);
print(fig1,fullfile(FIG_DIR,'Figure_Globe_1880_2023_f2024_2100_method2_w100'),'-dpng');
print(fig1,fullfile(FIG_DIR,'Figure_Globe_1880_2023_f2024_2100_method2_w100'),'-dpdf');

% ------------------------------------------------------------------------
% 6. Figure with two reference periods only: 1880-1900 and 1986-2005
% ------------------------------------------------------------------------

fig2 = figure(2); clf;
set(fig2,'Color','w');

% q05
subplot(3,1,1);
plot_reference_pair(aum_temp_q05_p1,aum_temp_q05_p2,h,C);
set(gca,'XTick',1:T,'XTickLabel',time,'FontSize',12);
ylim([0 4]);
title('q05');

% q50
subplot(3,1,2);
plot_reference_pair(aum_temp_q50_p1,aum_temp_q50_p2,h,C);
set(gca,'XTick',1:T,'XTickLabel',time,'FontSize',12);
ylim([0 4]);
title('q50');

% q95
subplot(3,1,3);
plot_reference_pair(aum_temp_q95_p1,aum_temp_q95_p2,h,C);
set(gca,'XTick',1:T,'XTickLabel',time,'FontSize',12);
ylim([0 4]);
title('q95');

annotation('textbox',[0.554166666666666 0.00229170305676858 0.04 0.04], ...
    'String',{'1986-2005'},'FitBoxToText','off','EdgeColor',[1 0 0]);
annotation('textbox',[0.452604166666666 0.00398982511923703 0.04 0.04], ...
    'String',{'1880-1900'},'FitBoxToText','off','EdgeColor',[0 0 1]);

set(gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0]);
print(fig2,fullfile(FIG_DIR,'Figure_Globe_1880_2023_f2024_2100_method2_w100_bis'),'-dpng');
print(fig2,fullfile(FIG_DIR,'Figure_Globe_1880_2023_f2024_2100_method2_w100_bis'),'-dpdf');

fprintf('Long-run density forecast figures saved in:\n%s\n', FIG_DIR);

% ========================================================================
% Local helper functions
% ========================================================================

function plot_reference_set(aum_p1,aum_p2,aum_p3,h,C)
    for ii = 1:h
        plot(aum_p1(:,ii),'Marker','o','MarkerFaceColor',C{1},'Color',C{1});
        hold on;
    end
    plot(aum_p1(:,h+1),'Color',C{1});

    for ii = 1:h
        plot(aum_p2(:,ii),'Marker','o','MarkerFaceColor',C{2},'Color',C{2});
        hold on;
    end
    plot(aum_p2(:,h+1),'Color',C{2});

    for ii = 1:h
        plot(aum_p3(:,ii),'Marker','o','MarkerFaceColor',C{3},'Color',C{3});
        hold on;
    end
    plot(aum_p3(:,h+1),'Color',C{3});
end

function plot_reference_pair(aum_p1,aum_p2,h,C)
    for ii = 1:h
        plot(aum_p1(:,ii),'Marker','o','MarkerFaceColor',C{1},'Color',C{1});
        hold on;
    end
    plot(aum_p1(:,h+1),'Color',C{1});

    for ii = 1:h
        plot(aum_p2(:,ii),'Marker','o','MarkerFaceColor',C{2},'Color',C{2});
        hold on;
    end
    plot(aum_p2(:,h+1),'Color',C{2});
end
