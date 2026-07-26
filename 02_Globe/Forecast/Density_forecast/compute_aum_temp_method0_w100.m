% ========================================================================
% compute_aum_temp_method0_w100.m
% ------------------------------------------------------------------------
% Compute temperature increases implied by Method 0 density forecasts and
% construct the long-run temperature-increase figures used in the IJF paper.
%
% Method 0:
% Pareto-superior models are selected independently for each quantile.
%
% Globe, CRU / CRUTEM monthly quantiles, 1880--2023
% Window: w = 100
% Long-run forecast horizon: year 2100
%
% This version is portable and writes outputs to:
%   Density_forecast/Figures
%   Density_forecast/Tables
%   Density_forecast/Results
% ========================================================================

clear; clc;
warning('off');

% ------------------------------------------------------------------------
% 1. Portable paths
% ------------------------------------------------------------------------

THIS_DIR     = fileparts(mfilename('fullpath'));
FORECAST_DIR = fileparts(THIS_DIR);
GLOBE_DIR    = fileparts(FORECAST_DIR);
ROOT_DIR     = fileparts(GLOBE_DIR);

FUN_GENERAL_DIR  = fullfile(ROOT_DIR,'functions');
FUN_FORECAST_DIR = fullfile(ROOT_DIR,'functions_for');

COMP_RES_DIR  = fullfile(FORECAST_DIR,'run_forecast_competition','Results_for');
INTRO_RES_DIR = fullfile(GLOBE_DIR,'Introduction','Results');

FIG_DIR = fullfile(THIS_DIR,'Figures');
TAB_DIR = fullfile(THIS_DIR,'Tables');
RES_DIR = fullfile(THIS_DIR,'Results');

if ~exist(FIG_DIR,'dir'), mkdir(FIG_DIR); end
if ~exist(TAB_DIR,'dir'), mkdir(TAB_DIR); end
if ~exist(RES_DIR,'dir'), mkdir(RES_DIR); end

addpath(genpath(FUN_GENERAL_DIR));
addpath(genpath(FUN_FORECAST_DIR));

% ------------------------------------------------------------------------
% 2. Load inputs
% ------------------------------------------------------------------------

FILE_RDOS_FOR = fullfile(RES_DIR,'RDOS_FOR_allQ_method0_w100_Globe_1880_2023.mat');
FILE_QUANT    = fullfile(INTRO_RES_DIR,'QUANTILES_monthly_Globe_1880_2023.mat');
FILE_GW       = fullfile(COMP_RES_DIR,'GW_rdos_for25_w100_all_quantiles.mat');
FILE_BIC      = fullfile(COMP_RES_DIR,'BIC.dat');

assert(exist(FILE_RDOS_FOR,'file')==2, ...
    'Missing file: %s. Run main_density_forecast_Pareto_allQ_w100_method0.m first.', FILE_RDOS_FOR);
assert(exist(FILE_QUANT,'file')==2, ...
    'Missing file: %s. Run the Introduction quantile construction first.', FILE_QUANT);
assert(exist(FILE_GW,'file')==2, ...
    'Missing file: %s. Run the forecast competition first.', FILE_GW);
assert(exist(FILE_BIC,'file')==2, ...
    'Missing file: %s. Run the forecast competition first.', FILE_BIC);

load(FILE_RDOS_FOR,'RDOS_FOR');
load(FILE_QUANT,'QUANTILES_monthly');
load(FILE_GW,'RDOS');
BIC_all = load(FILE_BIC);

QUANTILES = QUANTILES_monthly.Globe;

% ------------------------------------------------------------------------
% 3. Settings
% ------------------------------------------------------------------------

Q = 9:19;                         % q05, q10, ..., q95
q_labels = {'q05','q10','q20','q30','q40','q50', ...
            'q60','q70','q80','q90','q95'};

period_labels = {'1880-1900','1986-2002','1995-2014'};
h_labels = [1,10,25,77];          % 77 corresponds to forecast to 2100
f_long = 77;                      % 2023 + 77 = 2100
n_models = 14;

% Historical reference means for the selected quantiles
MATRIX_mean = [ ...
    mean(QUANTILES(1:21,Q),1)', ...
    mean(QUANTILES(107:126,Q),1)', ...
    mean(QUANTILES(116:135,Q),1)' ];

% ------------------------------------------------------------------------
% 4. Temperature increases implied by density forecasts
% ------------------------------------------------------------------------

RDOS_m0 = squeeze(RDOS_FOR);

aum_temp      = NaN(length(Q),length(period_labels),length(h_labels));
aum_temp_ci05 = NaN(length(Q),length(period_labels),length(h_labels));
aum_temp_ci95 = NaN(length(Q),length(period_labels),length(h_labels));

% h = 1, 10, 25
for ih = 1:3
    aum_temp(:,:,ih)      = squeeze(RDOS_m0(:,1,ih)) - MATRIX_mean;
    aum_temp_ci05(:,:,ih) = squeeze(RDOS_m0(:,2,ih)) - MATRIX_mean;
    aum_temp_ci95(:,:,ih) = squeeze(RDOS_m0(:,3,ih)) - MATRIX_mean;
end

% ------------------------------------------------------------------------
% 5. Long-run forecast to 2100
% ------------------------------------------------------------------------

FOR_comb    = NaN(length(Q),1);
CI_for_comb = NaN(length(Q),2);

for iq = 1:length(Q)

    selected_models = RDOS.GW.models(Q(iq),1:n_models); % Method 0

    BIC = BIC_all(Q(iq),selected_models==1);
    m = n_models - sum(isnan(selected_models));

    Weig_sbic = NaN(1,m);
    for k = 1:m
        Weig_sbic(k) = exp(-1/2*BIC(k)) / sum(exp(-1/2*BIC(1,:)));
    end

    Z = QUANTILES(:,Q(iq));

    [FORECAST,~,RES] = compute_forecast_models(Z,f_long);

    FORECAST2 = FORECAST(selected_models==1);
    RES = RES(:,selected_models==1);
    RES = rmmissing(RES);

    var_comb = compute_cov_comb(RES,Weig_sbic);

    FOR_comb(iq,1)    = Weig_sbic * FORECAST2;
    CI_for_comb(iq,:) = FOR_comb(iq) + norminv([0.05 0.95]) * sqrt(var_comb);

end

aum_temp(:,:,4)      = FOR_comb - MATRIX_mean;
aum_temp_ci05(:,:,4) = CI_for_comb(:,1) - MATRIX_mean;
aum_temp_ci95(:,:,4) = CI_for_comb(:,2) - MATRIX_mean;

AUM_temp.data      = aum_temp;
AUM_temp.ci05      = aum_temp_ci05;
AUM_temp.ci95      = aum_temp_ci95;
AUM_temp.structure = 'quantiles x periods x horizons';
AUM_temp.periods   = period_labels;
AUM_temp.quantiles = q_labels;
AUM_temp.H         = h_labels;

save(fullfile(RES_DIR,'AUM_temp_Globe_1880_2023_method0_w100.mat'),'AUM_temp');

% ------------------------------------------------------------------------
% 6. Figures: temperature increase by 2100
% ------------------------------------------------------------------------

fig1 = figure(1);
bar(aum_temp(:,:,4));
title('Temperature increase by 2100');
legend(period_labels,'Location','best');
set(gca,'XTick',1:length(Q));
set(gca,'XTickLabel',q_labels,'FontSize',12,'XTickLabelRotation',90);
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[20 16]);
set(gcf,'PaperPosition',[0 0 20 16]);
set(findall(gcf,'type','axes'),'FontSize',14);

print(fig1,fullfile(FIG_DIR,'Figure_Globe_1880_2023_f2100_method0_w100'),'-dpng');
print(fig1,fullfile(FIG_DIR,'Figure_Globe_1880_2023_f2100_method0_w100'),'-dpdf');
print(fig1,fullfile(FIG_DIR,'Figure_Globe_1880_2023_f2100_method0_w100'),'-deps');

fig2 = figure(2);
aum_temp2 = aum_temp;
aum_temp2(:,3,:) = [];
bar(aum_temp2(:,:,4));
title('Temperature increase by 2100');
legend(period_labels(1:2),'Location','best');
set(gca,'XTick',1:length(Q));
set(gca,'XTickLabel',q_labels,'FontSize',12,'XTickLabelRotation',90);
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[20 16]);
set(gcf,'PaperPosition',[0 0 20 16]);
set(findall(gcf,'type','axes'),'FontSize',14);

print(fig2,fullfile(FIG_DIR,'Figure_Globe_1880_2023_f2100_method0_w100_bis'),'-dpng');
print(fig2,fullfile(FIG_DIR,'Figure_Globe_1880_2023_f2100_method0_w100_bis'),'-dpdf');
print(fig2,fullfile(FIG_DIR,'Figure_Globe_1880_2023_f2100_method0_w100_bis'),'-deps');

% ------------------------------------------------------------------------
% 7. Mean temperature increases
% ------------------------------------------------------------------------

data_mean = [ ...
    mean(QUANTILES(1:21,1)), ...
    mean(QUANTILES(107:126,1)), ...
    mean(QUANTILES(116:135,1)) ];

q_mean = QUANTILES(:,1);
F_mean = [1,10,25,78];

FOR_comb_mean = NaN(length(F_mean),1);

% The original code uses the last selected model set and BIC weights from
% the long-run quantile loop. This is retained for exact replication.
for i = 1:length(F_mean)

    f = F_mean(i);
    [FORECAST,~,~] = compute_forecast_models(q_mean,f);
    FORECAST2 = FORECAST(selected_models==1);

    FOR_comb_mean(i) = Weig_sbic * FORECAST2;

end

Mean_increase = table( ...
    F_mean(:), ...
    FOR_comb_mean - data_mean(1), ...
    FOR_comb_mean - data_mean(2), ...
    FOR_comb_mean - data_mean(3), ...
    'VariableNames',{'Horizon','Increase_vs_1880_1900', ...
    'Increase_vs_1986_2002','Increase_vs_1995_2014'});

writetable(Mean_increase,fullfile(TAB_DIR,'Mean_temperature_increase_method0_w100.csv'));

disp('Increase of the mean');
disp(Mean_increase);
