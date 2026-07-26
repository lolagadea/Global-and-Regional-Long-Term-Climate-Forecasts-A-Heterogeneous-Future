% ========================================================================
% run_forecast_Africa.m
% ------------------------------------------------------------------------
% Forecast competition for Africa temperature distribution characteristics.
% Replication package IJF - Regional Forecast module.
% ========================================================================

clear; clc;
clear global;
warning('off');

restoredefaultpath;
rehash toolboxcache;

% ------------------------------------------------------------------------
% 1. Define portable paths
% ------------------------------------------------------------------------

THIS_DIR     = fileparts(mfilename('fullpath'));   % .../03_Regional/Forecast/Africa
FORECAST_DIR = fileparts(THIS_DIR);                % .../03_Regional/Forecast
REGIONAL_DIR = fileparts(FORECAST_DIR);            % .../03_Regional
ROOT_DIR     = fileparts(REGIONAL_DIR);            % .../CODES_IJF

FUN_GENERAL_DIR  = fullfile(ROOT_DIR,'functions');
FUN_FORECAST_DIR = fullfile(ROOT_DIR,'functions_for');

FIG_DIR = fullfile(THIS_DIR,'Figures');
RES_DIR = fullfile(THIS_DIR,'Results_for');
TAB_DIR = fullfile(THIS_DIR,'Tables');

if ~isfolder(FIG_DIR), mkdir(FIG_DIR); end
if ~isfolder(RES_DIR), mkdir(RES_DIR); end
if ~isfolder(TAB_DIR), mkdir(TAB_DIR); end

addpath(genpath(FUN_GENERAL_DIR));
addpath(genpath(FUN_FORECAST_DIR));

% ------------------------------------------------------------------------
% 2. Load regional quantiles
% ------------------------------------------------------------------------

DATA_FILE = fullfile(REGIONAL_DIR,'Introduction','Results', ...
    'QUANTILES_monthly_1960_2023.mat');

assert(isfile(DATA_FILE), ...
    'Input file QUANTILES_monthly_1960_2023.mat not found.');

load(DATA_FILE,'QUANTILES_monthly','name_labels');

Y = QUANTILES_monthly.Africa;
[t,n] = size(Y);

years = 1960:2023;
names = name_labels;

% ------------------------------------------------------------------------
% 3. Forecast settings
% ------------------------------------------------------------------------

h = [1, 10, 25];   % keep 50 if GW tests require for50
w = [25, 50];

models = { ...
    'mean', ...
    'linear-trend', ...
    'pol-trend', ...
    'pol-trend-av-sl', ...
    'pol-trend-log', ...
    'struct-breaks', ...
    'pol-trend-arp', ...
    'pol-trend-arp-av-sl', ...
    'arp', ...
    'rw', ...
    'rwd', ...
    'ima', ...
    'arfima', ...
    'arp20', ...
    'combined', ...
    'combined*', ...
    'combined**'};

m = length(models);

% ------------------------------------------------------------------------
% 4. Initialize output structure
% ------------------------------------------------------------------------

RDOS = struct();

RDOS.Xf = NaN(t,n,m,length(w),length(h));

RDOS.param.k  = NaN(n,t,m,length(w),length(h));
RDOS.param.p  = NaN(n,t,m,length(w),length(h));
RDOS.param.sb = NaN(n,t,length(w),length(h));

RDOS.sbic = NaN(n,t,m,length(w),length(h));

RDOS.region = 'Africa';
RDOS.models = models;
RDOS.names  = names;
RDOS.years  = years;
RDOS.W      = w;
RDOS.H      = h;

RDOS.readme = { ...
    'RDOS.Xf dimensions are: time x characteristics x models x window x horizon'; ...
    'RDOS.sbic dimensions are: characteristics x time x models x window x horizon'; ...
    'RDOS.param.k stores selected polynomial degrees'; ...
    'RDOS.param.p stores selected AR orders'; ...
    'RDOS.param.sb stores selected structural-break years'; ...
    'Forecast horizons: 1, 10, 25, 50'; ...
    'Rolling windows: 25, 50'; ...
    'Models 15-17 are placeholders for combined forecasts constructed downstream'};

% ------------------------------------------------------------------------
% 5. Rolling forecast competition
% ------------------------------------------------------------------------
% Some auxiliary functions save figures or .dat files without an explicit
% path. To keep the replication package organized, we temporarily move to
% Figures while preserving the original working directory.

ORIG_DIR = pwd;
cd(FIG_DIR);
cleanupObj = onCleanup(@() cd(ORIG_DIR));

for j = 1:length(w)

    for i = 1:length(h)

        fprintf('\nRunning Africa: window = %d, horizon = %d\n', w(j), h(i));

        if w(j) + h(i) >= t
            warning('Sample size is smaller than w+h: t=%d, w=%d, h=%d', ...
                t, w(j), h(i));

            RDOS.warning = sprintf( ...
                'Sample size is smaller than w+h: t=%d, w=%d, h=%d', ...
                t, w(j), h(i));
            continue
        end

        idx_store = w(j) + h(i):t;

        [Xf, sbic] = mean_roll_forecast_h(Y,years,names,h(i),w(j));
        RDOS.Xf(:,:,1,j,i) = Xf;
        RDOS.sbic(:,idx_store,1,j,i) = sbic;

        [Xf, sbic] = linear_trend_roll_forecast_h(Y,years,names,h(i),w(j));
        RDOS.Xf(:,:,2,j,i) = Xf;
        RDOS.sbic(:,idx_store,2,j,i) = sbic;

        [Xf,K, sbic] = pol_trend_roll_forecast_h(Y,years,names,h(i),w(j));
        RDOS.Xf(:,:,3,j,i) = Xf;
        RDOS.param.k(:,idx_store,3,j,i) = K;
        RDOS.sbic(:,idx_store,3,j,i) = sbic;

        [Xf,K, sbic] = pol_trend_average_slope_roll_forecast_h(Y,years,names,h(i),w(j));
        RDOS.Xf(:,:,4,j,i) = Xf;
        RDOS.param.k(:,idx_store,4,j,i) = K;
        RDOS.sbic(:,idx_store,4,j,i) = sbic;

        [Xf,K, sbic] = pol_trend_log_roll_forecast_h(Y,years,names,h(i),w(j));
        RDOS.Xf(:,:,5,j,i) = Xf;
        RDOS.param.k(:,idx_store,5,j,i) = K;
        RDOS.sbic(:,idx_store,5,j,i) = sbic;

        [Xf,TB,YEARS_breaks,sbic] = struct_break_roll_forecast_h(Y,years,names,h(i),w(j)); %#ok<ASGLU>
        RDOS.Xf(:,:,6,j,i) = Xf;
        RDOS.param.sb(:,idx_store,j,i) = YEARS_breaks;
        RDOS.sbic(:,idx_store,6,j,i) = sbic;

        [Xf,K,P,sbic] = pol_trend_arp_roll_forecast_h(Y,years,names,h(i),w(j));
        RDOS.Xf(:,:,7,j,i) = Xf;
        RDOS.param.k(:,idx_store,7,j,i) = K;
        RDOS.param.p(:,idx_store,7,j,i) = P;
        RDOS.sbic(:,idx_store,7,j,i) = sbic;

        [Xf,K,P,sbic] = pol_trend_arp_average_slope_roll_forecast_h(Y,years,names,h(i),w(j));
        RDOS.Xf(:,:,8,j,i) = Xf;
        RDOS.param.k(:,idx_store,8,j,i) = K;
        RDOS.param.p(:,idx_store,8,j,i) = P;
        RDOS.sbic(:,idx_store,8,j,i) = sbic;

        [Xf,P,sbic] = arp_roll_forecast_h(Y,years,names,h(i),w(j));
        RDOS.Xf(:,:,9,j,i) = Xf;
        RDOS.param.p(:,idx_store,9,j,i) = P;
        RDOS.sbic(:,idx_store,9,j,i) = sbic;

        [Xf,sbic] = rw_roll_forecast_h(Y,years,names,h(i),w(j));
        RDOS.Xf(:,:,10,j,i) = Xf;
        RDOS.sbic(:,idx_store,10,j,i) = sbic;

        [Xf,sbic] = rw_drift_roll_forecast_h(Y,years,names,h(i),w(j));
        RDOS.Xf(:,:,11,j,i) = Xf;
        RDOS.sbic(:,idx_store,11,j,i) = sbic;

        [Xf,sbic] = ima_roll_forecast_h(Y,years,names,h(i),w(j));
        RDOS.Xf(:,:,12,j,i) = Xf;
        RDOS.sbic(:,idx_store,12,j,i) = sbic;

        [Xf,sbic] = arfima_roll_forecast_h(Y,years,names,h(i),w(j));
        RDOS.Xf(:,:,13,j,i) = Xf;
        RDOS.sbic(:,idx_store,13,j,i) = sbic;

        [Xf,sbic] = arp20_roll_forecast_h(Y,years,names,h(i),w(j));
        RDOS.Xf(:,:,14,j,i) = Xf;
        RDOS.sbic(:,idx_store,14,j,i) = sbic;

    end

end

% ------------------------------------------------------------------------
% 6. Save results
% ------------------------------------------------------------------------

OUT_FILE = fullfile(RES_DIR,'RDOS_Africa_1960_2023.mat');
save(OUT_FILE,'RDOS');

% Move auxiliary .dat outputs generated by legacy functions to Results_for
datFiles = dir(fullfile(FIG_DIR,'*.dat'));
for ii = 1:numel(datFiles)
    movefile(fullfile(FIG_DIR,datFiles(ii).name), ...
             fullfile(RES_DIR,datFiles(ii).name));
end

fprintf('\nForecast competition completed for Africa.\n');
fprintf('Results saved in:\n%s\n',OUT_FILE);

clear global