% =========================================================================
% main_quantiles_Globe_monthly.m
% -------------------------------------------------------------------------
% Replication package:
% Global and regional long-term climate forecasts: A heterogeneous future
%
% Purpose:
%   This script constructs the monthly unconditional distributional
%   characteristics for the Globe using CRUTEM5 station-month units.
%
% Input:
%   STATIONS.mat        from 01_Data_construction_CRUTEM
%   index_crazy.dat     from 01_Data_construction_CRUTEM
%
% Output:
%   Results/QUANTILES_monthly_Globe_1880_2023.mat
%   Results/QUANTILES_monthly_Globe_1960_2023.mat
%   Tables/Selected_stations_Globe_1880_2023.xlsx
%   Tables/Selected_stations_Globe_1960_2023.xlsx
%   Figures/Figure_stations_Globe_1880_2023.*
%   Figures/Figure_stations_Globe_1960_2023.*
%
% Notes:
%   - The script is fully portable: all paths are defined relative to this
%     file and to the root folder of the replication package.
%   - The original sample selection and quantile construction are preserved.
%   - The script runs two Globe samples:
%       (i)  1880-2023, used for the long historical Globe analysis.
%       (ii) 1960-2023, used for comparison with the regional analysis.
%   - For each sample, the effective set of stations is selected separately.
%
% =========================================================================

clear; clc;
warning('off','MATLAB:xlswrite:AddSheet');
warning('off','MATLAB:print:ContentTypeImageSuggested');

% -------------------------------------------------------------------------
% 1. Define portable paths
% -------------------------------------------------------------------------

THISFILE = mfilename('fullpath');
THISDIR  = fileparts(THISFILE);

% Root folder of the replication package: CODES_IJF
ROOT = fileparts(fileparts(THISDIR));

% -------------------------------------------------------------------------
% Main folders
% -------------------------------------------------------------------------

DATA_DIR = fullfile(ROOT, '01_Data_construction_CRUTEM');

% General auxiliary functions
FUNC_DIR = fullfile(ROOT, 'functions');

% Auxiliary functions specific to Globe/Introduction
FUNC_INTRO_DIR = fullfile(THISDIR, 'functions');

% Output folders
FIG_DIR = fullfile(THISDIR, 'Figures');
TAB_DIR = fullfile(THISDIR, 'Tables');
RES_DIR = fullfile(THISDIR, 'Results');

% -------------------------------------------------------------------------
% Create output folders if needed
% -------------------------------------------------------------------------

if ~isfolder(FIG_DIR), mkdir(FIG_DIR); end
if ~isfolder(TAB_DIR), mkdir(TAB_DIR); end
if ~isfolder(RES_DIR), mkdir(RES_DIR); end

% -------------------------------------------------------------------------
% Check required folders
% -------------------------------------------------------------------------

assert(isfolder(DATA_DIR), ...
    ['Data folder not found: ', DATA_DIR]);

assert(isfolder(FUNC_DIR), ...
    ['General auxiliary functions folder not found: ', FUNC_DIR]);

assert(isfolder(FUNC_INTRO_DIR), ...
    ['Introduction auxiliary functions folder not found: ', FUNC_INTRO_DIR]);

% -------------------------------------------------------------------------
% Add auxiliary functions to MATLAB path
% -------------------------------------------------------------------------

addpath(genpath(FUNC_DIR));
addpath(genpath(FUNC_INTRO_DIR));
% -------------------------------------------------------------------------
% 2. Load data
% -------------------------------------------------------------------------

STATIONS_FILE = fullfile(DATA_DIR, 'STATIONS.mat');
CRAZY_FILE    = fullfile(DATA_DIR, 'index_crazy.dat');

assert(isfile(STATIONS_FILE), ...
    ['STATIONS.mat not found: ', STATIONS_FILE]);

assert(isfile(CRAZY_FILE), ...
    ['index_crazy.dat not found: ', CRAZY_FILE]);

load(STATIONS_FILE, 'STATIONS');
load(CRAZY_FILE, 'index_crazy');

% -------------------------------------------------------------------------
% 3. Prepare station-month data
% -------------------------------------------------------------------------

S = STATIONS.temp(:,:,2:13);
S(:,index_crazy,:) = [];

labels = [STATIONS.Latitude, STATIONS.Longitude, STATIONS.Height];

ID = str2num(cell2mat(STATIONS.ID)); %#ok<ST2NM>
ID(index_crazy,:) = [];

labels(index_crazy,:) = [];

years = STATIONS.years;

name_station = STATIONS.name;
country      = STATIONS.country;

[t,n,m] = size(S); 

names = name_labels; 
q = 19;

% -------------------------------------------------------------------------
% 4. Define samples to be constructed
% -------------------------------------------------------------------------
% Original indexing:
%   y1 = 31   -> starting year 1880
%   y1 = 101  -> starting year 1950
%   y1 = 111  -> starting year 1960
%   y1 = 121  -> starting year 1970
%
%   y2 = 174  -> ending year 2023
%   y2 = 173  -> ending year 2022

samples = {
    '1880_2023', 31,  174
    '1960_2023', 111, 174
    '1880_1960', 31,  111
};

% -------------------------------------------------------------------------
% 5. Build monthly quantiles for each Globe sample
% -------------------------------------------------------------------------

for s = 1:size(samples,1)

    sample_name = samples{s,1}; %#ok<NASGU>
    y1 = samples{s,2};
    y2 = samples{s,3};

    fprintf('\n============================================================\n');
    fprintf('Constructing Globe monthly quantiles\n');
    fprintf('Sample: %d-%d\n', years(y1), years(y2));
    fprintf('============================================================\n');

    index = NaN;

    [QUANTILES_monthly_Globe, num_s, A, INDEX] = ...
        function_build_quantiles_months_method1(S, years, y1, y2, index); %#ok<ASGLU>

    fprintf('Number of selected month-station units for the Globe: %d\n', num_s);

    % ---------------------------------------------------------------------
    % Identify effective stations for the current sample
    % ---------------------------------------------------------------------

    B = reshape(INDEX, n, 12);

    selected_stations = zeros(n,1);

    for i = 1:n
        if nansum(B(i,:)) > 0
            selected_stations(i) = 1;
        end
    end

    SE = [ID, labels];

    ID_selected = ID(selected_stations == 1,:);
    SE2 = [ID_selected, labels(selected_stations == 1,:)];

    name_station2 = name_station(selected_stations == 1,:);
    country2      = country(selected_stations == 1,:);

    fprintf('Number of effective stations in the Globe: %d\n', size(SE2,1));

    % ---------------------------------------------------------------------
    % Save selected-station information to Excel
    % ---------------------------------------------------------------------

    name_file_excel = fullfile(TAB_DIR, ...
        strcat('Selected_stations_Globe_', ...
        num2str(years(y1)), '_', num2str(years(y2)), '.xlsx'));

    xlswrite(name_file_excel, ...
        {'station','latitude','longitude','height','name-station','country'}, ...
        'all','a1:f1');

    xlswrite(name_file_excel, SE, 'all','a2');
    xlswrite(name_file_excel, name_station, 'all','e2');
    xlswrite(name_file_excel, country, 'all','f2');

    xlswrite(name_file_excel, ...
        {'station','latitude','longitude','height','name-station','country'}, ...
        'selected','a1:f1');

    xlswrite(name_file_excel, SE2, 'selected','a2');
    xlswrite(name_file_excel, name_station2, 'selected','e2');
    xlswrite(name_file_excel, country2, 'selected','f2');

    % ---------------------------------------------------------------------
    % Make station map
    % ---------------------------------------------------------------------

    make_map_function_Globe(SE2);

    fig_base = fullfile(FIG_DIR, ...
        strcat('Figure_stations_Globe_', ...
        num2str(years(y1)), '_', num2str(years(y2))));

    exportgraphics(gcf, strcat(fig_base,'.png'), ...
        'Resolution',300, ...
        'BackgroundColor',[0.729411780834198 0.831372559070587 0.95686274766922]);

    exportgraphics(gcf, strcat(fig_base,'.pdf'), ...
        'ContentType','vector', ...
        'BackgroundColor',[0.729411780834198 0.831372559070587 0.95686274766922]);

    exportgraphics(gcf, strcat(fig_base,'.eps'), ...
        'ContentType','vector', ...
        'BackgroundColor',[0.729411780834198 0.831372559070587 0.95686274766922]);

    close(gcf);

    % ---------------------------------------------------------------------
    % Save quantile structure
    % ---------------------------------------------------------------------

    QUANTILES_monthly = [];
    QUANTILES_monthly.Globe = QUANTILES_monthly_Globe;
    QUANTILES_monthly.years = years(y1:y2);

    name_file = fullfile(RES_DIR, ...
        strcat('QUANTILES_monthly_Globe_', ...
        num2str(years(y1)), '_', num2str(years(y2)), '.mat'));

    save(name_file, 'QUANTILES_monthly');

    fprintf('Quantiles file saved: %s\n', name_file);
    fprintf('Excel file saved    : %s\n', name_file_excel);

end

fprintf('\n============================================================\n');
fprintf('All Globe monthly quantile files successfully created.\n');
fprintf('Results folder: %s\n', RES_DIR);
fprintf('Tables folder : %s\n', TAB_DIR);
fprintf('Figures folder: %s\n', FIG_DIR);
fprintf('============================================================\n');