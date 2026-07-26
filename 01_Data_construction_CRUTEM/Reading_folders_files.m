% =========================================================================
% Reading_folders_files.m
% -------------------------------------------------------------------------
% Replication package:
% Global and regional long-term climate forecasts: A heterogeneous future
%
% Purpose:
%   This script reads the raw CRUTEM5 station NetCDF files and constructs the
%   STATIONS structure used throughout the forecasting project.
%
% Input:
%   Raw CRUTEM5 station files stored in:
%
%       ./CRUTEM.5.0.2.0.station_files/station_files/
%
% Output:
%   STATIONS.mat
%
% Notes:
%   - The script is fully portable: all paths are defined relative to the
%     location of this file.
%   - The internal structure of STATIONS is kept unchanged relative to the
%     original replication codes.
%   - The sample is initialized from 1850 to 2024 because this is the time
%     span covered by the raw CRUTEM5 station files used in this revision.
%
% =========================================================================

clear; clc;

% -------------------------------------------------------------------------
% 1. Define portable paths
% -------------------------------------------------------------------------

THISFILE = mfilename('fullpath');
THISDIR  = fileparts(THISFILE);

RAW_DIR = fullfile(THISDIR, ...
    'CRUTEM.5.0.2.0.station_files', ...
    'station_files');

assert(isfolder(RAW_DIR), ...
    ['Raw CRUTEM station folder not found. Expected folder: ', RAW_DIR]);

fprintf('\n============================================================\n');
fprintf('Reading CRUTEM5 station files\n');
fprintf('============================================================\n');
fprintf('Script folder : %s\n', THISDIR);
fprintf('Raw data folder: %s\n\n', RAW_DIR);

% -------------------------------------------------------------------------
% 2. Find all folders containing station files
% -------------------------------------------------------------------------

folders  = cell(0);
nFolders = 0;

allFolders = dir(RAW_DIR);

for k = 1:numel(allFolders)

    if allFolders(k).isdir && ~ismember(allFolders(k).name, {'.','..'})
        nFolders = nFolders + 1;
        folders{nFolders,1} = fullfile(RAW_DIR, allFolders(k).name);
    end

end

fprintf('Number of subfolders found: %d\n', nFolders);

% -------------------------------------------------------------------------
% 3. Find all NetCDF station files inside these folders
% -------------------------------------------------------------------------

files        = cell(0);
nameStations = [];
nFiles       = 0;

for k = 1:nFolders

    allFiles = dir(folders{k});

    for l = 1:numel(allFiles)

        curFile = fullfile(folders{k}, allFiles(l).name);

        if ~allFiles(l).isdir
            nFiles = nFiles + 1;
            files{nFiles,1} = curFile;
            nameStations = [nameStations; string(allFiles(l).name)];
        end

    end

end

fprintf('Number of station files found: %d\n\n', nFiles);

assert(nFiles > 0, ...
    'No station files were found. Please check the raw CRUTEM folder.');

% -------------------------------------------------------------------------
% 4. Initialize STATIONS arrays
% -------------------------------------------------------------------------

START_YEAR = 1850;
END_YEAR_FULL = 2024;

T     = END_YEAR_FULL - START_YEAR + 1;
YEARS = START_YEAR:END_YEAR_FULL;

% S stores:
%   column 1  : year
%   columns 2-13 : monthly temperatures, January to December
S = NaN(T, nFiles, 13);

LAT    = NaN(nFiles,1);
LONG   = NaN(nFiles,1);
HEIGHT = NaN(nFiles,1);

ID              = cell(nFiles,1);
FIRST_GOOD_YEAR = cell(nFiles,1);
END_YEAR        = cell(nFiles,1);
NAME            = cell(nFiles,1);
COUNTRY         = cell(nFiles,1);

% -------------------------------------------------------------------------
% 5. Read information for each station
% -------------------------------------------------------------------------

for k = 1:nFiles

    fprintf('Reading station %d of %d\n', k, nFiles);

    % ---------------------------------------------------------------------
    % Read NetCDF variables
    % ---------------------------------------------------------------------
    temp      = ncread(files{k}, 'tas');
    time      = ncread(files{k}, 'time');
    latitude  = ncread(files{k}, 'latitude');
    longitude = ncread(files{k}, 'longitude');

    % ---------------------------------------------------------------------
    % Read station metadata
    % ---------------------------------------------------------------------
    first_good_year = ncreadatt(files{k}, '/', 'first_good_year');
    country         = ncreadatt(files{k}, '/', 'country');
    name            = ncreadatt(files{k}, '/', 'name');
    end_year        = ncreadatt(files{k}, '/', 'cru_data_to');
    height          = ncreadatt(files{k}, '/', 'height');
    id              = ncreadatt(files{k}, '/', 'stationid');

    % ---------------------------------------------------------------------
    % Time index
    % ---------------------------------------------------------------------
    % CRUTEM time is expressed relative to 1850. The original code uses
    % datevec(time) and then shifts years by 1850. We preserve this logic
    % to keep the generated STATIONS structure identical to the original one.
    dates  = datevec(time);
    years  = dates(:,1) + START_YEAR;
    months = dates(:,2);

    if months(1) ~= 1
        error('The first monthly observation is not January in file: %s', files{k});
    end

    % Position in the full 1850-2024 array
    pos_first_year = dates(1,1)   + 1;
    pos_end_year   = dates(end,1) + 1;

    t = length(time);

    % Each station file contains monthly data. Reshape into annual rows:
    % [year, Jan, Feb, ..., Dec]
    S(pos_first_year:pos_end_year, k, :) = ...
        [unique(years), reshape(temp, 12, t/12)'];

    % Store metadata
    LAT(k)    = latitude;
    LONG(k)   = longitude;
    HEIGHT(k) = height;

    ID{k}              = id;
    NAME{k}            = name;
    COUNTRY{k}         = country;
    FIRST_GOOD_YEAR{k} = first_good_year;
    END_YEAR{k}        = end_year;

end

% -------------------------------------------------------------------------
% 6. Build STATIONS structure
% -------------------------------------------------------------------------

STATIONS.temp = S;
STATIONS.years = YEARS;

STATIONS.Longitude = LONG;
STATIONS.Latitude  = LAT;
STATIONS.Height    = HEIGHT;

STATIONS.ID = ID;
STATIONS.first_good_year = FIRST_GOOD_YEAR;
STATIONS.end_year = END_YEAR;
STATIONS.name = NAME;
STATIONS.country = COUNTRY;

% -------------------------------------------------------------------------
% 7. Save output
% -------------------------------------------------------------------------

OUTFILE = fullfile(THISDIR, 'STATIONS.mat');

save(OUTFILE, 'STATIONS', '-v7.3');

fprintf('\n============================================================\n');
fprintf('STATIONS structure successfully created.\n');
fprintf('Output file: %s\n', OUTFILE);
fprintf('Number of station files processed: %d\n', nFiles);
fprintf('Years covered: %d-%d\n', START_YEAR, END_YEAR_FULL);
fprintf('============================================================\n');