% =========================================================================
% detect_crazy_stations.m
% -------------------------------------------------------------------------
% Replication package:
% Global and regional long-term climate forecasts: A heterogeneous future
%
% Purpose:
%   This script identifies CRUTEM5 stations with invalid or missing
%   geographical metadata.
%
% Input:
%   STATIONS.mat
%
% Output:
%   index_crazy.dat
%
% Notes:
%   - The script is fully portable: all paths are defined relative to the
%     location of this file.
%   - The original criteria for detecting problematic stations are preserved.
%   - index_crazy.dat contains the indices of stations with invalid longitude,
%     latitude, or height information.
%
% =========================================================================

clear; clc;

% -------------------------------------------------------------------------
% 1. Define portable paths
% -------------------------------------------------------------------------

THISFILE = mfilename('fullpath');
THISDIR  = fileparts(THISFILE);

INFILE  = fullfile(THISDIR, 'STATIONS.mat');
OUTFILE = fullfile(THISDIR, 'index_crazy.dat');

assert(isfile(INFILE), ...
    ['STATIONS.mat not found. Expected file: ', INFILE]);

% -------------------------------------------------------------------------
% 2. Load station database
% -------------------------------------------------------------------------

load(INFILE, 'STATIONS');

years  = STATIONS.years';
LAT    = STATIONS.Latitude;
LONG   = STATIONS.Longitude;
HEIGHT = STATIONS.Height;

n = length(LAT);
t = length(years);

% -------------------------------------------------------------------------
% 3. Detect stations with invalid geographical metadata
% -------------------------------------------------------------------------
% Original criteria:
%   - invalid longitude codes: -199.9000, -999.9000, 999.9000
%   - invalid latitude code  : -99.9000
%   - invalid height codes   : -999, -9999

index_long = find(LONG == -199.9000 | LONG == -999.9000 | LONG == 999.9000);
index_lat  = find(LAT  == -99.9000);
index_height = [find(HEIGHT == -999); find(HEIGHT == -9999)];

index = [index_long; index_lat; index_height];
index_crazy = unique(index);

% -------------------------------------------------------------------------
% 4. Save output
% -------------------------------------------------------------------------

save(OUTFILE, 'index_crazy', '-ASCII');

fprintf('\n============================================================\n');
fprintf('Problematic station metadata detected.\n');
fprintf('Output file: %s\n', OUTFILE);
fprintf('Number of stations checked: %d\n', n);
fprintf('Number of problematic stations: %d\n', length(index_crazy));
fprintf('============================================================\n');