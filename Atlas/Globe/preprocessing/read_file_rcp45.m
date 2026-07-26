% ========================================================================
% read_file_rcp45.m
% ------------------------------------------------------------------------
% Reads NetCDF data and constructs global temperature series.
% Portable version: input NetCDF files are read from ../data_nc and
% TEMP_*.mat outputs are saved into ../results.
% ========================================================================

clear; clc;
warning('off');

% ------------------------------------------------------------------------
% 1. Portable paths
% ------------------------------------------------------------------------
THIS_DIR = fileparts(mfilename('fullpath'));

% Script location:
% CODES_IJF/Atlas/Globe/preprocessing
ATLAS_DIR = fileparts(THIS_DIR);

DATA_DIR    = fullfile(ATLAS_DIR, 'data_nc');
RESULTS_DIR = fullfile(ATLAS_DIR, 'results');

if ~exist(DATA_DIR, 'dir')
    error('Data folder not found: %s', DATA_DIR);
end

if ~exist(RESULTS_DIR, 'dir')
    mkdir(RESULTS_DIR);
end

% ------------------------------------------------------------------------
% 2. Input and output files
% ------------------------------------------------------------------------
NC_FILE = fullfile(DATA_DIR, 't_CMIP5_rcp45_mon_200601-210012.nc');
OUT_FILE = fullfile(RESULTS_DIR, 'TEMP_2006_2100_CMIP5_Globe_rcp45.mat');

assert(exist(NC_FILE, 'file') == 2, ...
    'NetCDF file not found: %s', NC_FILE);

% ------------------------------------------------------------------------
% 3. Original processing code
% ------------------------------------------------------------------------
ncdisp(NC_FILE)
gcm_model  = ncread(NC_FILE,'gcm_model');
gcm_variant  = ncread(NC_FILE,'gcm_variant');
t  = ncread(NC_FILE,'t');
time  = ncread(NC_FILE,'time');
lat= ncread(NC_FILE,'lat');
lon= ncread(NC_FILE,'lon');
m = ncreadatt(NC_FILE,'t','missing_value');
%Adjust the time
dt= datetime(time, 'ConvertFrom', 'datenum');
dt2=datevec(dt);
years=dt2(:,1)+1850; months=dt2(:,2);
years=[years(2:end);2100]; months=[months(2:end);12];

[nlon,nlat,nt,nmodels]=size(t);
N=nlon*nlat*nt*nmodels;

for i=1:nlon
    for j=1:nlat
        for k=1:nt
            for v=1:nmodels
    if t(i,j,k,v)==m
        t(i,j,k,v)=NaN;
    end
            end
        end
    end
end

w=cos(lat*pi/180);
w=w/sum(w);

tw=NaN(nlon,nt,nmodels);
for k=1:nt
    for v=1:nmodels
        for i=1:nlon
            a=squeeze(t(i,:,k,v));
            aw=a*w;
            tw(i,k,v)=aw;
        end
    end
end


for k=1:nt
    b=reshape(tw(:,k,:),nlon*nmodels,1);
    temp_month(k,1)=nanmean(b);
end

years2=unique(years);
nyears=length(years2);
temp_years=NaN(nyears,1);
for i=1:nyears
fprintf('Year block %d / %d\n', i, nyears);
    temp_years(i,:)=mean(temp_month((i-1)*12+1:i*12,:));
end

TEMP.data.montly=temp_month;
TEMP.data.yearly=temp_years;
TEMP.years.montly=years;
TEMP.years.years=years2;
TEMP.months=months;
TEMP.models=gcm_model;
TEMP.variant=gcm_variant;
TEMP.lat=lat;
TEMP.lon=lon;

save(OUT_FILE, 'TEMP', '-v7.3');




fprintf('Saved:\n%s\n', OUT_FILE);
