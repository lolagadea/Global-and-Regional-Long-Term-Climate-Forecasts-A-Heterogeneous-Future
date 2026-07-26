% ========================================================================
% main_locate_quantiles_stations_Globe.m
% ------------------------------------------------------------------------
% Geographic location of stations contributing to q05, q50, and q95
% of the annual cross-sectional temperature distribution, 1880--2023.
%
% Output:
%   Figures/Figure_locate_quantiles_1880_2023.*
% ========================================================================

clear; clc;
warning('off');

% ------------------------------------------------------------------------
% 1. Portable paths
% ------------------------------------------------------------------------

THIS_DIR  = fileparts(mfilename('fullpath'));
GLOBE_DIR = fileparts(THIS_DIR);
ROOT_DIR  = fileparts(GLOBE_DIR);

FUN_DIR       = fullfile(ROOT_DIR,'functions');
LOCAL_FUN_DIR = fullfile(THIS_DIR,'functions');
DATA_DIR      = fullfile(ROOT_DIR,'01_Data_construction_CRUTEM');

FIGURES_DIR = fullfile(THIS_DIR,'Figures');

addpath(genpath(FUN_DIR));
addpath(genpath(LOCAL_FUN_DIR));

if ~isfolder(FIGURES_DIR), mkdir(FIGURES_DIR); end

% ------------------------------------------------------------------------
% 2. Load station data
% ------------------------------------------------------------------------

load(fullfile(DATA_DIR,'STATIONS.mat'));
load(fullfile(DATA_DIR,'index_crazy.dat'));

S = STATIONS.temp(:,:,2:13);
S(:,index_crazy,:) = [];

labels = [STATIONS.Latitude, STATIONS.Longitude, STATIONS.Height];

ID = str2num(cell2mat(STATIONS.ID)); %#ok<ST2NM>
ID(index_crazy,:) = [];

labels(index_crazy,:) = [];

years = STATIONS.years;

% ------------------------------------------------------------------------
% 3. Sample definition
% ------------------------------------------------------------------------

y1 = 31;    % 1880
y2 = 174;   % 2023
t  = y2-y1+1;

fprintf('\nLocating Globe quantile stations: %d--%d\n', years(y1), years(y2));

% ------------------------------------------------------------------------
% 4. Build Globe quantiles and selected station index
% ------------------------------------------------------------------------

index = NaN;

[~, num_s, A, INDEX] = ...
    function_build_quantiles_months_method1(S,years,y1,y2,index); %#ok<ASGLU>

fprintf('Selected month-station units, Globe: %d\n', num_s);

[~,n,~] = size(S);

B = reshape(INDEX,n,12);
selected_stations = nansum(B,2) > 0;

ID_selected = ID(selected_stations,:);
labels_sel  = labels(selected_stations,:);

pos = [ID_selected, labels_sel];

% Avoid negative heights in map/export routines
pos(pos(:,4)<0,4) = 0.001;

S2 = S(y1:y2,selected_stations,:);
Z  = nanmean(S2,3);

ns = size(pos,1);

fprintf('Effective stations, Globe: %d\n', ns);

% ------------------------------------------------------------------------
% 5. Locate q05 stations
% ------------------------------------------------------------------------

q05 = round(ns*5/100);
G   = NaN(t,q05,4);

for i = 1:t
    z = Z(i,:);
    [~,ix] = sort(z,'ascend');
    G(i,:,:) = pos(ix(1:q05),:);
end

COORD_q05 = [ ...
    reshape(G(:,:,1),t*q05,1), ...
    reshape(G(:,:,2),t*q05,1), ...
    reshape(G(:,:,3),t*q05,1), ...
    reshape(G(:,:,4),t*q05,1)];

[~,ixu] = unique(COORD_q05(:,1));
coord_q05 = COORD_q05(ixu,:);

% ------------------------------------------------------------------------
% 6. Locate q95 stations
% ------------------------------------------------------------------------

q95 = round(ns*5/100);
G   = NaN(t,q95,4);

for i = 1:t
    z = Z(i,:);
    [~,ix] = sort(z,'descend');
    G(i,:,:) = pos(ix(1:q95),:);
end

COORD_q95 = [ ...
    reshape(G(:,:,1),t*q95,1), ...
    reshape(G(:,:,2),t*q95,1), ...
    reshape(G(:,:,3),t*q95,1), ...
    reshape(G(:,:,4),t*q95,1)];

[~,ixu] = unique(COORD_q95(:,1));
coord_q95 = COORD_q95(ixu,:);

% ------------------------------------------------------------------------
% 7. Locate q50 stations
% ------------------------------------------------------------------------

med = median(Z,2);

st_median = [];

for i = 1:t
    for j = 1:ns
        if Z(i,j) < med(i)+2 && Z(i,j) > med(i)-2
            st_median = [st_median j]; %#ok<AGROW>
        end
    end
end

st_median = unique(st_median);
coord_q50 = pos(st_median,:);

% ------------------------------------------------------------------------
% 8. Combined figure
% ------------------------------------------------------------------------

figure('Color','w');

make_map_subplot_quantiles(coord_q05,coord_q50,coord_q95);

set(gcf,'PaperSize',[21.0 29.7], ...
        'PaperPosition',[0 0 21.0 29.7]);

fig_base = fullfile(FIGURES_DIR,'Figure_locate_quantiles_1880_2023');

print(gcf,fig_base,'-dpdf');
print(gcf,fig_base,'-dpng');
print(gcf,fig_base,'-deps');

fprintf('\nSaved figure in:\n%s\n', FIGURES_DIR);