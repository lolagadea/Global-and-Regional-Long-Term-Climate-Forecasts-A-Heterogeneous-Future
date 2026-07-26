% =========================================================================
% make_figure_density.m
% -------------------------------------------------------------------------
% Replication package:
% Global and regional long-term climate forecasts: A heterogeneous future
%
% Purpose:
%   Generate the 3D kernel-density figure for raw global temperature
%   month-station units over 1880-2023.
%
% Input:
%   STATIONS.mat        from 01_Data_construction_CRUTEM
%   index_crazy.dat     from 01_Data_construction_CRUTEM
%
% Output:
%   Figures/Fig_density_Globe_1880_2023.*
%
% Notes:
%   - This script reproduces Figure 3 of the paper.
%   - Densities are computed year by year over retained month-station units.
%   - The selected units must be observed throughout 1880-2023.
%
% =========================================================================

clear; clc;

warning('off','MATLAB:print:ContentTypeImageSuggested');

% -------------------------------------------------------------------------
% 1. Define portable paths
% -------------------------------------------------------------------------

THISFILE = mfilename('fullpath');
THISDIR  = fileparts(THISFILE);

ROOT = fileparts(fileparts(THISDIR));

DATA_DIR = fullfile(ROOT, '01_Data_construction_CRUTEM');
FIG_DIR  = fullfile(THISDIR, 'Figures');

if ~isfolder(FIG_DIR), mkdir(FIG_DIR); end

STATIONS_FILE = fullfile(DATA_DIR, 'STATIONS.mat');
CRAZY_FILE    = fullfile(DATA_DIR, 'index_crazy.dat');

assert(isfile(STATIONS_FILE), ...
    ['STATIONS.mat not found: ', STATIONS_FILE]);

assert(isfile(CRAZY_FILE), ...
    ['index_crazy.dat not found: ', CRAZY_FILE]);

load(STATIONS_FILE, 'STATIONS');
load(CRAZY_FILE, 'index_crazy');

% -------------------------------------------------------------------------
% 2. Prepare raw station-month data
% -------------------------------------------------------------------------

S = STATIONS.temp(:,:,2:13);     % years x stations x months
S(:,index_crazy,:) = [];

years = STATIONS.years(:);

y1 = find(years == 1880);
y2 = find(years == 2023);

S = S(y1:y2,:,:);
years_sample = years(y1:y2);

[T,N,M] = size(S);

% Convert station-month units into columns
A = reshape(S, T, N*M);

% Keep only station-month units observed for all years
ok = all(~isnan(A),1);
A  = A(:,ok);

fprintf('\n============================================================\n');
fprintf('Generating density figure: Globe 1880-2023\n');
fprintf('Selected month-station units: %d\n', size(A,2));
fprintf('============================================================\n');

% -------------------------------------------------------------------------
% 3. Compute yearly kernel densities
% -------------------------------------------------------------------------

% Temperature grid. This range reproduces the broad support used in the
% original figure and allows the kernel tails to be visible.
temp_grid = linspace(-75, 35, 100);

DENS = NaN(T, length(temp_grid));

for t = 1:T

    z = A(t,:);
    z = z(~isnan(z));

    % Remove physically implausible raw values, if any
    z(z < -60 | z > 60) = [];

    DENS(t,:) = ksdensity(z, temp_grid, ...
        'Kernel','epanechnikov');

end

% -------------------------------------------------------------------------
% 4. Plot 3D density surface
% -------------------------------------------------------------------------

[X,Y] = meshgrid(years_sample, temp_grid);

fig = figure('Color','w','Units','pixels','Position',[100 100 1000 720]);

surf(X, Y, DENS', ...
    'EdgeColor',[0.15 0.15 0.15], ...
    'LineWidth',0.15);

shading interp;
colormap(jet);

view([-35 28]);

grid on;
box on;

xlabel('years','FontSize',12);
ylabel('temperature in degrees Celsius (month-station units)','FontSize',12);
zlabel('density','FontSize',12);

set(gca,'FontSize',11);
set(gca,'YLim',[-75 35]);
set(gca,'ZLim',[0 0.07]);

% -------------------------------------------------------------------------
% 5. Export figure
% -------------------------------------------------------------------------

exportgraphics(fig, ...
    fullfile(FIG_DIR,'Fig_density_Globe_1880_2023.pdf'), ...
    'ContentType','vector', ...
    'BackgroundColor','white');

exportgraphics(fig, ...
    fullfile(FIG_DIR,'Fig_density_Globe_1880_2023.eps'), ...
    'ContentType','vector');

exportgraphics(fig, ...
    fullfile(FIG_DIR,'Fig_density_Globe_1880_2023.png'), ...
    'Resolution',600, ...
    'BackgroundColor','white');

fprintf('Density figure successfully created.\n');
fprintf('Figures folder: %s\n', FIG_DIR);
fprintf('============================================================\n');