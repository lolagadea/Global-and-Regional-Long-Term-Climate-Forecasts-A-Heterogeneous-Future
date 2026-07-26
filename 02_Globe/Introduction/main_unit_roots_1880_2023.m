% =========================================================================
% main_unit_roots_1880_2023.m
% -------------------------------------------------------------------------
% Replication package:
% Global and regional long-term climate forecasts: A heterogeneous future
%
% Purpose:
%   This script computes ADF unit root tests for the Globe distributional
%   characteristics over the 1880-2023 sample.
%
% Input:
%   Results/QUANTILES_monthly_Globe_1880_2023.mat
%
% Output:
%   Results/ADF_Globe_1880_2023.dat
%   Results/PVALUE_UR_Globe_1880_2023.dat
%   Results/LAGS_UR_Globe_1880_2023.dat
%   Tables/Table_UR_1880_2023.tex
%
% Notes:
%   - The econometric procedure is unchanged relative to the original code.
%   - ADF tests include intercept and trend.
%   - The lag length is selected by SBIC using calcula_adf_t_sbic.m.
%
% =========================================================================

clear; clc;

% -------------------------------------------------------------------------
% 1. Define portable paths
% -------------------------------------------------------------------------

THISFILE = mfilename('fullpath');
THISDIR  = fileparts(THISFILE);

% Root folder of the replication package: CODES_IJF
ROOT = fileparts(fileparts(THISDIR));

% General auxiliary functions
FUNC_DIR = fullfile(ROOT, 'functions');

% Local input/output folders
RES_DIR = fullfile(THISDIR, 'Results');
TAB_DIR = fullfile(THISDIR, 'Tables');

if ~isfolder(RES_DIR), mkdir(RES_DIR); end
if ~isfolder(TAB_DIR), mkdir(TAB_DIR); end

assert(isfolder(FUNC_DIR), ...
    ['General auxiliary functions folder not found: ', FUNC_DIR]);

addpath(genpath(FUNC_DIR));

% -------------------------------------------------------------------------
% 2. Load Globe quantiles
% -------------------------------------------------------------------------

INFILE = fullfile(RES_DIR, 'QUANTILES_monthly_Globe_1880_2023.mat');

assert(isfile(INFILE), ...
    ['Quantiles file not found: ', INFILE]);

load(INFILE, 'QUANTILES_monthly');

Y = QUANTILES_monthly.Globe;

[t,n] = size(Y);

% -------------------------------------------------------------------------
% 3. Compute ADF unit root tests
% -------------------------------------------------------------------------

ADF  = NaN(n,1);
P    = NaN(n,1);
LAG  = NaN(n,1);

pmax = round(12*(t/100)^0.25);

fprintf('\n============================================================\n');
fprintf('Computing ADF unit root tests: Globe 1880-2023\n');
fprintf('Number of characteristics: %d\n', n);
fprintf('Maximum lag length: %d\n', pmax);
fprintf('============================================================\n');

for i = 1:n

    y = Y(:,i);

    % ADF with intercept and trend; lag length selected by SBIC
    [adf,p] = calcula_adf_t_sbic(y,pmax);

    ADF(i,1) = adf;

    [~,pValue,~,~,~] = adftest(y, ...
        'model','TS', ...
        'lags',p);

    if pValue == 0.001
        pValue = 0;
    end

    P(i,1)   = pValue;
    LAG(i,1) = p;

end

% -------------------------------------------------------------------------
% 4. Save numerical results
% -------------------------------------------------------------------------

save(fullfile(RES_DIR,'ADF_Globe_1880_2023.dat'), ...
    'ADF','-ASCII');

save(fullfile(RES_DIR,'PVALUE_UR_Globe_1880_2023.dat'), ...
    'P','-ASCII');

save(fullfile(RES_DIR,'LAGS_UR_Globe_1880_2023.dat'), ...
    'LAG','-ASCII');

% -------------------------------------------------------------------------
% 5. Create LaTeX table
% -------------------------------------------------------------------------

rowLabels    = name_labels;
columnLabels = ['Characteristic&','ADF-SBIC&', 'p-value&', 'lags&'];

file = fullfile(TAB_DIR, 'Table_UR_1880_2023.tex');

fid = fopen(file,'w');

k = 3;

txt = ['\\begin{table}\\caption{Unit root tests}\\label{UR}' ...
       '\\begin{center}\\scalebox{1}{\\begin{tabular}{l' ...
       '*{' int2str(k+1) '}{c}} \\hline \\hline \n'];

fprintf(fid,txt);

txt = [columnLabels, '\\\\ \\hline \n'];
fprintf(fid,txt);

formatSpec = '%s %5.4f %5.3f %5.0f';

for i = 1:length(ADF)

    txt = [rowLabels{i} '&' ...
           num2str(ADF(i), '%5.2f') '&' ...
           num2str(P(i),   '%5.3f') '&' ...
           num2str(LAG(i), '%5.0f')];

    fprintf(fid,formatSpec,txt);
    fprintf(fid, '\\\\ \n');

end

txt = '\\hline \\hline \\end{tabular}}\\end{center}\\end{table}\n';

fprintf(fid,txt);
fclose(fid);

fprintf('\n============================================================\n');
fprintf('Unit root tests successfully computed.\n');
fprintf('Results folder: %s\n', RES_DIR);
fprintf('Table file    : %s\n', file);
fprintf('============================================================\n');