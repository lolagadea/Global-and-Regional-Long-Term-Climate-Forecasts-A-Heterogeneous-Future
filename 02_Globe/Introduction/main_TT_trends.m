% =========================================================================
% main_TT_trends.m
% -------------------------------------------------------------------------
% Replication package:
% Global and regional long-term climate forecasts: A heterogeneous future
%
% Purpose:
%   Estimate linear trends for Globe distributional characteristics using
%   HAC standard errors.
%
% Input:
%   Results/QUANTILES_monthly_Globe_1880_2023.mat
%
% Output:
%   Results/BETA_Globe_1880_2023.dat
%   Results/PVALUE_Globe_1880_2023.dat
%   Tables/Table_TT_cross_1880_2023.tex
%   Figures/Figure4_joint_Globe_1880_2023.*
%   Figures/Figure4a_characteristics_Globe_1880_2023.*
%   Figures/Figure4b_quantiles_Globe_1880_2023.*
%
% Notes:
%   - The econometric calculation is unchanged relative to the original code.
%   - Trends are estimated from Ct = alpha + beta*t + u_t.
%   - HAC t-ratios are computed using ols_hac.m.
%
% =========================================================================

clear; clc;

% -------------------------------------------------------------------------
% Suppress non-critical export warnings
% -------------------------------------------------------------------------

warning('off','MATLAB:print:ContentTypeImageSuggested');

% -------------------------------------------------------------------------
% 1. Define portable paths
% -------------------------------------------------------------------------

THISFILE = mfilename('fullpath');
THISDIR  = fileparts(THISFILE);

% Root folder of the replication package: CODES_IJF
ROOT = fileparts(fileparts(THISDIR));

FUNC_DIR = fullfile(ROOT, 'functions');

RES_DIR = fullfile(THISDIR, 'Results');
TAB_DIR = fullfile(THISDIR, 'Tables');
FIG_DIR = fullfile(THISDIR, 'Figures');

if ~isfolder(RES_DIR), mkdir(RES_DIR); end
if ~isfolder(TAB_DIR), mkdir(TAB_DIR); end
if ~isfolder(FIG_DIR), mkdir(FIG_DIR); end

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
years = QUANTILES_monthly.years;
years = years(:);

[t,q] = size(Y);

fprintf('\n============================================================\n');
fprintf('Estimating trend tests: Globe 1880-2023\n');
fprintf('Number of characteristics: %d\n', q);
fprintf('Sample: %d-%d\n', years(1), years(end));
fprintf('============================================================\n');

% -------------------------------------------------------------------------
% 3. Estimate linear trends with HAC standard errors
% -------------------------------------------------------------------------

X = [ones(t,1), (1:t)'];

BETA   = NaN(q,1);
TRATIO = NaN(q,1);
PVALUE = NaN(q,1);

for i = 1:q

    y = Y(:,i);

    [b,t_ratio] = ols_hac(y,X);

    BETA(i)   = b(2);
    TRATIO(i) = t_ratio(2);
    PVALUE(i) = 2*(1 - tcdf(abs(TRATIO(i)), t-2));

end

% -------------------------------------------------------------------------
% 4. Save numerical results
% -------------------------------------------------------------------------

save(fullfile(RES_DIR,'BETA_Globe_1880_2023.dat'), ...
    'BETA','-ASCII');

save(fullfile(RES_DIR,'PVALUE_Globe_1880_2023.dat'), ...
    'PVALUE','-ASCII');

save(fullfile(RES_DIR,'TRATIO_Globe_1880_2023.dat'), ...
    'TRATIO','-ASCII');

% -------------------------------------------------------------------------
% 5. Create LaTeX table
% -------------------------------------------------------------------------

rowLabels = name_labels;
columnLabels = ['Characteristic&','Coeff&','p-value'];

file = fullfile(TAB_DIR, 'Table_TT_cross_1880_2023.tex');

fid = fopen(file,'w');

txt = ['\\begin{table}\\caption{Trend test (CRU data, 1880--2023)}' ...
       '\\label{Tab-TT-cross-1880-2023}\\begin{center}' ...
       '\\scalebox{1}{\\begin{tabular}{l*{2}{c}} \\hline \\hline \n'];

fprintf(fid,txt);

txt = [columnLabels, '\\\\ \\hline \n'];
fprintf(fid,txt);

for i = 1:length(BETA)

    txt = [rowLabels{i} '&' ...
           num2str(BETA(i),   '%5.4f') '&' ...
           num2str(PVALUE(i), '%5.4f')];

    fprintf(fid,'%s',txt);
    fprintf(fid,'\\\\ \n');

end

txt = '\\hline \\hline \\end{tabular}}\\end{center}\\end{table}\n';

fprintf(fid,txt);
fclose(fid);

% -------------------------------------------------------------------------
% 6. Figure 4A: distributional characteristics
% -------------------------------------------------------------------------

figW = 7.0;
figH = 5.0;

f4a = figure('Color','w','Units','inches','Position',[1 1 figW figH]);

tl4a = tiledlayout(f4a,4,2,'TileSpacing','compact','Padding','tight');

names = {'mean','max','min','std','iqr','range','kur','skw'};
time  = years(:);
T     = numel(time);
tickx = round(linspace(1,T,8));

for i = 1:8

    ax = nexttile(tl4a);

    plot(ax, time, Y(:,i), 'LineWidth',1.1);
    grid(ax,'on');
    box(ax,'on');

    title(ax, upper(names{i}), ...
        'FontWeight','normal', ...
        'FontSize',9);

    set(ax,'FontSize',8);

    xticks(ax, time(tickx));
    xtickangle(ax,0);

    if mod(i,2) == 1
        ylabel(ax,'Temp (°C)','FontSize',8);
    else
        ax.YTickLabel = [];
    end

end

xlabel(tl4a,'Year','FontSize',9);

title(tl4a, ...
    'Distributional characteristics of global temperature (CRU 1880–2023)', ...
    'FontSize',11);

set(f4a,'ToolBar','none','MenuBar','none');

axall = findall(f4a,'Type','axes');

for k = 1:numel(axall)
    try
        axtoolbar(axall(k),'Visible','off');
    catch
    end
    if exist('disableDefaultInteractivity','file')
        disableDefaultInteractivity(axall(k));
    end
end

exportgraphics(f4a, ...
    fullfile(FIG_DIR,'Figure4a_characteristics_Globe_1880_2023.pdf'), ...
    'ContentType','vector','BackgroundColor','white');

exportgraphics(f4a, ...
    fullfile(FIG_DIR,'Figure4a_characteristics_Globe_1880_2023.eps'), ...
    'ContentType','vector');

exportgraphics(f4a, ...
    fullfile(FIG_DIR,'Figure4a_characteristics_Globe_1880_2023.png'), ...
    'Resolution',600,'BackgroundColor','white');

% -------------------------------------------------------------------------
% 7. Figure 4B: quantiles
% -------------------------------------------------------------------------

figW = 7.0;
figH = 3.8;

f4b = figure('Color','w','Units','inches','Position',[1 1 figW figH]);

tl4b = tiledlayout(f4b,1,1,'TileSpacing','compact','Padding','tight');

ax = nexttile(tl4b);
hold(ax,'on');
grid(ax,'on');
box(ax,'on');
set(ax,'FontSize',9);

allQ = Y(:,9:19);

% Background lines: q10-q90
h_bg = plot(ax, time, allQ(:,2:10), 'LineWidth',0.8);
set(h_bg, 'HandleVisibility','off');

% Highlight selected quantiles
h05 = plot(ax, time, allQ(:,1),  'LineWidth',1.6, 'LineStyle','--');
h50 = plot(ax, time, allQ(:,6),  'LineWidth',1.8, 'LineStyle','-');
h95 = plot(ax, time, allQ(:,11), 'LineWidth',1.6, 'LineStyle',':');

% Proxy handle for context quantiles
h_ctx = plot(ax, NaN, NaN, 'LineWidth',0.9);

lg = legend(ax, [h_ctx h05 h50 h95], ...
    {'q10..q90 (context)','q05','q50','q95'}, ...
    'Location','southoutside', ...
    'Orientation','horizontal', ...
    'Box','off', ...
    'FontSize',9);

lg.NumColumns = 4;
lg.Layout.Tile = 'south';

xticks(ax, time(tickx));
xtickangle(ax,0);

xlabel(ax,'Year','FontSize',10);
ylabel(ax,'Temp (°C)','FontSize',10);

title(tl4b, ...
    'Quantiles of global temperature (CRU 1880–2023)', ...
    'FontSize',11);

set(f4b,'ToolBar','none','MenuBar','none');

exportgraphics(f4b, ...
    fullfile(FIG_DIR,'Figure4b_quantiles_Globe_1880_2023.pdf'), ...
    'ContentType','vector','BackgroundColor','white');

exportgraphics(f4b, ...
    fullfile(FIG_DIR,'Figure4b_quantiles_Globe_1880_2023.eps'), ...
    'ContentType','vector');

exportgraphics(f4b, ...
    fullfile(FIG_DIR,'Figure4b_quantiles_Globe_1880_2023.png'), ...
    'Resolution',600,'BackgroundColor','white');

% -------------------------------------------------------------------------
% 8. Joint Figure 4: characteristics and quantiles
% -------------------------------------------------------------------------

figW = 7.0;
figH = 6.1;

fB = figure('Color','w','Units','inches','Position',[1 1 figW figH]);

tlB = tiledlayout(fB,5,2,'TileSpacing','compact','Padding','tight');

for i = 1:8

    ax = nexttile(tlB);

    plot(ax, time, Y(:,i), 'LineWidth',1.1);
    grid(ax,'on');
    box(ax,'on');

    title(ax, upper(names{i}), ...
        'FontWeight','normal', ...
        'FontSize',9);

    set(ax,'FontSize',8);

    xticks(ax, time(tickx));
    xtickangle(ax,0);

    if mod(i,2) == 1
        ylabel(ax,'Temp (°C)','FontSize',8);
    else
        ax.YTickLabel = [];
    end

end

axQ = nexttile(tlB,[1 2]);

q_idx = 9:19;

plot(axQ, time, Y(:,q_idx), 'LineWidth',1.0);

grid(axQ,'on');
box(axQ,'on');
set(axQ,'FontSize',8);

xticks(axQ, time(tickx));
xtickangle(axQ,0);

xlabel(axQ,'Year','FontSize',9);
ylabel(axQ,'Temp (°C)','FontSize',9);

q_labels = {'q05','q10','q20','q30','q40','q50','q60','q70','q80','q90','q95'};

lg = legend(axQ, q_labels, ...
    'Orientation','horizontal', ...
    'Box','off', ...
    'FontSize',8);

lg.NumColumns = 6;
lg.Layout.Tile = 'south';

title(tlB, ...
    'Distributional characteristics and quantiles of global temperature (CRU 1880–2023)', ...
    'FontSize',11);

set(fB,'ToolBar','none','MenuBar','none');

axall = findall(fB,'Type','axes');

for k = 1:numel(axall)
    try
        axtoolbar(axall(k),'Visible','off');
    catch
    end
    if exist('disableDefaultInteractivity','file')
        disableDefaultInteractivity(axall(k));
    end
end

exportgraphics(fB, ...
    fullfile(FIG_DIR,'Figure4_joint_Globe_1880_2023.pdf'), ...
    'ContentType','vector','BackgroundColor','white');

exportgraphics(fB, ...
    fullfile(FIG_DIR,'Figure4_joint_Globe_1880_2023.eps'), ...
    'ContentType','vector');

exportgraphics(fB, ...
    fullfile(FIG_DIR,'Figure4_joint_Globe_1880_2023.png'), ...
    'Resolution',600,'BackgroundColor','white');

fprintf('\n============================================================\n');
fprintf('Trend tests successfully computed.\n');
fprintf('Results folder: %s\n', RES_DIR);
fprintf('Table file    : %s\n', file);
fprintf('Figures folder: %s\n', FIG_DIR);
fprintf('============================================================\n');