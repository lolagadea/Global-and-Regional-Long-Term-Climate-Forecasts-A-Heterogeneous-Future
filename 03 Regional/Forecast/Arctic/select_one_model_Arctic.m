% ========================================================================
% select_one_model_Arctic.m
% ------------------------------------------------------------------------
% Select first-best / second-best Pareto-superior models for Arctic.
% ========================================================================

clear; clc;
clear global;
warning('off');

restoredefaultpath;
rehash toolboxcache;

% ------------------------------------------------------------------------
% 1. Portable paths
% ------------------------------------------------------------------------

THIS_DIR     = fileparts(mfilename('fullpath'));
FORECAST_DIR = fileparts(THIS_DIR);
REGIONAL_DIR = fileparts(FORECAST_DIR);
ROOT_DIR     = fileparts(REGIONAL_DIR);

FUN_DIR     = fullfile(ROOT_DIR,'functions');
FUN_FOR_DIR = fullfile(ROOT_DIR,'functions_for');

addpath(genpath(FUN_DIR));
addpath(genpath(FUN_FOR_DIR));

FIG_DIR = fullfile(THIS_DIR,'Figures');
RES_DIR = fullfile(THIS_DIR,'Results_for');

if ~isfolder(FIG_DIR), mkdir(FIG_DIR); end
if ~isfolder(RES_DIR), mkdir(RES_DIR); end

% ------------------------------------------------------------------------
% 2. Settings
% ------------------------------------------------------------------------

mod = 14;

models = { ...
    'mean', ...
    'linear-trend', ...
    'pol-trend', ...
    'pol-trend-av-sl', ...
    'pol-trend-log', ...
    'struct-break', ...
    'pol-trend-arp', ...
    'pol-trend-arp-av-sl', ...
    'arp', ...
    'rw', ...
    'rwd', ...
    'ima', ...
    'arfima', ...
    'arp20'};

Q = 11;

quantiles = { ...
    'q05','q10','q20','q30','q40','q50', ...
    'q60','q70','q80','q90','q95'};

F = [1,10,25];
W = [25];

M  = NaN(mod,length(W),length(F));
M0 = NaN(mod,length(W),length(F));
M1 = NaN(mod,length(W),length(F));
M2 = NaN(mod,length(W),length(F));

% ------------------------------------------------------------------------
% 3. Select one model
% ------------------------------------------------------------------------

ifig = 1;

for i = 1:length(F)

    f = F(i);

    for j = 1:length(W)

        w = W(j);

        f2 = strcat('for',num2str(f));
        w2 = strcat('w',num2str(w));

        GW_FILE = fullfile(RES_DIR, ...
            strcat('GW_rdos_',f2,'_',w2,'_all_quantiles.mat'));

        assert(isfile(GW_FILE), ...
            'GW file not found: %s', GW_FILE);

        load(GW_FILE,'RDOS');

        m = nansum(RDOS.GW.models(9:19,:)); % select only quantiles

        % ------------------------------------------------------------
        % Heatmap
        % ------------------------------------------------------------

        figure(ifig)

        MAT = RDOS.GW.models(9:19,:);

        my_heatmap(MAT,[],[],[], ...
            'TickAngle',45, ...
            'Colormap','red', ...
            'MinColorValue',0, ...
            'MaxColorValue',1, ...
            'ColorLevels',64);

        set(gca,'XTick',1:mod);
        set(gca,'XTickLabel',models,'Fontsize',12);

        set(gca,'YTick',1:Q);
        set(gca,'YTickLabel',quantiles, ...
            'Fontsize',12, ...
            'XTickLabelRotation',90);

        title(strcat( ...
            'Decision matrix for Pareto-superior models', ...
            '{ }(',w2,'{ }',f2,')'), ...
            'Fontsize',12);

        set(gcf, ...
            'Papersize',[29.7 21.0], ...
            'PaperPosition',[0 0 29.7 21.0]);

        print(figure(ifig), ...
            fullfile(FIG_DIR, ...
            strcat('Figure_heatmap_one_model_',w2,'_',f2)), ...
            '-dpdf');

        print(figure(ifig), ...
            fullfile(FIG_DIR, ...
            strcat('Figure_heatmap_one_model_',w2,'_',f2)), ...
            '-deps');

        print(figure(ifig), ...
            fullfile(FIG_DIR, ...
            strcat('Figure_heatmap_one_model_',w2,'_',f2)), ...
            '-dpng');

        ifig = ifig + 1;

        % ------------------------------------------------------------
        % First-best / second-best models
        % ------------------------------------------------------------

        M1(:,j,i) = m;

        pos = find(m == max(m));

        M2(pos,j,i) = 1;

        for mm = 1:mod

            if M1(mm,j,i) == Q
                M0(mm,j,i) = 1;
            end

        end

        if sum(isnan(M0(:,j,i))) < mod
            M(:,j,i) = M0(:,j,i);   % first best
        else
            M(:,j,i) = M2(:,j,i);   % second best
        end

    end

end

% ------------------------------------------------------------------------
% 4. Save output
% ------------------------------------------------------------------------

RDOS_one_model.select     = M;
RDOS_one_model.firstbest  = M0;
RDOS_one_model.count      = M1;
RDOS_one_model.secondbest = M2;
RDOS_one_model.F          = F;
RDOS_one_model.W          = W;

RDOS_one_model.structure = { ...
    'models', ...
    'window', ...
    'horizon'};

save(fullfile(RES_DIR,'RDOS_one_model.mat'),'RDOS_one_model');

fprintf('\nRDOS_one_model saved successfully.\n');

clear global;