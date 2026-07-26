% ========================================================================
% select_one_model.m
% ------------------------------------------------------------------------
% Select models that are Pareto-superior across quantiles.
%
% Outputs:
%   Results_for/RDOS_one_model.mat
%   Figures_model_selection/Figure_heatmap_one_model_*.pdf/.eps/.png
% ========================================================================

clear; clc;
warning('off');

THIS_DIR = fileparts(mfilename('fullpath'));
FORECAST_DIR = fileparts(THIS_DIR);
GLOBE_DIR = fileparts(FORECAST_DIR);
ROOT_DIR = fileparts(GLOBE_DIR);

FUN_GENERAL_DIR  = fullfile(ROOT_DIR,'functions');
FUN_FORECAST_DIR = fullfile(ROOT_DIR,'functions_for');

RES_DIR = fullfile(THIS_DIR,'Results_for');
FIG_MODEL_DIR = fullfile(THIS_DIR,'Figures_model_selection');

if ~isfolder(RES_DIR), mkdir(RES_DIR); end
if ~isfolder(FIG_MODEL_DIR), mkdir(FIG_MODEL_DIR); end

addpath(genpath(FUN_GENERAL_DIR));
addpath(genpath(FUN_FORECAST_DIR));
addpath(genpath(RES_DIR));

mod = 14;

models = {'mean','linear-trend','pol-trend','pol-trend-av-sl', ...
    'pol-trend-log','struct-breaks','pol-trend-arp', ...
    'pol-trend-arp-av-sl','arp','rw','rwd','ima','arfima','arp20'};

Q = 11;

quantiles = {'q05','q10','q20','q30','q40','q50', ...
    'q60','q70','q80','q90','q95'};

W = [25,50];
F = [1,10,25];

M  = NaN(mod,length(W),length(F));
M0 = NaN(mod,length(W),length(F));
M1 = NaN(mod,length(W),length(F));
M2 = NaN(mod,length(W),length(F));

ifig = 1;

for i = 1:length(F)
    
    f = F(i);
    f2 = strcat('for',num2str(f));
    
    for j = 1:length(W)
        
        w = W(j);
        w2 = strcat('w',num2str(w));
        
        GW_FILE = fullfile(RES_DIR, ...
            strcat('GW_rdos_',f2,'_',w2,'_all_quantiles.mat'));
        
        assert(isfile(GW_FILE), ...
            ['Missing GW results file: ', GW_FILE]);
        
        load(GW_FILE,'RDOS');
        
        % Select only quantiles q05--q95
        MAT = RDOS.GW.models(9:19,:);
        model_count = nansum(MAT);
        
        figure(ifig);
        
        my_heatmap(MAT,[],[],[], ...
            'TickAngle',45, ...
            'Colormap','red', ...
            'MinColorValue',0, ...
            'MaxColorValue',1, ...
            'ColorLevels',64);
        
        set(gca,'XTick',1:1:mod);
        set(gca,'XTickLabel',models,'Fontsize',12);
        set(gca,'YTick',1:1:Q);
        set(gca,'YTickLabel',quantiles, ...
            'Fontsize',12, ...
            'XTickLabelRotation',90);
        
        title(strcat('Decision matrix for Pareto-superior models', ...
            '{ }','(',w2,'{ }',f2,')'), ...
            'Fontsize',12);
        
        set(gcf,'PaperUnits','centimeters');
        set(gcf,'PaperSize',[18 14]);
        set(gcf,'PaperPosition',[0 0 18 14]);
        set(gca,'FontSize',14);
        
        fig_base = fullfile(FIG_MODEL_DIR, ...
            strcat('Figure_heatmap_one_model_',w2,'_',f2));
        
        print(gcf,fig_base,'-dpdf');
        print(gcf,fig_base,'-deps');
        print(gcf,fig_base,'-dpng');
        
        ifig = ifig + 1;
        
        M1(:,j,i) = model_count;
        
        pos = find(model_count == max(model_count));
        M2(pos,j,i) = 1;
        
        for m = 1:mod
            if M1(m,j,i) == Q
                M0(m,j,i) = 1;
            end
        end
        
        if sum(isnan(M0(:,j,i))) < mod
            M(:,j,i) = M0(:,j,i);   % first best
        else
            M(:,j,i) = M2(:,j,i);   % second best
        end
        
    end
end

RDOS_one_model.select     = M;
RDOS_one_model.firstbest  = M0;
RDOS_one_model.count      = M1;
RDOS_one_model.secondbest = M2;
RDOS_one_model.F          = F;
RDOS_one_model.W          = W;
RDOS_one_model.models     = models;
RDOS_one_model.quantiles  = quantiles;
RDOS_one_model.structure  = {'models','window','horizon'};

save(fullfile(RES_DIR,'RDOS_one_model.mat'),'RDOS_one_model');

fprintf('\nOne-model selection completed.\n');
fprintf('Results saved in:\n%s\n', RES_DIR);
fprintf('Figures saved in:\n%s\n', FIG_MODEL_DIR);