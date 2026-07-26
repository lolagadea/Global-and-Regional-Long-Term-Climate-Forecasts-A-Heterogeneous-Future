% ========================================================================
% main_comb_all_models_q05.m
% ------------------------------------------------------------------------
% Long-term forecasts for the global q05 using all individual models and
% BIC-weighted forecast combinations.
% ========================================================================

clear; clc;
warning('off');

% ------------------------------------------------------------------------
% 1. Portable paths
% ------------------------------------------------------------------------

THIS_DIR = fileparts(mfilename('fullpath'));
FORECAST_DIR = fileparts(THIS_DIR);
GLOBE_DIR = fileparts(FORECAST_DIR);
ROOT_DIR = fileparts(GLOBE_DIR);

FUN_GENERAL_DIR  = fullfile(ROOT_DIR,'functions');
FUN_FORECAST_DIR = fullfile(ROOT_DIR,'functions_for');

COMP_RES_DIR = fullfile(FORECAST_DIR, ...
    'run_forecast_competition','Results_for');

FIG_DIR = fullfile(THIS_DIR,'Figures');
TAB_DIR = fullfile(THIS_DIR,'Tables');
RES_DIR = fullfile(THIS_DIR,'Results');

if ~isfolder(FIG_DIR), mkdir(FIG_DIR); end
if ~isfolder(TAB_DIR), mkdir(TAB_DIR); end
if ~isfolder(RES_DIR), mkdir(RES_DIR); end

addpath(genpath(FUN_GENERAL_DIR));
addpath(genpath(FUN_FORECAST_DIR));
addpath(genpath(COMP_RES_DIR));

% ------------------------------------------------------------------------
% 2. Load data and BIC weights
% ------------------------------------------------------------------------

DATA_FILE = fullfile(ROOT_DIR, ...
    '02_Globe','Introduction','Results', ...
    'QUANTILES_monthly_Globe_1880_2023.mat');

assert(isfile(DATA_FILE), ...
    'Input file QUANTILES_monthly_Globe_1880_2023.mat not found.');

load(DATA_FILE,'QUANTILES_monthly');

Weights_BIC = load(fullfile(COMP_RES_DIR,'Weights_BIC.dat'));
BIC         = load(fullfile(COMP_RES_DIR,'BIC.dat'));

q_index = 9;

Y = QUANTILES_monthly.Globe(:,q_index);

[t,n] = size(Y);

pmax   = 12;
maxpol = 12;
mod    = 14;

S = 1;
F = [1,25,50,100];

RDOS = NaN(mod,length(S),length(F));
CI   = NaN(mod,length(S),length(F),2);

comb0 = NaN(length(S),length(F));
comb1 = NaN(length(S),length(F));
comb2 = NaN(length(S),length(F));
comb3 = NaN(length(S),length(F));

RES_models = NaN(t,mod,length(S),length(F));

% ------------------------------------------------------------------------
% 3. Forecasts
% ------------------------------------------------------------------------

for u = 1:length(S)
    
    y = Y(:,S(u));
    
    for v = 1:length(F)
        
        f = F(v);
        
        % Mean model
        [b,t_nw,se_nw,res,r2,varres,varBhat,y_hat] = ...
            ols_hac_forecast(y,ones(t,1));
        
        xf_mean = b;
        ci_mean = xf_mean + norminv([0.05 0.95])*sqrt(varres);
        
        RDOS(1,u,v) = xf_mean;
        CI(1,u,v,:) = ci_mean;
        RES_models(:,1,u,v) = res;
        
        % Linear trend model
        k = 1;
        [beta,t_ratio,y_hat,varres,varbeta] = estima_pol_trend_hac(y,k);
        res = y-y_hat;
        
        trend_f = [1,t+f];
        xf_linear_trend = beta'*trend_f';
        ci_linear_trend = xf_linear_trend + norminv([0.05 0.95])*sqrt(varres);
        
        RDOS(2,u,v) = xf_linear_trend;
        CI(2,u,v,:) = ci_linear_trend;
        RES_models(:,2,u,v) = res;
        
        % Polynomial trend model
        k = estima_select_pol_trend_hac_sbic(y,maxpol);
        [beta,t_ratio,y_hat,varres,varbeta] = estima_pol_trend_hac(y,k);
        res = y-y_hat;
        slope = compute_slope(beta,t);
        
        if k > 0
            trend_f = [];
            for ii = 0:k
                trend_f = [trend_f,(t+f)^ii];
            end
            
            xf_pol_trend = beta'*trend_f';
            ci_pol_trend = xf_pol_trend + norminv([0.05 0.95])*sqrt(varres);
            
            beta2 = [beta(1),slope];
            trend_f2 = [1,(t+f)];
            xf2_pol_trend = beta2*trend_f2';
            ci2_pol_trend = xf2_pol_trend + norminv([0.05 0.95])*sqrt(varres);
            
        elseif k == 0
            xf_pol_trend = beta;
            xf2_pol_trend = beta;
            ci_pol_trend = xf_pol_trend + norminv([0.05 0.95])*sqrt(varres);
            ci2_pol_trend = xf2_pol_trend + norminv([0.05 0.95])*sqrt(varres);
        end
        
        RDOS(3,u,v) = xf_pol_trend;
        RDOS(4,u,v) = xf2_pol_trend;
        CI(3,u,v,:) = ci_pol_trend;
        CI(4,u,v,:) = ci2_pol_trend;
        RES_models(:,3,u,v) = res;
        RES_models(:,4,u,v) = res;
        
        % Polynomial trend in logs
        k = estima_select_pol_trend_log_hac_sbic(y,maxpol);
        [beta,t_ratio,y_hat,varres,varbeta] = estima_pol_trend_log_hac(y,k);
        res = y-y_hat;
        
        if k > 0
            trend_f = [];
            for ii = 0:k
                trend_f = [trend_f,log(t+f)^ii];
            end
            xf_pol_trend_log = beta'*trend_f';
        elseif k == 0
            xf_pol_trend_log = beta;
        end
        
        ci_pol_trend_log = xf_pol_trend_log + norminv([0.05 0.95])*sqrt(varres);
        
        RDOS(5,u,v) = xf_pol_trend_log;
        CI(5,u,v,:) = ci_pol_trend_log;
        RES_models(:,5,u,v) = res;
        
        % Structural-break model
        model = 3;
        criteria = 2;
        eps = 0.15;
        kmax = fix(12*(t/100)^(1/4));
        
        [wald,cv,tb] = qfgls_simus(y,kmax,model,criteria,eps);
        
        if wald > cv(2)
            [beta_break,t_nw,se_nw,res,r2,varres,varBhat,y_hat] = ...
                estima_break_model(y,tb,model);
            
            xf_struct_breaks = beta_break(1)+beta_break(2)+beta_break(4)*(t+f);
            ci_struct_breaks = xf_struct_breaks + norminv([0.05 0.95])*sqrt(varres);
        else
            k = 1;
            [beta,t_ratio,y_hat,varres,varbeta] = estima_pol_trend_hac(y,k);
            res = y-y_hat;
            trend_f = [1,t+f];
            xf_struct_breaks = beta'*trend_f';
            ci_struct_breaks = xf_struct_breaks + norminv([0.05 0.95])*sqrt(varres);
        end
        
        RDOS(6,u,v) = xf_struct_breaks;
        CI(6,u,v,:) = ci_struct_breaks;
        RES_models(:,6,u,v) = res;
        
        % Polynomial trend + AR(p)
        [p,k] = estima_select_pol_trend_arp_hac_sbic(y,maxpol,pmax);
        [beta,t_ratio,y_hat,r2,varres,varBhat] = ...
            estima_pol_trend_arp_sbic_direct_h(y,k,p,f);
        
        res = y(p+f+1:end)-y_hat;
        
        if p > 0 && k > 0
            z = [];
            for ww = 0:p-1
                z = [z,y(t-ww)];
            end
            
            trend_f = [];
            for ii = 0:k
                trend_f = [trend_f,(t+f)^ii];
            end
            
            Z = [trend_f';z'];
            xf_pol_trend_arp = beta'*Z;
            
        elseif p == 0 && k == 0
            xf_pol_trend_arp = beta;
            
        elseif p == 0 && k > 0
            trend_f = [];
            for ii = 0:k
                trend_f = [trend_f,(t+f)^ii];
            end
            xf_pol_trend_arp = beta'*trend_f';
            
        elseif p > 0 && k == 0
            z = [];
            for ww = 0:p-1
                z = [z,y(t-ww)];
            end
            z = [1,z];
            xf_pol_trend_arp = z*beta;
        end
        
        ci_pol_trend_arp = xf_pol_trend_arp + norminv([0.05 0.95])*sqrt(varres);
        
        RDOS(7,u,v) = xf_pol_trend_arp;
        CI(7,u,v,:) = ci_pol_trend_arp;
        RES_models(p+f+1:t,7,u,v) = res;
        
        % Polynomial trend + AR(p), average-slope version
        [p,k] = estima_select_pol_trend_arp_hac_sbic(y,maxpol,pmax);
        [beta,t_ratio,y_hat,r2,varres,varBhat] = ...
            estima_pol_trend_arp_sbic_direct_h(y,k,p,f);
        
        res = y(p+f+1:end)-y_hat;
        
        if p > 0 && k > 0
            z = [];
            for ww = 0:p-1
                z = [z,y(t-ww)];
            end
            z = z';
            
            slope = compute_slope(beta(1:k+1),t);
            beta2 = [beta(1);slope;beta(k+2:end)];
            trend_f = [1,(t+f)];
            Z = [trend_f';z];
            xf_pol_trend_arp_av_sl = beta2'*Z;
            
        elseif p == 0 && k == 0
            xf_pol_trend_arp_av_sl = beta;
            
        elseif p == 0 && k > 0
            slope = compute_slope(beta(1:k+1),t);
            trend_f = [1,(t+f)];
            xf_pol_trend_arp_av_sl = [beta(1),slope]*trend_f';
            
        elseif p > 0 && k == 0
            z = 1;
            for vv = 0:p-1
                z = [z;y(t-vv)];
            end
            xf_pol_trend_arp_av_sl = beta'*z;
        end
        
        ci_pol_trend_arp_av_sl = xf_pol_trend_arp_av_sl + ...
            norminv([0.05 0.95])*sqrt(varres);
        
        RDOS(8,u,v) = xf_pol_trend_arp_av_sl;
        CI(8,u,v,:) = ci_pol_trend_arp_av_sl;
        RES_models(p+f+1:t,8,u,v) = res;
        
        % AR(p)
        [~,~,p] = estima_select_arp_sbic(y,pmax);
        [beta,y_hat,varres] = estima_arp_direct_h(y,p,f);
        res = y(p+f+1:end)-y_hat;
        
        if p > 0
            z = [];
            for ww = 0:p-1
                z = [z,y(t-ww)];
            end
            z = [1,z];
            xf_arp = z*beta;
        elseif p == 0
            xf_arp = beta;
        end
        
        ci_arp = xf_arp + norminv([0.05 0.95])*sqrt(varres);
        
        RDOS(9,u,v) = xf_arp;
        CI(9,u,v,:) = ci_arp;
        RES_models(p+f+1:t,9,u,v) = res;
        
        % Random walk
        xf_rw = y(end);
        res = my_diff(y,f);
        varres = res'*res/(length(y)-1);
        ci_rw = xf_rw + norminv([0.05 0.95])*sqrt(varres);
        
        RDOS(10,u,v) = y(t);
        CI(10,u,v,:) = ci_rw;
        RES_models(f+1:t,10,u,v) = res;
        
        % Random walk with drift
        dy = diff(y,1);
        reg = ones(t-1,1);
        alpha = ols_hac_forecast(dy,reg);
        
        res = my_diff(y,f)-alpha*f;
        varres = res'*res/(length(y)-2);
        xf_rwd = y(t)+alpha*f;
        ci_rwd = xf_rwd + norminv([0.05 0.95])*sqrt(varres);
        
        RDOS(11,u,v) = xf_rwd;
        CI(11,u,v,:) = ci_rwd;
        RES_models(f+1:t,11,u,v) = res;
        
        % IMA model
        dy = diff(y,1);
        mu = mean(dy);
        
        ym = y-mean(y);
        dym = diff(ym,1);
        
        ima = armax(dym,[0 1]);
        theta = -ima.c(2);
        
        eps_ima = y(t);
        for kk = 1:t-1
            eps_ima = eps_ima + theta^(kk-1)*(theta-1)*y(t-kk);
        end
        
        if theta < 0.97
            eps_ima = eps_ima - mu/(1-theta);
            xf_ima = mu*f+y(t)-theta*eps_ima;
            res = my_diff(y,f)-mu*f+theta*eps_ima;
            
        elseif theta > 0.99
            xf_ima = y(1)+(t+f)*mu;
            res = my_diff(y,f)-mu*(t+f);
            
        else
            eps_ima = eps_ima - mu/(1-theta);
            xf_ima = y(1)+mu*f+y(t)-theta*eps_ima;
            res = my_diff(y,f)-mu*f+theta*eps_ima;
        end
        
        varres = res'*res/(length(y)-2);
        ci_ima = xf_ima + norminv([0.05 0.95])*sqrt(varres);
        
        RDOS(12,u,v) = xf_ima;
        CI(12,u,v,:) = ci_ima;
        RES_models(f+1:t,12,u,v) = res;
        
        % ARFIMA model
        [d,ar,sigma2,res,y_hat] = estima_arfima(y,0);
        xf = arfima_forecast(y,f,d,[],[],mean(y),sigma2);
        xf_arfima = xf(end);
        varres = res'*res;
        ci_arfima = xf_arfima + norminv([0.05 0.95])*sqrt(varres);
        
        RDOS(13,u,v) = xf_arfima;
        CI(13,u,v,:) = ci_arfima;
        RES_models(:,13,u,v) = res;
        
        % AR(20)
        p = 20;
        [beta,y_hat,varres] = estima_arp_direct_h(y,p,f);
        res = y(p+f+1:end)-y_hat;
        
        z = [];
        for ww = 0:p-1
            z = [z,y(t-ww)];
        end
        z = [1,z];
        xf_arp20 = z*beta;
        
        ci_arp20 = xf_arp20 + norminv([0.05 0.95])*sqrt(varres);
        
        RDOS(14,u,v) = xf_arp20;
        CI(14,u,v,:) = ci_arp20;
        RES_models(p+f+1:t,14,u,v) = res;
        
        % Forecast combinations
        
        comb0(u,v) = Weights_BIC(q_index,:)*squeeze(RDOS(:,u,v));
        Weig0 = Weights_BIC(q_index,:);
        
        % Combination 1: remove two extreme forecasts
        BIC0 = BIC(q_index,:);
        A = squeeze(RDOS(:,u,v));
        [~,maxpos] = max(A);
        [~,minpos] = min(A);
        A([minpos,maxpos]) = [];
        pos_comb1 = [minpos,maxpos];
        
        BIC1 = BIC0;
        BIC1([minpos,maxpos]) = [];
        
        Weig1 = NaN(1,length(BIC1));
        for ll = 1:length(BIC1)
            Weig1(ll) = exp(-0.5*BIC1(ll))/sum(exp(-0.5*BIC1));
        end
        
        comb1(u,v) = Weig1*A;
        
        % Combination 2: remove four extreme forecasts
        [~,maxpos] = max(A);
        [~,minpos] = min(A);
        A([minpos,maxpos]) = [];
        pos_comb2 = [minpos,maxpos];
        
        BIC2 = BIC1;
        BIC2([minpos,maxpos]) = [];
        
        Weig2 = NaN(1,length(BIC2));
        for ll = 1:length(BIC2)
            Weig2(ll) = exp(-0.5*BIC2(ll))/sum(exp(-0.5*BIC2));
        end
        
        comb2(u,v) = Weig2*A;
        
        % Combination 3: remove six extreme forecasts
        [~,maxpos] = max(A);
        [~,minpos] = min(A);
        A([minpos,maxpos]) = [];
        pos_comb3 = [minpos,maxpos];
        
        BIC3 = BIC2;
        BIC3([minpos,maxpos]) = [];
        
        Weig3 = NaN(1,length(BIC3));
        for ll = 1:length(BIC3)
            Weig3(ll) = exp(-0.5*BIC3(ll))/sum(exp(-0.5*BIC3));
        end
        
        comb3(u,v) = Weig3*A;
        
        % Confidence intervals for combinations
        CI_models(:,u,v,:) = squeeze(CI(:,u,v,:));
        
        RES = squeeze(RES_models(:,:,u,v));
        
        var_comb0 = compute_cov_comb(RES,Weig0);
        CI_comb0(u,v,:) = comb0(u,v) + ...
            [norminv(0.05)*sqrt(var_comb0),norminv(0.95)*sqrt(var_comb0)];
        
        RES1 = RES;
        RES1(:,pos_comb1) = [];
        var_comb1 = compute_cov_comb(RES1,Weig1);
        CI_comb1(u,v,:) = comb1(u,v) + ...
            [norminv(0.05)*sqrt(var_comb1),norminv(0.95)*sqrt(var_comb1)];
        
        RES2 = RES1;
        RES2(:,pos_comb2) = [];
        var_comb2 = compute_cov_comb(RES2,Weig2);
        CI_comb2(u,v,:) = comb2(u,v) + ...
            [norminv(0.05)*sqrt(var_comb2),norminv(0.95)*sqrt(var_comb2)];
        
        RES3 = RES2;
        RES3(:,pos_comb3) = [];
        var_comb3 = compute_cov_comb(RES3,Weig3);
        CI_comb3(u,v,:) = comb3(u,v) + ...
            [norminv(0.05)*sqrt(var_comb3),norminv(0.95)*sqrt(var_comb3)];
        
    end
    
    COMB0(u,:) = comb0;
    COMB1(u,:) = comb1;
    COMB2(u,:) = comb2;
    COMB3(u,:) = comb3;
    
    FORECAST_table = squeeze(RDOS);
    FORECAST_table = [FORECAST_table;COMB0;COMB1;COMB2;COMB3];
    
    CI_models = squeeze(CI_models);
    CI_fit_table = [CI_models;CI_comb0;CI_comb1;CI_comb2;CI_comb3];
    
    % --------------------------------------------------------------------
    % 4. LaTeX table
    % --------------------------------------------------------------------
    
    models = {'mean','linear-trend','pol-trend','pol-trend-av-sl', ...
        'pol-trend-log','struct-breaks','pol-trend-arp', ...
        'pol-trend-arp-av-sl','arp','rw','rwd','ima','arfima', ...
        'arp20','combined0','combined1','combined2','combined3'};
    
    TAB_FILE = fullfile(TAB_DIR, ...
        'table_forecast_comb_all_q05_Globe_1880_2023.tex');
    
    fid = fopen(TAB_FILE,'w');
    
    fprintf(fid,'\\begin{tabular}{lrrrr}\n');
    fprintf(fid,'\\toprule\n');
    fprintf(fid,'Models/horizon & h=1 & h=25 & h=50 & h=100 \\\\\n');
    fprintf(fid,'\\midrule\n');
    
    for ii = 1:length(models)
        fprintf(fid,'%s',models{ii});
        for jj = 1:length(F)
            fprintf(fid,' & %5.2f',FORECAST_table(ii,jj));
        end
        fprintf(fid,' \\\\\n');
        
        fprintf(fid,'');
        for jj = 1:length(F)
            fprintf(fid,' & (%5.2f,%5.2f)', ...
                CI_fit_table(ii,jj,1), CI_fit_table(ii,jj,2));
        end
        fprintf(fid,' \\\\\n');
    end
    
    fprintf(fid,'\\bottomrule\n');
    fprintf(fid,'\\end{tabular}\n');
    
    fclose(fid);
    
   
    
    % --------------------------------------------------------------------
    % 6. Save RDOS
    % --------------------------------------------------------------------
    
    RDOS_FILE = fullfile(RES_DIR,'RDOS_for_q05_Globe_1880_2023.mat');

if isfile(RDOS_FILE)
    Sload = load(RDOS_FILE);
    
    if isfield(Sload,'RDOS')
        RDOS_out = Sload.RDOS;
    else
        RDOS_out = struct();
    end
else
    RDOS_out = struct();
end
    
    RDOS_out.comb_all.H = F;
    RDOS_out.comb_all.models = models;
    RDOS_out.comb_all.forecast = FORECAST_table;
    RDOS_out.comb_all.ci = CI_fit_table;
    
    RDOS = RDOS_out;

    save(RDOS_FILE,'RDOS');
    
    % --------------------------------------------------------------------
    % 7. Figure
    % --------------------------------------------------------------------
    
    FORECAST = [COMB0;COMB1;COMB2;COMB3];
    FORECAST = FORECAST';
    
    figure(1)
    
    C = {'r','g','y',[0.5 0 0.5],'b'};
    
    for pp = 1:4
        
        subplot(2,2,pp)
        
        f = F(pp);
        years = 1880:2023+f;
        mcomb = 4;
        Z = NaN(t+F(end),mcomb);
        
        for ii = 1:mcomb
            Z(t+f,ii) = FORECAST(pp,ii);
        end
        
        plot([y;NaN(f,1)],'Color','b')
        hold on
        
        for ii = 1:mcomb
            plot(Z(:,ii), ...
                'Marker','o', ...
                'MarkerFaceColor',C{ii}, ...
                'Color',C{ii})
            hold on
        end
        
        hold off
        
        ylim([min(min(FORECAST(pp,:)),min(y))-1, ...
              max(max(FORECAST(pp,:)),max(y))+1])
        
        legend('original','Combined','Combined1','Combined2','Combined3', ...
            'Location','best');
        
        xticks = [1 20 40 60 80 100 120 t t+f];
        xticks = xticks(xticks <= length(years));
        
        set(gca,'XTick',xticks,'FontSize',12);
        set(gca,'XTickLabel',years(xticks),'FontSize',12);
        
        title(strcat('Forecast horizon',{' '},'h=',num2str(f)))
        
    end
    
    set(gcf,'PaperSize',[29.7 21.0], ...
        'PaperPosition',[0 0 29.7 21.0])
    
    figfile = fullfile(FIG_DIR,'Figure_for_q05_comb_all_Globe_1880_2023');
    
    print(gcf,figfile,'-dpdf');
    print(gcf,figfile,'-dpng');
    print(gcf,figfile,'-deps');
    
end

fprintf('\nCombined forecasts for the global q05 completed.\n');
fprintf('Tables saved in:\n%s\n', TAB_DIR);
fprintf('Figures saved in:\n%s\n', FIG_DIR);
fprintf('Results saved in:\n%s\n', RES_DIR);