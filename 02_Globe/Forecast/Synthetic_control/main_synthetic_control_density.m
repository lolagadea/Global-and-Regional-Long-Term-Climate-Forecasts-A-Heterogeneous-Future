% ========================================================================
% main_synthetic_control_density.m
% ------------------------------------------------------------------------
% Quasi-synthetic control exercise for global temperature quantiles.
%
% The script uses the pre-1961 sample (1880--1960) as the training period
% and produces counterfactual density forecasts for 1961--2023. These
% counterfactual quantiles are compared with the observed temperature
% quantiles in the full sample.
%
% Main output:
%   Figures/Figure_SC_density_Globe_1880_1960.*
%   Figures/Figure_SC_density_Globe_1880_1960_diff.*
%   Results/SC_density_Globe_1880_1960.mat
%
% Notes:
%   This script reproduces the density-based quasi-synthetic control figure
%   reported in the appendix of the paper.
% ========================================================================

clear; clc;
warning('off');

% ------------------------------------------------------------------------
% 1. Portable paths
% ------------------------------------------------------------------------

THIS_DIR = fileparts(mfilename('fullpath'));

% Synthetic_control is inside Forecast
FORECAST_DIR = fileparts(THIS_DIR);

GLOBE_DIR = fileparts(FORECAST_DIR);

ROOT_DIR = fileparts(GLOBE_DIR);

FUN_GENERAL_DIR  = fullfile(ROOT_DIR, 'functions');
FUN_FORECAST_DIR = fullfile(ROOT_DIR, 'functions_for');

INTRO_RES_DIR = fullfile(GLOBE_DIR, ...
    'Introduction', 'Results');

FIG_DIR = fullfile(THIS_DIR, 'Figures');
RES_DIR = fullfile(THIS_DIR, 'Results');
TAB_DIR = fullfile(THIS_DIR, 'Tables');

if ~isfolder(FIG_DIR), mkdir(FIG_DIR); end
if ~isfolder(RES_DIR), mkdir(RES_DIR); end
if ~isfolder(TAB_DIR), mkdir(TAB_DIR); end

addpath(genpath(FUN_GENERAL_DIR));
addpath(genpath(FUN_FORECAST_DIR));

% ------------------------------------------------------------------------
% 2. Load data
% ------------------------------------------------------------------------

DATA_FILE = fullfile(INTRO_RES_DIR, ...
    'QUANTILES_monthly_Globe_1880_2023.mat');

assert(isfile(DATA_FILE), ...
    'Input file not found: %s', DATA_FILE);

load(DATA_FILE, 'QUANTILES_monthly', 'name_labels');

% Full sample selected quantiles q05--q95
Yr = QUANTILES_monthly.Globe(:,9:19);

% Training sample: 1880--1960
Y = QUANTILES_monthly.Globe(1:81,9:19);

[t,n] = size(Y);

assert(t == 81, 'Training sample should contain 81 observations: 1880--1960.');
assert(n == 11, 'Expected 11 selected quantiles: q05--q95.');

% Quantile labels
if exist('name_labels','var') && numel(name_labels) >= 19
    Q = name_labels(9:19);
else
    Q = {'q05','q10','q20','q30','q40','q50','q60','q70','q80','q90','q95'};
end

% ------------------------------------------------------------------------
% 3. Forecast settings
% ------------------------------------------------------------------------

pmax   = 12;
maxpol = 12;
mod    = 14;

S = 1:n;       % selected quantiles q05--q95
F = 1:63;      % horizons 1961--2023

RDOS       = NaN(mod,length(S),length(F));
CI         = NaN(mod,length(S),length(F),2);
RES_models = NaN(t,mod,length(S),length(F));

% In the original code, model 11 is the selected Pareto model used for the
% density-based synthetic control figure.
selected_model_Pareto = 11;

% ------------------------------------------------------------------------
% 4. Forecast models estimated on the training sample
% ------------------------------------------------------------------------

for u = 1:length(S)

    fprintf('Quantile %d of %d\n', u, length(S));
    y = Y(:,S(u));

    for v = 1:length(F)

        f = F(v);

        % ----------------------------------------------------------------
        % 1. Mean model
        % ----------------------------------------------------------------
        [b,~,~,res,~,varres,~,~] = ols_hac_forecast(y,ones(t,1));
        xf_mean = b;
        ci_mean = xf_mean + norminv([0.05 0.95])*sqrt(varres);

        RDOS(1,u,v) = xf_mean;
        CI(1,u,v,:) = ci_mean;
        RES_models(:,1,u,v) = res;

        % ----------------------------------------------------------------
        % 2. Linear trend model
        % ----------------------------------------------------------------
        k = 1;
        [beta,~,y_hat,varres,~] = estima_pol_trend_hac(y,k);
        res = y - y_hat;
        trend_f = [1,t+f];
        xf_linear_trend = beta'*trend_f';
        ci_linear_trend = xf_linear_trend + norminv([0.05 0.95])*sqrt(varres);

        RDOS(2,u,v) = xf_linear_trend;
        CI(2,u,v,:) = ci_linear_trend;
        RES_models(:,2,u,v) = res;

        % ----------------------------------------------------------------
        % 3. Polynomial trend model
        % ----------------------------------------------------------------
        k = estima_select_pol_trend_hac_sbic(y,maxpol);
        [beta,~,y_hat,varres,~] = estima_pol_trend_hac(y,k);
        res = y - y_hat;
        slope = compute_slope(beta,t);

        if k > 0
            trend_f = [];
            for ii = 0:k
                trend_f = [trend_f,(t+f)^ii]; %#ok<AGROW>
            end

            xf_pol_trend = beta'*trend_f';
            ci_pol_trend = xf_pol_trend + norminv([0.05 0.95])*sqrt(varres);

            beta2 = [beta(1),slope];
            trend_f2 = [1,(t+f)];
            xf2_pol_trend = beta2*trend_f2';
            ci2_pol_trend = xf2_pol_trend + norminv([0.05 0.95])*sqrt(varres);

        elseif k == 0
            xf_pol_trend  = beta;
            xf2_pol_trend = beta;
            ci_pol_trend  = xf_pol_trend  + norminv([0.05 0.95])*sqrt(varres);
            ci2_pol_trend = xf2_pol_trend + norminv([0.05 0.95])*sqrt(varres);
        end

        RDOS(3,u,v) = xf_pol_trend;
        RDOS(4,u,v) = xf2_pol_trend;
        CI(3,u,v,:) = ci_pol_trend;
        CI(4,u,v,:) = ci2_pol_trend;
        RES_models(:,3,u,v) = res;
        RES_models(:,4,u,v) = res;

        % ----------------------------------------------------------------
        % 5. Polynomial trend in logs
        % ----------------------------------------------------------------
        k = estima_select_pol_trend_log_hac_sbic(y,maxpol);
        [beta,~,y_hat,varres,~] = estima_pol_trend_log_hac(y,k);
        res = y - y_hat;

        if k > 0
            trend_f = [];
            for ii = 0:k
                trend_f = [trend_f,log(t+f)^ii]; %#ok<AGROW>
            end
            xf_pol_trend_log = beta'*trend_f';
        elseif k == 0
            xf_pol_trend_log = beta;
        end

        ci_pol_trend_log = xf_pol_trend_log + norminv([0.05 0.95])*sqrt(varres);

        RDOS(5,u,v) = xf_pol_trend_log;
        CI(5,u,v,:) = ci_pol_trend_log;
        RES_models(:,5,u,v) = res;

        % ----------------------------------------------------------------
        % 6. Structural-break model
        % ----------------------------------------------------------------
        model = 3;
        criteria = 2;
        eps = 0.15;
        kmax = fix(12*(t/100)^(1/4));

        [wald,cv,tb] = qfgls_simus(y,kmax,model,criteria,eps);

        if wald > cv(2)
            [beta_break,~,~,res,~,varres,~,~] = estima_break_model(y,tb,model);
            xf_struct_breaks = beta_break(1) + beta_break(2) + beta_break(4)*(t+f);
            ci_struct_breaks = xf_struct_breaks + norminv([0.05 0.95])*sqrt(varres);
        else
            k = 1;
            [beta,~,y_hat,varres,~] = estima_pol_trend_hac(y,k);
            res = y - y_hat;
            trend_f = [1,t+f];
            xf_struct_breaks = beta'*trend_f';
            ci_struct_breaks = xf_struct_breaks + norminv([0.05 0.95])*sqrt(varres);
        end

        RDOS(6,u,v) = xf_struct_breaks;
        CI(6,u,v,:) = ci_struct_breaks;
        RES_models(:,6,u,v) = res;

        % ----------------------------------------------------------------
        % 7. Polynomial trend + AR(p)
        % ----------------------------------------------------------------
        [p,k] = estima_select_pol_trend_arp_hac_sbic(y,maxpol,pmax);
        [beta,~,y_hat,~,varres,~] = estima_pol_trend_arp_sbic_direct_h(y,k,p,f);
        res = y(p+f+1:end) - y_hat;

        if p > 0 && k > 0
            z = [];
            for jj = 0:p-1
                z = [z,y(t-jj)]; %#ok<AGROW>
            end

            trend_f = [];
            for ii = 0:k
                trend_f = [trend_f,(t+f)^ii]; %#ok<AGROW>
            end

            Z = [trend_f';z'];
            xf_pol_trend_arp = beta'*Z;

        elseif p == 0 && k == 0
            xf_pol_trend_arp = beta;

        elseif p == 0 && k > 0
            trend_f = [];
            for ii = 0:k
                trend_f = [trend_f,(t+f)^ii]; %#ok<AGROW>
            end
            xf_pol_trend_arp = beta'*trend_f';

        elseif p > 0 && k == 0
            z = [];
            for jj = 0:p-1
                z = [z,y(t-jj)]; %#ok<AGROW>
            end
            z = [1,z];
            xf_pol_trend_arp = z*beta;
        end

        ci_pol_trend_arp = xf_pol_trend_arp + norminv([0.05 0.95])*sqrt(varres);
        RDOS(7,u,v) = xf_pol_trend_arp;
        CI(7,u,v,:) = ci_pol_trend_arp;
        RES_models(p+f+1:t,7,u,v) = res;

        % ----------------------------------------------------------------
        % 8. Polynomial trend + AR(p), average slope
        % ----------------------------------------------------------------
        [p,k] = estima_select_pol_trend_arp_hac_sbic(y,maxpol,pmax);
        [beta,~,y_hat,~,varres,~] = estima_pol_trend_arp_sbic_direct_h(y,k,p,f);
        res = y(p+f+1:end) - y_hat;

        if p > 0 && k > 0
            z = [];
            for jj = 0:p-1
                z = [z,y(t-jj)]; %#ok<AGROW>
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
            for jj = 0:p-1
                z = [z;y(t-jj)]; %#ok<AGROW>
            end
            xf_pol_trend_arp_av_sl = beta'*z;
        end

        ci_pol_trend_arp_av_sl = xf_pol_trend_arp_av_sl + norminv([0.05 0.95])*sqrt(varres);
        RDOS(8,u,v) = xf_pol_trend_arp_av_sl;
        CI(8,u,v,:) = ci_pol_trend_arp_av_sl;
        RES_models(p+f+1:t,8,u,v) = res;

        % ----------------------------------------------------------------
        % 9. AR(p)
        % ----------------------------------------------------------------
        [~,~,p] = estima_select_arp_sbic(y,pmax);
        [beta,y_hat,varres] = estima_arp_direct_h(y,p,f);
        res = y(p+f+1:end) - y_hat;

        if p > 0
            z = [];
            for jj = 0:p-1
                z = [z,y(t-jj)]; %#ok<AGROW>
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

        % ----------------------------------------------------------------
        % 10. Random walk
        % ----------------------------------------------------------------
        xf_rw = y(end);
        res = my_diff(y,f);
        varres = res'*res/(length(y)-1);
        ci_rw = xf_rw + norminv([0.05 0.95])*sqrt(varres);

        RDOS(10,u,v) = y(t);
        CI(10,u,v,:) = ci_rw;
        RES_models(f+1:t,10,u,v) = res;

        % ----------------------------------------------------------------
        % 11. Random walk with drift
        % ----------------------------------------------------------------
        dy = diff(y,1);
        reg = ones(t-1,1);
        alpha = ols_hac(dy,reg);
        res = my_diff(y,f) - alpha*f;
        varres = res'*res/(length(y)-2);
        xf_rwd = y(t) + alpha*f;
        ci_rwd = xf_rwd + norminv([0.05 0.95])*sqrt(varres);

        RDOS(11,u,v) = xf_rwd;
        CI(11,u,v,:) = ci_rwd;
        RES_models(f+1:t,11,u,v) = res;

        % ----------------------------------------------------------------
        % 12. IMA model
        % ----------------------------------------------------------------
        dy = diff(y,1);
        mu = mean(dy);
        ym = y - mean(y);
        dym = diff(ym,1);
        ima = armax(dym,[0 1]);
        theta = -ima.c(2);
        eps_ima = y(t);

        for kk = 1:t-1
            eps_ima = eps_ima + theta^(kk-1)*(theta-1)*y(t-kk);
        end

        if theta < 0.97
            eps_ima = eps_ima - mu/(1-theta);
            xf_ima = mu*f + y(t) - theta*eps_ima;
            res = my_diff(y,f) - mu*f + theta*eps_ima;
        elseif theta > 0.99
            xf_ima = y(1) + (t+f)*mu;
            res = my_diff(y,f) - mu*(t+f);
        else
            eps_ima = eps_ima - mu/(1-theta);
            xf_ima = y(1) + mu*f + y(t) - theta*eps_ima;
            res = my_diff(y,f) - mu*f + theta*eps_ima;
        end

        varres = res'*res/(length(y)-2);
        ci_ima = xf_ima + norminv([0.05 0.95])*sqrt(varres);

        RDOS(12,u,v) = xf_ima;
        CI(12,u,v,:) = ci_ima;
        RES_models(f+1:t,12,u,v) = res;

        % ----------------------------------------------------------------
        % 13. ARFIMA model
        % ----------------------------------------------------------------
        [d,~,sigma2,res,~] = estima_arfima(y,0);
        xf = arfima_forecast(y,f,d,[],[],mean(y),sigma2);
        xf_arfima = xf(end);
        varres = res'*res;
        ci_arfima = xf_arfima + norminv([0.05 0.95])*sqrt(varres);

        RDOS(13,u,v) = xf_arfima;
        CI(13,u,v,:) = ci_arfima;
        RES_models(:,13,u,v) = res;

        % ----------------------------------------------------------------
        % 14. AR(20)
        % ----------------------------------------------------------------
        p = 20;
        [beta,y_hat,varres] = estima_arp_direct_h(y,p,f);
        res = y(p+f+1:end) - y_hat;

        z = [];
        for jj = 0:p-1
            z = [z,y(t-jj)]; %#ok<AGROW>
        end
        z = [1,z];
        xf_arp20 = z*beta;

        ci_arp20 = xf_arp20 + norminv([0.05 0.95])*sqrt(varres);
        RDOS(14,u,v) = xf_arp20;
        CI(14,u,v,:) = ci_arp20;
        RES_models(p+f+1:t,14,u,v) = res;

        % ----------------------------------------------------------------
        % Discard zeros
        % ----------------------------------------------------------------
        for mm = 1:mod
            if RDOS(mm,u,v) == 0
                RDOS(mm,u,v) = NaN;
                CI(mm,u,v,:) = [NaN, NaN];
                RES_models(:,mm,u,v) = NaN(t,1);
            end
        end

    end
end

% ------------------------------------------------------------------------
% 5. Synthetic-control summaries
% ------------------------------------------------------------------------

true_2023 = Yr(end,:)';
pred_2023 = RDOS(selected_model_Pareto,:,end)';

true_value = Yr(82:end,:);
for_value = squeeze(RDOS(selected_model_Pareto,:,:))';

q = length(Q);
D = NaN(1,q);

for ii = 1:q
    D(ii) = mean(true_value(:,ii) - for_value(:,ii), 'omitnan');
end

% Save results
OUT_FILE = fullfile(RES_DIR, 'SC_density_Globe_1880_1960.mat');
save(OUT_FILE, 'RDOS', 'CI', 'RES_models', 'selected_model_Pareto', ...
    'true_2023', 'pred_2023', 'true_value', 'for_value', 'D', 'Q');

% ------------------------------------------------------------------------
% 6. Figures
% ------------------------------------------------------------------------

figure(1); clf;
A = [true_2023, pred_2023];
barh(A);
legend('true value','predicted value','Location','best');
set(gca,'YTick',1:length(Q),'FontSize',12);
set(gca,'YTickLabel',Q,'FontSize',12);
title('Synthetic control exercise (density forecast)');

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[20 16]);
set(gcf,'PaperPosition',[0 0 20 16]);
set(findall(gcf,'type','axes'),'FontSize',14);

print(gcf, fullfile(FIG_DIR,'Figure_SC_density_Globe_1880_1960'), '-dpdf');
print(gcf, fullfile(FIG_DIR,'Figure_SC_density_Globe_1880_1960'), '-dpng');
print(gcf, fullfile(FIG_DIR,'Figure_SC_density_Globe_1880_1960'), '-depsc');

figure(2); clf;
bar(D);
set(gca, 'XTickLabel',Q, 'XTick',1:q);
title('Observed minus counterfactual differences, 1961--2023');

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[20 16]);
set(gcf,'PaperPosition',[0 0 20 16]);
set(findall(gcf,'type','axes'),'FontSize',14);

print(gcf, fullfile(FIG_DIR,'Figure_SC_density_Globe_1880_1960_diff'), '-dpdf');
print(gcf, fullfile(FIG_DIR,'Figure_SC_density_Globe_1880_1960_diff'), '-dpng');
print(gcf, fullfile(FIG_DIR,'Figure_SC_density_Globe_1880_1960_diff'), '-depsc');

fprintf('\nSynthetic control density exercise completed.\n');
fprintf('Results saved in:\n%s\n', OUT_FILE);
fprintf('Figures saved in:\n%s\n', FIG_DIR);
