function [QUANTILES, num_s, A, INDEX] = function_build_quantiles_months_method1(S, years, y1, y2, index)
% =========================================================================
% function_build_quantiles_months_method1.m
% -------------------------------------------------------------------------
% Purpose:
%   Construct annual distributional characteristics from CRUTEM station-month
%   units using Method 1.
%
% Method 1:
%   For each station-month unit, the function checks whether observations are
%   available for all years in the selected sample. If the station-month unit
%   is complete, it is retained. Annual distributional characteristics are
%   then computed cross-sectionally over the retained station-month units.
%
% Inputs:
%   S      : 3D array with dimensions years x stations x months
%   years  : vector of years
%   y1     : index of first year in the selected sample
%   y2     : index of last year in the selected sample
%   index  : optional station index used to select a subset of stations;
%            use NaN when no subset is imposed
%
% Outputs:
%   QUANTILES : matrix with annual distributional characteristics:
%               mean, max, min, std, iqr, range, kurtosis, skewness,
%               q05, q10, q20, q30, q40, q50, q60, q70, q80, q90, q95
%   num_s     : number of retained station-month units
%   A         : annual data matrix for retained station-month units
%   INDEX     : indicator for retained station-month units
%
% Notes:
%   - The original sample-selection rule is preserved.
%   - A station-month unit is retained only if it is observed for all years
%     in the selected sample.
%   - The function uses array_NaN.m to check completeness.
%
% =========================================================================

% -------------------------------------------------------------------------
% 1. Initial dimensions and optional station subset
% -------------------------------------------------------------------------

[t,n,m] = size(S); %#ok<ASGLU>

if isnan(index) == 0
    S2 = S(:,index,:);      % Optional subset of stations, e.g. polar areas
    n2 = size(S2,2);
else
    S2 = S;
    n2 = n;
end

% -------------------------------------------------------------------------
% 2. Select station-month units observed throughout the sample
% -------------------------------------------------------------------------

num_years = y2 - y1 + 1;

A     = NaN(num_years, n2*m);
INDEX = NaN(n2*m,1);
num_s = 0;

% This loop runs over station-month units. A unit is retained if it has
% valid information for all years in the selected sample.

for i = 1:n2*m

    s = S2(y1:y2,i);

    % Original alternatives, not used in the paper:
    % s = fillmissing(s,'linear','MaxGap',1);
    % s = fillmissing(s,'movmean',5,'MaxGap',1);

    if array_NaN(s) == 1
        num_s = num_s + 1;
        INDEX(i,1) = 1;
        A(:,i) = S2(y1:y2,i);
    end

end

% -------------------------------------------------------------------------
% 3. Initialize annual distributional characteristics
% -------------------------------------------------------------------------

m1      = NaN(num_years,1);
maximun = NaN(num_years,1);
minimun = NaN(num_years,1);
vol     = NaN(num_years,1);
iqr     = NaN(num_years,1);
rank    = NaN(num_years,1);
kur     = NaN(num_years,1);
skw     = NaN(num_years,1);

q05 = NaN(num_years,1);
q10 = NaN(num_years,1);
q20 = NaN(num_years,1);
q30 = NaN(num_years,1);
q40 = NaN(num_years,1);
q50 = NaN(num_years,1);
q60 = NaN(num_years,1);
q70 = NaN(num_years,1);
q80 = NaN(num_years,1);
q90 = NaN(num_years,1);
q95 = NaN(num_years,1);

N = NaN(num_years,1); %#ok<NASGU>

% -------------------------------------------------------------------------
% 4. Compute annual cross-sectional characteristics
% -------------------------------------------------------------------------

for k = 1:num_years

    z = A(k,:);
    z = z(~isnan(z));

    N(k) = length(z);

    if isempty(z) == 0

        m1(k)      = mean(z);
        maximun(k) = max(z);
        minimun(k) = min(z);
        vol(k)     = std(z);
        iqr(k)     = prctile(z,75) - prctile(z,25);
        rank(k)    = maximun(k) - minimun(k);
        skw(k)     = skewness(z);
        kur(k)     = kurtosis(z);

        q05(k) = prctile(z,5);
        q10(k) = prctile(z,10);
        q20(k) = prctile(z,20);
        q30(k) = prctile(z,30);
        q40(k) = prctile(z,40);
        q50(k) = prctile(z,50);
        q60(k) = prctile(z,60);
        q70(k) = prctile(z,70);
        q80(k) = prctile(z,80);
        q90(k) = prctile(z,90);
        q95(k) = prctile(z,95);

    end

end

% -------------------------------------------------------------------------
% 5. Collect outputs
% -------------------------------------------------------------------------

QUANTILES = [ ...
    m1, maximun, minimun, vol, iqr, rank, kur, skw, ...
    q05, q10, q20, q30, q40, q50, q60, q70, q80, q90, q95];

end