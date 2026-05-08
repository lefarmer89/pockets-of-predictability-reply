function [cwBestMat, alphaAnnualBest, tStatAlphMatBest, SRBestAnnual] = ...
    adaptivePanelsForSplit(cfg)
% ADAPTIVEPANELSFORSPLIT  Per-day adaptive-shrinkage Panel A (CW) and
% Panel B (economic significance) statistics for a date-restricted
% evaluation window. Feeds the Adaptive (min RMSE) / Adaptive (max
% alpha) panels of Table 1 with cfg.sampleSplit = 'full'.
%
%   [cwBestMat, alphaAnnualBest, tStatAlphMatBest, SRBestAnnual] = ...
%       adaptivePanelsForSplit(cfg)
%
% Inputs:
%   cfg  config struct. cfg.sampleSplit ∈ {'full','pre1989','post1989'}.
%        cfg.shrinkageTarget defaults to 'zero' (bench=0 filter, the
%        published methodology).
%
% Outputs (each indexed by spec ii ∈ 1..9 and axis a ∈ {1=RMSE, 2=Alpha}):
%   cwBestMat[ii, c, a]      -- 9 x 3 x 2: Panel A CW stats. c ∈ {1=full,
%                               2=in-pocket, 3=out-of-pocket}.
%   alphaAnnualBest[a, ii]   -- 2 x 9: annualized alpha in percent
%   tStatAlphMatBest[a, ii]  -- 2 x 9: HAC t-stat of alpha
%   SRBestAnnual[a, ii]      -- 2 x 9: annualized Sharpe ratio
%
% The per-day adaptive selection uses full-sample expanding-window
% rankings (matching the published methodology); only the regression /
% aggregation is restricted to the date window. This is the cheap path
% described in computeAdaptiveSampleSplit.m -- correct because at any t
% the rank-1 model uses only data 1:t-1, so the existing rankings are
% valid OOS rankings to apply on the date-restricted sub-sample.
%
% Reads:
%   results/hyperparameters/forecastResults_2_2.5yE_1yDM_1S_pm_HyperR1.mat
%   results/hyperparameters/OOSResults_2_HyperR1_ExpandingC.mat
%
% Memory peak ~25 GB (yF2All filtered + alphMatExpanding loads + transients).

if nargin < 1 || isempty(cfg); cfg = default_config(); end
if ~isfield(cfg, 'sampleSplit') || isempty(cfg.sampleSplit)
    cfg.sampleSplit = 'full';
end
if ~ismember(cfg.sampleSplit, {'full', 'pre1989', 'post1989'})
    error('adaptivePanelsForSplit:badSplit', ...
        'cfg.sampleSplit must be ''full'', ''pre1989'', or ''post1989''.');
end

paths = cfg.paths;

% cfg.shrinkageTarget controls which slice of the 9720-combo grid the
% Adaptive selectors search over:
%   'zero' (default)   — bench=0 (zero-shrinkage target, 4860 combos,
%                        the published Adaptive panels of Table 1)
%   'prevailingMean'   — bench=1 (PM-shrinkage target, 4860 combos)
%   'any'              — both (full 9720 combos, Adaptive selectors
%                        compete across both shrinkage targets)
[paramCombsFull, idx] = constructHyperparameterGrid();
if isfield(cfg, 'shrinkageTarget') && strcmp(cfg.shrinkageTarget, 'prevailingMean')
    keepInd = find(paramCombsFull(:, idx.benchmark) == 1);
elseif isfield(cfg, 'shrinkageTarget') && strcmp(cfg.shrinkageTarget, 'any')
    keepInd = (1:size(paramCombsFull, 1))';
else
    keepInd = find(paramCombsFull(:, idx.benchmark) == 0);
end

%% Load H1: forecasts, pocket indicators, dates
robFolder = fullfile(paths.results, 'hyperparameters');
fileName1 = fullfile(robFolder, 'forecastResults_2_2.5yE_1yDM_1S_pm_HyperR1.mat');

S1 = load(fileName1, 'yActual', 'yF1PMMat', 'yF2All', 'pocketIndAll', 'dateVec', 'riskFree');
yActual      = S1.yActual ./ 100;
yF1PMMat     = S1.yF1PMMat ./ 100;
yF2All       = S1.yF2All(:, :, keepInd) ./ 100;
pocketIndAll = S1.pocketIndAll(:, :, keepInd);
dateVec_full = S1.dateVec;
riskFree     = S1.riskFree ./ 100;     % decimal form, matches yActual
clear S1;

trim = max(sum(isnan(yF2All)), [], 'all');
yActual      = yActual(trim+1:end);
yF1PMMat     = yF1PMMat(trim+1:end, :);
yF2All       = yF2All(trim+1:end, :, :);
pocketIndAll = pocketIndAll(trim+1:end, :, :);
pocketIndAll(isnan(pocketIndAll)) = 0;
pocketIndAll = logical(pocketIndAll);
dateVec      = dateVec_full(trim+1:end);
riskFree     = riskFree(trim+1:end);
T = numel(yActual);

% Transaction-cost adjustment in basis points (0 = published baseline,
% 5 / 10 = Table A.3 panels B / C). Threaded through panelBEcon below.
if isfield(cfg, 'adjCostBps') && ~isempty(cfg.adjCostBps)
    adjCostBps = cfg.adjCostBps;
else
    adjCostBps = 0;
end

%% Date mask
cutoff = datetime(1989, 1, 1);
switch cfg.sampleSplit
    case 'pre1989';  dateMask = dateVec <  cutoff;
    case 'post1989'; dateMask = dateVec >= cutoff;
    case 'full';     dateMask = true(size(dateVec));
    otherwise; error('unreachable');
end

%% Combinations
combPocket = squeeze(any(pocketIndAll(:, 1:4, :), 2));
yComb2 = combineInPocketAverageLocal(yF2All, pocketIndAll, yF1PMMat);
yComb3 = squeeze(mean(yF2All(:, 1:4, :), 2, 'omitmissing'));
yF2InPocket = yF2All .* pocketIndAll + yF1PMMat .* (~pocketIndAll);
yComb1 = squeeze(mean(yF2InPocket(:, 1:4, :), 2, 'omitmissing'));

%% Build yF2AllS9 (in-pocket-masked + 3 combos) for rmseAll
yF2AllS9 = NaN(T, 9, size(yF2All, 3));
yF2AllS9(:, 1:6, :) = yF2InPocket;
yF2AllS9(:, 7, :)   = yComb1;
yF2AllS9(:, 8, :)   = yComb2;
yF2AllS9(:, 9, :)   = yComb3;
clear yF2InPocket;

rmseAll = sqrt(cumsum(((yActual - yF2AllS9).^2) ./ ((1:T)')));

%% Load H2 alphMatExpanding / tStatAlphMatExpanding (filtered)
fileName2 = fullfile(robFolder, 'OOSResults_2_HyperR1_ExpandingC.mat');

S2 = load(fileName2, 'alphMatExpanding');
alphMatExpanding = S2.alphMatExpanding(:, :, keepInd);
clear S2;
S2 = load(fileName2, 'tStatAlphMatExpanding');
tStatAlphMatExpanding = S2.tStatAlphMatExpanding(:, :, keepInd);
clear S2;
T_h2 = size(alphMatExpanding, 1);
if T_h2 ~= T
    alphMatExpanding      = alphMatExpanding(end-T+1:end, :, :);
    tStatAlphMatExpanding = tStatAlphMatExpanding(end-T+1:end, :, :);
end

numSpecifications = size(alphMatExpanding, 2);

%% Per-day best indices (full-sample expanding-window rankings)
% Deterministic tie-break: when two combos tie at FP-precision (~1e-15),
% min/max returns the first index, which depends on the dim-3 ordering
% set by parfor scheduling in dailyEmpiricsHyperparameters{1,2}. We
% break ties by combo-index (lowest combo wins) so the choice is
% reproducible across runs regardless of worker count. The 1e-13
% perturbation is well above FP-noise but well below any real
% alpha / RMSE / t-stat difference (alpha O(0.01-0.1), RMSE O(0.1),
% t-stats O(1)).
nCombos = size(rmseAll, 3);
tieBreakAsc  =  1e-13 * reshape(1:nCombos, 1, 1, []);
tieBreakDesc = -1e-13 * reshape(1:nCombos, 1, 1, []);
[~, bestRMSEInd]      = min(rmseAll              + tieBreakAsc,  [], 3, 'omitmissing');
[~, bestAlphInd]      = max(alphMatExpanding      + tieBreakDesc, [], 3, 'omitmissing');
[~, bestTStatAlphInd] = max(tStatAlphMatExpanding + tieBreakDesc, [], 3, 'omitmissing');
bestRMSEInd      = squeeze(bestRMSEInd);
bestAlphInd      = squeeze(bestAlphInd);
bestTStatAlphInd = squeeze(bestTStatAlphInd);
clear rmseAll alphMatExpanding tStatAlphMatExpanding;

%% Panel A: CW stats, date-restricted aggregation
cwBestMat = panelACW( ...
    yActual, yF1PMMat, yF2All, yComb1, yComb2, yComb3, ...
    pocketIndAll, combPocket, ...
    bestRMSEInd, bestAlphInd, T, numSpecifications, dateMask);

%% Panel B: per-day-adaptive portfolio time series, date-restricted regression
[alphaAnnualBest, tStatAlphMatBest, SRBestAnnual] = panelBEcon( ...
    yF2AllS9, yF1PMMat(:, 1), yActual, riskFree, ...
    bestRMSEInd, bestAlphInd, ...
    T, numSpecifications, dateMask, adjCostBps);

end


% ====================================================================
%  Local helpers
% ====================================================================

function yComb2 = combineInPocketAverageLocal(yF2, pocketIndAll, yF1PMMat)
T = size(yF2, 1);
N = size(yF2, 3);
yComb2 = NaN(T, N);
for ii = 1:T
    for jj = 1:N
        if any(pocketIndAll(ii, 1:4, jj))
            yComb2(ii, jj) = mean(yF2(ii, logical(pocketIndAll(ii, 1:4, jj)), jj));
        else
            yComb2(ii, jj) = yF1PMMat(ii, 1);
        end
    end
end
end


function cwBestMat = panelACW( ...
    yActual, yF1PMMat, yF2All, yComb1, yComb2, yComb3, ...
    pocketIndAll, combPocket, ...
    bestRMSEInd, bestAlphInd, T, numSpecifications, dateMask)
% Per-spec, per-axis HAC t-stats for the per-day best combo's CW
% differential, restricted to dateMask. axes: 1=RMSE, 2=Alpha.
cwBestMat = NaN(numSpecifications - 1, 3, 2);
yF1pm1 = yF1PMMat(:, 1);
bestInds = {bestRMSEInd, bestAlphInd};

for ii = 1:(numSpecifications - 1)
    if ii <= 6
        v = ii;
        yv = squeeze(yF2All(:, v, :));
        srcCW = (yActual - yF1PMMat(:, v)).^2 - ...
                ((yActual - yv).^2 - (yF1PMMat(:, v) - yv).^2);
        srcPocket = squeeze(pocketIndAll(:, v, :));
    else
        switch ii
            case 7; yC = yComb1;
            case 8; yC = yComb2;
            case 9; yC = yComb3;
        end
        srcCW = (yActual - yF1pm1).^2 - ...
                ((yActual - yC).^2 - (yF1pm1 - yC).^2);
        srcPocket = combPocket;
    end

    for ax = 1:2
        bInd = bestInds{ax};
        fCWBest       = NaN(T, 1);
        pocketIndBest = false(T, 1);
        for t = 2:T
            fCWBest(t)       = srcCW(t, bInd(t-1, ii));
            pocketIndBest(t) = srcPocket(t, bInd(t-1, ii));
        end
        % Restrict to dateMask before stats
        fCW_d   = fCWBest(dateMask);
        pInd_d  = pocketIndBest(dateMask);
        rowCW = computeDMCWStats(fCW_d, pInd_d, 'cw');
        cwBestMat(ii, :, ax) = rowCW(1:3);
    end
end
end


function [alphaAnnualBest, tStatAlphMatBest, SRBestAnnual] = panelBEcon( ...
    yF2AllS9, yF1pm1, yActual, riskFree, ...
    bestRMSEInd, bestAlphInd, ...
    T, numSpecifications, dateMask, adjCostBps)
% Per-spec, per-axis: build the per-day-adaptive portfolio time series
% (using OOS c-scaling), restrict to dateMask, regress on yActual, then
% annualize alpha and Sharpe ratio. Mirrors H2's bestComboRegPanel
% pipeline but with the date-restriction applied at the regression step.
%
% Transaction-cost adjustment matches constructPortfolioFactors.m's
% daily (geometric self-financing) form so Table A.3's Adaptive columns
% scale with adjCostBps the same way the Fixed column does. When
% adjCostBps == 0, the subtraction is a no-op, recovering the published
% baseline.

if nargin < 10 || isempty(adjCostBps); adjCostBps = 0; end
adjCost = adjCostBps / 100 / 100;

weightRestriction = 1;
evalStartInd = 252;

% Warmup: benchmark-spec portfolio (constant across combos in spec=end).
cBench = sqrt(var(yActual, 'omitmissing') / ...
              var(yF1pm1 .* yActual, 'omitmissing'));
pwBench = cBench * yF1pm1;
if weightRestriction
    pwBench(pwBench < 0) = 0;
    pwBench(pwBench > 2) = 2;
end
yPTBench = pwBench .* yActual;

alphaAnnualBest  = NaN(2, numSpecifications - 1);
tStatAlphMatBest = NaN(2, numSpecifications - 1);
SRBestAnnual     = NaN(2, numSpecifications - 1);
SRmkt = mean(yActual(dateMask), 'omitmissing') / std(yActual(dateMask), 'omitmissing');
bestInds = {bestRMSEInd, bestAlphInd};

for ii = 1:(numSpecifications - 1)
    yF2_ii = squeeze(yF2AllS9(:, ii, :));   % T x N

    for ax = 1:2
        bInd = bestInds{ax};
        ypt  = NaN(T, 1);
        pw   = NaN(T, 1);                       % portfolio weight series (for adjTerm)
        ypt(1:evalStartInd) = yPTBench(1:evalStartInd);
        pw(1:evalStartInd)  = pwBench(1:evalStartInd);
        for t = evalStartInd+1:T
            yF2col = yF2_ii(:, bInd(t-1, ii));
            cIn = sqrt(var(yActual(1:t-1), 'omitmissing') / ...
                       var(yF2col(1:t-1) .* yActual(1:t-1), 'omitmissing'));
            pwIn = cIn * yF2col(1:t);
            if weightRestriction
                pwIn(pwIn < 0) = 0;
                pwIn(pwIn > 2) = 2;
            end
            ypt(t) = pwIn(t) * yActual(t);
            pw(t)  = pwIn(t);
        end

        % Transaction-cost adjTerm: |pw(t) - drifted pw(t-1)| via the
        % daily geometric self-financing form used in
        % constructPortfolioFactors.m line 72-78. Subtract from ypt
        % before regressing.
        if adjCost > 0
            pwPrev = pw(1:end-1);
            yPrev  = yActual(1:end-1);
            rfPrev = riskFree(1:end-1);
            drift = pwPrev .* exp(yPrev) ./ ...
                    (pwPrev .* exp(yPrev) + (1 - pwPrev) .* exp(rfPrev));
            drift = min(max(drift, 0), 2);
            adjTerm = [NaN; abs(pw(2:end) - drift)];
            ypt = ypt - adjTerm * adjCost;
        end

        % Restrict to dateMask before regression
        yA_d = yActual(dateMask);
        ypt_d = ypt(dateMask);
        s = regstats2Fast(ypt_d, yA_d);
        alphaAnnualBest(ax, ii)  = s.beta(1) * 252 * 100;
        tStatAlphMatBest(ax, ii) = s.hac.t(1);
        SR = sqrt(SRmkt^2 + s.beta(1)^2 / s.mse);
        SRBestAnnual(ax, ii) = SR * sqrt(252);
    end
end
end
