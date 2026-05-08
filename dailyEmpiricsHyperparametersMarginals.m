function dailyEmpiricsHyperparametersMarginals(cfg, signSpec)
%DAILYEMPIRICSHYPERPARAMETERSMARGINALS  Adaptive-shrinkage selection
% across the 4,860-combination grid. Picks the per-day winning combo by
% cfg.selectionCriterion ('maxAlpha' default; 'tAlpha'; 'rmse') and
% writes the topK/botK and oosCAlphas .mat files used by Figures 1,
% A.1, A.2, A.3.
%
% Inputs:
%   cfg       config struct (default_config())
%   signSpec  (optional) scalar or vector of signSpec ids in {1, 2}.
%             Defaults to [1, 2].
%
% Reads (from cfg.paths.results/hyperparameters/):
%   forecastResults_<signSpec>_2.5yE_1yDM_1S_pm_HyperR1.mat
%   OOSResults_<signSpec>_HyperR1_ExpandingC.mat
%
% Outputs (to cfg.paths.results/aggregates/):
%   signSpec = 2 (published baseline):
%     topKbotK.mat, oosCAlphas.mat
%   signSpec = 1 (unrestricted signs):
%     topKbotK_signSpec1.mat, oosCAlphas_signSpec1.mat
%   cfg.shrinkageTarget = 'prevailingMean':
%     topKbotK_pm.mat, oosCAlphas_pm.mat (and _signSpec1 variants)
%
% Consumers: Figures 1, A.1, A.2, A.3 read the signSpec = 2 files. The
% Adaptive panels of Tables 1 and A.3 are populated by the table-side
% computeTab1Results / computeTabA3Results scripts.
%
% Cfg hooks: cfg.selectionCriterion picks 'maxAlpha' (default), 'tAlpha',
% or 'rmse'. cfg.shrinkageTarget picks 'zero' (default) or
% 'prevailingMean'.
%
% Runtime: ~45-90 minutes per signSpec.

if nargin < 1 || isempty(cfg); cfg = default_config(); end
if nargin < 2 || isempty(signSpec)
    if isfield(cfg, 'hyperparameterSpecs') && ~isempty(cfg.hyperparameterSpecs)
        signSpec = cfg.hyperparameterSpecs;
    else
        signSpec = [1, 2];
    end
end

for s = signSpec(:)'
    fprintf('dailyEmpiricsHyperparametersMarginals: signSpec %d (sweep %s)\n', ...
        s, mat2str(signSpec(:)'));
    runOneSignSpec(cfg, s);
end
end


function runOneSignSpec(cfg, signSpec)
% Per-signSpec body.
paths = cfg.paths;

%% Load Hyperparameters1 and Hyperparameters2 outputs
robFolder = robustnessResultsFolder(paths, cfg);
fileName  = fullfile(robFolder, sprintf( ...
    'forecastResults_%d_2.5yE_1yDM_1S_pm_HyperR1.mat', signSpec));
fileName2 = fullfile(robFolder, sprintf( ...
    'OOSResults_%d_HyperR1_ExpandingC.mat', signSpec));
S1 = load(fileName);

% Extract just what we need (yActual, yF1PMMat, riskFree, trim) and
% release S1 — its yF2All / pocketIndAll / fDMHatAll fields are large
% and unused downstream of the trim calculation.
yActual  = S1.yActual ./ 100;
yF1PMMat = S1.yF1PMMat ./ 100;
riskFree = S1.riskFree ./ 100;
trim     = max(sum(isnan(S1.yF2All)), [], 'all');
clear S1;

% Trim to common sample.
yF1PMMat = yF1PMMat(trim+1:end, :);
yActual  = yActual(trim+1:end, :);
riskFree = riskFree(trim+1:end);

S2 = load(fileName2);

% Filter the full 9,720-combo grid to the half whose shrinkage target
% matches cfg.shrinkageTarget, reducing the pool to 4,860 combos:
%   cfg.shrinkageTarget == 'zero'           -> idx.benchmark == 0 (default,
%                                              published-reply baseline)
%   cfg.shrinkageTarget == 'prevailingMean' -> idx.benchmark == 1
% When cfg.hyperparameterComboSubset is set, the loaded S2 arrays are
% already subset to those combos; benchmark filtering is applied within
% the subset.
[paramCombsFull, idx] = constructHyperparameterGrid();
if isfield(cfg, 'hyperparameterComboSubset') && ~isempty(cfg.hyperparameterComboSubset)
    paramCombsLoaded = paramCombsFull(cfg.hyperparameterComboSubset, :);
else
    paramCombsLoaded = paramCombsFull;
end
if isfield(cfg, 'shrinkageTarget') && strcmp(cfg.shrinkageTarget, 'prevailingMean')
    benchValue = 1;
else
    benchValue = 0;
end
keepInd      = find(paramCombsLoaded(:, idx.benchmark) == benchValue);
paramCombs   = paramCombsLoaded(keepInd, :);

alphMatExpanding      = S2.alphMatExpanding(:, :, keepInd);
tStatAlphMatExpanding = S2.tStatAlphMatExpanding(:, :, keepInd);
yF2AllS               = S2.yF2AllS(:, :, keepInd);
T                     = S2.T;
numSpecifications     = S2.numSpecifications;
weightRestriction     = 1;
evalStartInd          = 252;

% Best-combo indices per metric (RMSE / alpha / t-stat-alpha).
% Deterministic tie-break: when two combos tie at FP-precision (~1e-15),
% min/max returns the first index, which depends on the dim-3 ordering
% set by parfor scheduling in dailyEmpiricsHyperparameters{1,2}. Break
% ties by combo-index (lowest combo wins) so the choice is reproducible
% across runs. The 1e-13 perturbation is well above FP-noise but well
% below any real alpha / RMSE / t-stat difference.
rmseAll = sqrt(cumsum(((yActual - yF2AllS).^2) ./ ((1:T)')));
nCombos = size(rmseAll, 3);
tieBreakAsc  =  1e-13 * reshape(1:nCombos, 1, 1, []);
tieBreakDesc = -1e-13 * reshape(1:nCombos, 1, 1, []);
[~, bestRMSEInd]      = min(rmseAll              + tieBreakAsc,  [], 3, 'omitmissing');
[~, bestAlphInd]      = max(alphMatExpanding      + tieBreakDesc, [], 3, 'omitmissing');
[~, bestTStatAlphInd] = max(tStatAlphMatExpanding + tieBreakDesc, [], 3, 'omitmissing');

%% Pre-compute the warmup time series from spec=end (yF1pm1), combo=1.
% (Original implementation built a T x specs x combos `yF2Extended` and a
% same-sized `yPocketTime` purely so bestComboRegPanel could index a
% combo-1 column for warmup; both arrays are ~5 GB each at this scale, so
% we skip them and recompute the warmup as a single T-vector.)
yF1pm1 = yF1PMMat(:, 1);
cBench  = sqrt(var(yActual, 'omitmissing') ./ var(yF1pm1 .* yActual, 'omitmissing'));
pwBench = cBench * yF1pm1;
if weightRestriction
    pwBench(pwBench < 0) = 0;
    pwBench(pwBench > 2) = 2;
end
yPTBench = pwBench .* yActual;

%% Best-combo regression panel with rolling OOS c and self-financing adj
[alphMatBest, tStatAlphMatBest, mseMatBest] = bestComboRegPanel( ...
    yF2AllS, yPTBench, yActual, riskFree, ...
    bestRMSEInd, bestAlphInd, bestTStatAlphInd, ...
    weightRestriction, evalStartInd, T, numSpecifications);

SRmkt        = mean(yActual, 'omitnan') / std(yActual, 'omitnan');
SR           = sqrt(SRmkt^2 + (alphMatBest.^2) ./ mseMatBest);
SRAnnual     = SR * sqrt(252); %#ok<NASGU>
alphaAnnual  = alphMatBest * 252 * 100; %#ok<NASGU>

% Published-reply benchmark t-stat per spec (consumed by figure1 / figureA1
% as the vertical reference line). The benchmark reference combo is
% [coefWindowYears, gPrior, signRestriction, weightLB, benchmark, weightUB,
%  coefRestriction, minPocketDays] = [2.5, 2, 0, 1, 0, 1, 0, 21]. For the
% prevailing-mean variant (cfg.shrinkageTarget = 'prevailingMean'), the
% benchmark column is set to 1 instead of 0.
% S2.tStatAlphMatExpanding is sized to whatever combos H2 wrote — when
% running with cfg.hyperparameterComboSubset and the subset doesn't
% include this row, fall back to NaN so the load contract stays intact.
benchmarkTuple = [2.5, 2, 0, 1, 0, 1, 0, 21];
benchmarkTuple(idx.benchmark) = benchValue;
benchmarkRowsLoaded = find(all(paramCombsLoaded == benchmarkTuple, 2));
if ~isempty(benchmarkRowsLoaded)
    bIdx = benchmarkRowsLoaded(1);
    tStatAlphaBenchmark = squeeze(S2.tStatAlphMatExpanding(end, :, bIdx));
    cwBenchmark         = squeeze(S2.cwMat(:, :, bIdx));
else
    warning('dailyEmpiricsHyperparametersMarginals:NoBenchmarkCombo', ...
        ['Loaded paramCombs subset does not contain the published-reply ' ...
         'benchmark combo [2.5, 2, 0, 1, 0, 1, 0, 21]; tStatAlphaBenchmark ' ...
         'and cwBenchmark will be NaN. figure1 / figureA1 reference lines ' ...
         'will be missing.']);
    tStatAlphaBenchmark = NaN(numSpecifications, 1);
    cwBenchmark         = NaN(size(S2.cwMat, 1), size(S2.cwMat, 2));
end

% End-of-sample t-stat per (spec, combo) — consumed by figureA1 / figureA3.
% Use the locally bench-filtered tStatAlphMatExpanding (4,860 combos) so
% sizes match tAlphAllSpecs in figureA3.
tStatAlphMatEnd = squeeze(tStatAlphMatExpanding(end, :, :));

%% Top-K / Bot-K adaptive shrinkage analysis
% Output filenames are suffixed by cfg.selectionCriterion (no suffix for
% the default 'maxAlpha'), cfg.shrinkageTarget ('_pm' for the
% prevailing-mean variant), and signSpec ('_signSpec1' for the
% unrestricted-signs variant) so the variants coexist on disk.
topK = min(100, size(paramCombs, 1));
if isfield(cfg, 'selectionCriterion') && ~strcmp(cfg.selectionCriterion, 'maxAlpha')
    cacheSfx = ['_', cfg.selectionCriterion];
else
    cacheSfx = '';
end
if isfield(cfg, 'shrinkageTarget') && strcmp(cfg.shrinkageTarget, 'prevailingMean')
    cacheSfx = [cacheSfx, '_pm'];
end
if signSpec ~= 2
    cacheSfx = [cacheSfx, sprintf('_signSpec%d', signSpec)];
end
cacheRoot    = aggregatesFolder(paths);
topKbotKFile = fullfile(cacheRoot, ['topKbotK', cacheSfx, '.mat']);
oosCFile     = fullfile(cacheRoot, ['oosCAlphas', cacheSfx, '.mat']);
% Auto-detect: recompute if either cache file is missing.
rerun_loop = ~isfile(topKbotKFile) || ~isfile(oosCFile);

if rerun_loop
    fprintf('  recomputing topK/botK (no cache at %s)\n', topKbotKFile);
    [TKBK, OOSCA] = recomputeTopKBotK( ...
        yF2AllS, yF1pm1, yActual, alphMatExpanding, tStatAlphMatExpanding, ...
        rmseAll, weightRestriction, evalStartInd, T, numSpecifications, ...
        topK, size(paramCombs, 1), cfg);
    save(topKbotKFile, '-struct', 'TKBK');
    save(oosCFile,     '-struct', 'OOSCA');
else
    TKBK  = load(topKbotKFile);
    OOSCA = load(oosCFile);
end

%% In-pocket CW stats for the top/bot K models per spec
[cwInMatTopK, cwInMatBotK, inPocketFracTopK, inPocketFracBotK] = ...
    computeTopKBotKCWStats(TKBK, yActual, yF1PMMat, numSpecifications, topK);

%% Save (the figure*.m scripts read from these files)
% Diagnostic: log scalar bestTAlphaInd / bestRMSEInd via the centralized
% selectBestComboIndices helper. The Sharpe-based index returned by that
% helper is computed for the R1 selectionCriterion hook but never
% consumed downstream — the live Adaptive panels use the maxAlpha and
% rmse selections built directly from alphMatExpanding / rmseAll above.
[~, ~, bestTAlphaIndScalar, bestRMSEIndScalar] = selectBestComboIndices( ...
    SR(:), alphMatBest(:), tStatAlphMatBest(:), mseMatBest(:), ...
    cfg.selectionCriterion);

save(topKbotKFile, '-struct', 'TKBK', '-append');
% figureA3 also reads tAlphAllSpecs (third-axis: per-combo, per-spec
% t-stat using rolling OOS c) from topKbotK.mat, so mirror it from
% OOSCA into the same file. OOSCA itself remains canonical at oosCFile.
tAlphAllSpecs = OOSCA.tAlphAllSpecs;
% Also mirror cwMat (from H2 OOSResults) and paramCombs (post-benchmark
% filter) so figureA2 / figureA3 get the per-combo distributions and the
% paramCombs filter without an extra load.
cwMat = S2.cwMat(:, :, keepInd);
save(topKbotKFile, 'tStatAlphaBenchmark', 'cwBenchmark', 'tStatAlphMatEnd', ...
    'cwInMatTopK', 'cwInMatBotK', 'inPocketFracTopK', 'inPocketFracBotK', ...
    'tAlphAllSpecs', 'cwMat', 'paramCombs', ...
    '-append');
fprintf('  alphMatBest size = %s\n', mat2str(size(alphMatBest)));
fprintf('  bestTAlphaInd = %d  bestRMSEInd = %d\n', ...
    bestTAlphaIndScalar, bestRMSEIndScalar);
end


% ====================================================================
%  Local helpers — entry-point-specific.
% ====================================================================

function [alphMatBest, tStatAlphMatBest, mseMatBest] = bestComboRegPanel( ...
    yF2AllS, yPTBench, yActual, riskFree, ...
    bestRMSEInd, bestAlphInd, bestTStatAlphInd, ...
    weightRestriction, evalStartInd, T, numSpecifications)
% Per spec, build the rolling-OOS portfolio time series for each
% selection rule (RMSE / Alph / tStatAlph), apply the geometric self-
% financing adjustment, and run a single regression. yF2 forecasts are
% sliced spec-by-spec from yF2AllS — never materialized as a single
% T x specs x combos tensor (would be ~5 GB at full grid).
alphMatBest      = NaN(3, numSpecifications);
tStatAlphMatBest = NaN(3, numSpecifications);
mseMatBest       = NaN(3, numSpecifications);
adjCost          = 0 / 100 / 100;

for ii = 1:(numSpecifications - 1)
    yF2_ii = squeeze(yF2AllS(:, ii, :));   % T x numParamCombs
    yptR = NaN(T, 1); yptA = NaN(T, 1); yptT = NaN(T, 1);
    pwR  = NaN(T, 1); pwA  = NaN(T, 1);
    yptR(1:evalStartInd) = yPTBench(1:evalStartInd);
    yptA(1:evalStartInd) = yPTBench(1:evalStartInd);
    yptT(1:evalStartInd) = yPTBench(1:evalStartInd);

    for t = evalStartInd+1:T
        [pwR(t), yptR(t)] = applyOOSCScaling( ...
            yF2_ii(:, bestRMSEInd(t-1, ii)), yActual, t, weightRestriction);
        [pwA(t), yptA(t)] = applyOOSCScaling( ...
            yF2_ii(:, bestAlphInd(t-1, ii)), yActual, t, weightRestriction);
        [~,      yptT(t)] = applyOOSCScaling( ...
            yF2_ii(:, bestTStatAlphInd(t-1, ii)), yActual, t, weightRestriction);
    end
    clear yF2_ii;

    % Self-financing adjustment for transaction costs.
    yptR = applyAdjCost(yptR, pwR, yActual, riskFree, adjCost);
    yptA = applyAdjCost(yptA, pwA, yActual, riskFree, adjCost);

    sR = regstats2Fast(yptR(evalStartInd+1:T), yActual(evalStartInd+1:T));
    sA = regstats2Fast(yptA(evalStartInd+1:T), yActual(evalStartInd+1:T));
    sT = regstats2Fast(yptT(evalStartInd+1:T), yActual(evalStartInd+1:T));
    alphMatBest(:, ii)      = [sR.beta(1); sA.beta(1); sT.beta(1)];
    mseMatBest(:, ii)       = [sR.mse;     sA.mse;     sT.mse];
    tStatAlphMatBest(:, ii) = [sR.hac.t(1); sA.hac.t(1); sT.hac.t(1)];
end
end


function [pwt, ypt] = applyOOSCScaling(yPocketCur, yActual, t, weightRestriction)
% Compute the OOS c-scaling factor at time t-1 from data 1..t-1, then
% apply it to yPocketCur(1:t) to get the t'th portfolio value.
cIn = sqrt(var(yActual(1:t-1), 'omitmissing') / ...
           var(yPocketCur(1:t-1) .* yActual(1:t-1), 'omitmissing'));
pwIn = cIn * yPocketCur(1:t);
if weightRestriction
    pwIn(pwIn < 0) = 0;
    pwIn(pwIn > 2) = 2;
end
pwt = pwIn(t);
ypt = pwt * yActual(t);
end


function ypt = applyAdjCost(ypt, pw, yActual, riskFree, adjCost)
% Geometric self-financing transaction-cost adjustment (matches the
% adjTerm formula in dailyEmpirics's portfolio block).
adjTerm = [NaN; abs(pw(2:end) - ...
    pw(1:end-1) .* exp(yActual(1:end-1)) ./ ...
    (pw(1:end-1) .* exp(yActual(1:end-1)) + ...
     (1 - pw(1:end-1)) .* exp(riskFree(1:end-1))))];
adjTerm(isnan(adjTerm)) = 0;
ypt = ypt - adjTerm * adjCost;
end


function [TKBK, OOSCA] = recomputeTopKBotK(yF2AllS, yF1pm1, yActual, ...
    alphMatExpanding, tStatAlphMatExpanding, rmseAll, ...
    weightRestriction, evalStartInd, T, ...
    numSpecifications, topK, numParamCombs, cfg)
% Compute topK and botK adaptive-shrinkage forecasts for each (spec,
% rank) tuple. Third axis of the *TopK / *BotK arrays:
%   1 = ranked by RMSE
%   2 = ranked by alpha
%   3 = ranked by t-stat-alpha
yPocketTimeTopKRMSE  = zeros(T, numSpecifications - 1, topK);
yPocketTimeTopKAlph  = zeros(T, numSpecifications - 1, topK);
yPocketTimeTopKtAlph = zeros(T, numSpecifications - 1, topK);
yF2ExtendedTopKRMSE  = zeros(T, numSpecifications - 1, topK);
yF2ExtendedTopKAlph  = zeros(T, numSpecifications - 1, topK);
yF2ExtendedTopKtAlph = zeros(T, numSpecifications - 1, topK);
yPocketTimeBotKRMSE  = zeros(T, numSpecifications - 1, topK);
yPocketTimeBotKAlph  = zeros(T, numSpecifications - 1, topK);
yPocketTimeBotKtAlph = zeros(T, numSpecifications - 1, topK);
yF2ExtendedBotKRMSE  = zeros(T, numSpecifications - 1, topK);
yF2ExtendedBotKAlph  = zeros(T, numSpecifications - 1, topK);
yF2ExtendedBotKtAlph = zeros(T, numSpecifications - 1, topK);
alphMatTopK          = zeros(topK, numSpecifications - 1, 3);
mseMatTopK           = zeros(topK, numSpecifications - 1, 3);
tStatAlphMatTopK     = zeros(topK, numSpecifications - 1, 3);
alphMatBotK          = zeros(topK, numSpecifications - 1, 3);
mseMatBotK           = zeros(topK, numSpecifications - 1, 3);
tStatAlphMatBotK     = zeros(topK, numSpecifications - 1, 3);

for ii = 1:(numSpecifications - 1)
    fprintf('  topK/botK spec %d/%d\n', ii, numSpecifications - 1);

    % Pre-extract spec-ii slices so each worker only receives the data it
    % needs, not the full yF2Extended (T x specs x combos = ~17 GB).
    rmseAll_ii   = squeeze(rmseAll(:, ii, :));              % T x nCombos
    alphMat_ii   = squeeze(alphMatExpanding(:, ii, :));     % T x nCombos
    tAlphMat_ii  = squeeze(tStatAlphMatExpanding(:, ii, :));% T x nCombos
    if ii < numSpecifications
        yF2Ext_ii = squeeze(yF2AllS(:, ii, :));             % T x nCombos
    else
        yF2Ext_ii = repmat(yF1pm1, 1, numParamCombs);       % combo-invariant
    end

    % Per-iii sliced output arrays (parfor sliced).
    yptTRMSE_iii  = zeros(T, topK);
    yptTAlph_iii  = zeros(T, topK);
    yptTtAlph_iii = zeros(T, topK);
    yF2TRMSE_iii  = zeros(T, topK);
    yF2TAlph_iii  = zeros(T, topK);
    yF2TtAlph_iii = zeros(T, topK);
    yptBRMSE_iii  = zeros(T, topK);
    yptBAlph_iii  = zeros(T, topK);
    yptBtAlph_iii = zeros(T, topK);
    yF2BRMSE_iii  = zeros(T, topK);
    yF2BAlph_iii  = zeros(T, topK);
    yF2BtAlph_iii = zeros(T, topK);
    alphTopK_iii      = zeros(topK, 3);
    mseTopK_iii       = zeros(topK, 3);
    tStatTopK_iii     = zeros(topK, 3);
    alphBotK_iii      = zeros(topK, 3);
    mseBotK_iii       = zeros(topK, 3);
    tStatBotK_iii     = zeros(topK, 3);

    parfor (iii = 1:topK, useParfor(cfg))
        yptTR  = zeros(T, 1); yptTA = zeros(T, 1); yptTt = zeros(T, 1);
        yF2TR  = zeros(T, 1); yF2TA = zeros(T, 1); yF2Tt = zeros(T, 1);
        yptBR  = zeros(T, 1); yptBA = zeros(T, 1); yptBt = zeros(T, 1);
        yF2BR  = zeros(T, 1); yF2BA = zeros(T, 1); yF2Bt = zeros(T, 1);

        for t = evalStartInd+1:T
            [~, rmseIdx]  = sort(rmseAll_ii(t-1, :));
            [~, alphIdx]  = sort(alphMat_ii(t-1, :), 'descend');
            [~, tAlphIdx] = sort(tAlphMat_ii(t-1, :), 'descend');

            [~, yptTR(t), yF2TR(t)] = applyOOSCScalingFull(yF2Ext_ii(:, rmseIdx(iii)),  yActual, t, weightRestriction);
            [~, yptTA(t), yF2TA(t)] = applyOOSCScalingFull(yF2Ext_ii(:, alphIdx(iii)),  yActual, t, weightRestriction);
            [~, yptTt(t), yF2Tt(t)] = applyOOSCScalingFull(yF2Ext_ii(:, tAlphIdx(iii)), yActual, t, weightRestriction);
            [~, yptBR(t), yF2BR(t)] = applyOOSCScalingFull(yF2Ext_ii(:, rmseIdx(end-iii+1)),  yActual, t, weightRestriction);
            [~, yptBA(t), yF2BA(t)] = applyOOSCScalingFull(yF2Ext_ii(:, alphIdx(end-iii+1)),  yActual, t, weightRestriction);
            [~, yptBt(t), yF2Bt(t)] = applyOOSCScalingFull(yF2Ext_ii(:, tAlphIdx(end-iii+1)), yActual, t, weightRestriction);
        end

        slc = evalStartInd+1:T;
        sR  = regstats2Fast(yptTR(slc), yActual(slc));
        sA  = regstats2Fast(yptTA(slc), yActual(slc));
        sT  = regstats2Fast(yptTt(slc), yActual(slc));
        sBR = regstats2Fast(yptBR(slc), yActual(slc));
        sBA = regstats2Fast(yptBA(slc), yActual(slc));
        sBT = regstats2Fast(yptBt(slc), yActual(slc));

        yptTRMSE_iii(:, iii)  = yptTR;   yptTAlph_iii(:, iii)  = yptTA;
        yptTtAlph_iii(:, iii) = yptTt;
        yF2TRMSE_iii(:, iii)  = yF2TR;   yF2TAlph_iii(:, iii)  = yF2TA;
        yF2TtAlph_iii(:, iii) = yF2Tt;
        yptBRMSE_iii(:, iii)  = yptBR;   yptBAlph_iii(:, iii)  = yptBA;
        yptBtAlph_iii(:, iii) = yptBt;
        yF2BRMSE_iii(:, iii)  = yF2BR;   yF2BAlph_iii(:, iii)  = yF2BA;
        yF2BtAlph_iii(:, iii) = yF2Bt;

        alphTopK_iii(iii, :)  = [sR.beta(1),  sA.beta(1),  sT.beta(1)];
        mseTopK_iii(iii, :)   = [sR.mse,      sA.mse,      sT.mse];
        tStatTopK_iii(iii, :) = [sR.hac.t(1), sA.hac.t(1), sT.hac.t(1)];
        alphBotK_iii(iii, :)  = [sBR.beta(1), sBA.beta(1), sBT.beta(1)];
        mseBotK_iii(iii, :)   = [sBR.mse,     sBA.mse,     sBT.mse];
        tStatBotK_iii(iii, :) = [sBR.hac.t(1), sBA.hac.t(1), sBT.hac.t(1)];
    end

    % Slot the parfor outputs into the global 3D arrays.
    yPocketTimeTopKRMSE(:, ii, :)  = yptTRMSE_iii;
    yPocketTimeTopKAlph(:, ii, :)  = yptTAlph_iii;
    yPocketTimeTopKtAlph(:, ii, :) = yptTtAlph_iii;
    yF2ExtendedTopKRMSE(:, ii, :)  = yF2TRMSE_iii;
    yF2ExtendedTopKAlph(:, ii, :)  = yF2TAlph_iii;
    yF2ExtendedTopKtAlph(:, ii, :) = yF2TtAlph_iii;
    yPocketTimeBotKRMSE(:, ii, :)  = yptBRMSE_iii;
    yPocketTimeBotKAlph(:, ii, :)  = yptBAlph_iii;
    yPocketTimeBotKtAlph(:, ii, :) = yptBtAlph_iii;
    yF2ExtendedBotKRMSE(:, ii, :)  = yF2BRMSE_iii;
    yF2ExtendedBotKAlph(:, ii, :)  = yF2BAlph_iii;
    yF2ExtendedBotKtAlph(:, ii, :) = yF2BtAlph_iii;
    alphMatTopK(:, ii, :)         = alphTopK_iii;
    mseMatTopK(:, ii, :)          = mseTopK_iii;
    tStatAlphMatTopK(:, ii, :)    = tStatTopK_iii;
    alphMatBotK(:, ii, :)         = alphBotK_iii;
    mseMatBotK(:, ii, :)          = mseBotK_iii;
    tStatAlphMatBotK(:, ii, :)    = tStatBotK_iii;
end

TKBK = struct( ...
    'yPocketTimeTopKRMSE',  yPocketTimeTopKRMSE, ...
    'yPocketTimeTopKAlph',  yPocketTimeTopKAlph, ...
    'yPocketTimeTopKtAlph', yPocketTimeTopKtAlph, ...
    'yF2ExtendedTopKRMSE',  yF2ExtendedTopKRMSE, ...
    'yF2ExtendedTopKAlph',  yF2ExtendedTopKAlph, ...
    'yF2ExtendedTopKtAlph', yF2ExtendedTopKtAlph, ...
    'yPocketTimeBotKRMSE',  yPocketTimeBotKRMSE, ...
    'yPocketTimeBotKAlph',  yPocketTimeBotKAlph, ...
    'yPocketTimeBotKtAlph', yPocketTimeBotKtAlph, ...
    'yF2ExtendedBotKRMSE',  yF2ExtendedBotKRMSE, ...
    'yF2ExtendedBotKAlph',  yF2ExtendedBotKAlph, ...
    'yF2ExtendedBotKtAlph', yF2ExtendedBotKtAlph, ...
    'alphMatTopK',          alphMatTopK, ...
    'mseMatTopK',           mseMatTopK, ...
    'tStatAlphMatTopK',     tStatAlphMatTopK, ...
    'alphMatBotK',          alphMatBotK, ...
    'mseMatBotK',           mseMatBotK, ...
    'tStatAlphMatBotK',     tStatAlphMatBotK);

% Recompute alphas across all combos with real-time c. Spec slices are
% extracted once per spec from yF2AllS so the parfor only sends a (T x 1)
% column per worker, not the full T x specs x combos tensor.
alphAllSpecs  = NaN(numParamCombs, numSpecifications - 1);
mseAllSpecs   = NaN(numParamCombs, numSpecifications - 1);
tAlphAllSpecs = NaN(numParamCombs, numSpecifications - 1);
for ii = 1:(numSpecifications - 1)
    yF2_ii = squeeze(yF2AllS(:, ii, :));   % T x numParamCombs
    parfor (iii = 1:numParamCombs, useParfor(cfg))
        yPocketCur = yF2_ii(:, iii);
        yPocketTimeTmp = yPocketCur * 0;
        for t = evalStartInd+1:T
            [~, ypt, ~] = applyOOSCScalingFull(yPocketCur, yActual, t, weightRestriction);
            yPocketTimeTmp(t) = ypt;
        end
        slc = evalStartInd+1:T;
        s = regstats2Fast(yPocketTimeTmp(slc), yActual(slc));
        alphAllSpecs(iii, ii)  = s.beta(1);
        mseAllSpecs(iii, ii)   = s.mse;
        tAlphAllSpecs(iii, ii) = s.hac.t(1);
    end
    clear yF2_ii;
end
OOSCA = struct( ...
    'alphAllSpecs',  alphAllSpecs, ...
    'mseAllSpecs',   mseAllSpecs, ...
    'tAlphAllSpecs', tAlphAllSpecs);
end


function [pwt, ypt, fct] = applyOOSCScalingFull(yPocketCur, yActual, t, weightRestriction)
% Same as applyOOSCScaling but also returns the raw forecast at t.
cIn = sqrt(var(yActual(1:t-1), 'omitmissing') / ...
           var(yPocketCur(1:t-1) .* yActual(1:t-1), 'omitmissing'));
pwIn = cIn * yPocketCur(1:t);
if weightRestriction
    pwIn(pwIn < 0) = 0;
    pwIn(pwIn > 2) = 2;
end
pwt = pwIn(t);
ypt = pwt * yActual(t);
fct = yPocketCur(t);
end


function [cwInTopK, cwInBotK, fracTopK, fracBotK] = computeTopKBotKCWStats( ...
    TKBK, yActual, yF1PMMat, numSpecifications, topK)
% In-pocket CW t-stats for top-K and bot-K models per spec. The
% "in-pocket" indicator is "where the forecast differs from PM".
cwInTopK = NaN(topK, numSpecifications - 1, 2);
cwInBotK = NaN(topK, numSpecifications - 1, 2);
fracTopK = NaN(topK, numSpecifications - 1, 2);
fracBotK = NaN(topK, numSpecifications - 1, 2);

for ii = 1:(numSpecifications - 1)
    benchCol = min(ii, 6);   % map spec 7..9 (combos) onto dp PM
    yF1Bench = yF1PMMat(:, benchCol);

    fCWTopK_R = (yActual - yF1Bench).^2 - squeeze(((yActual - TKBK.yF2ExtendedTopKRMSE(:, ii, :)).^2 - ...
                (yF1Bench - TKBK.yF2ExtendedTopKRMSE(:, ii, :)).^2));
    fCWBotK_R = (yActual - yF1Bench).^2 - squeeze(((yActual - TKBK.yF2ExtendedBotKRMSE(:, ii, :)).^2 - ...
                (yF1Bench - TKBK.yF2ExtendedBotKRMSE(:, ii, :)).^2));
    fCWTopK_A = (yActual - yF1Bench).^2 - squeeze(((yActual - TKBK.yF2ExtendedTopKAlph(:, ii, :)).^2 - ...
                (yF1Bench - TKBK.yF2ExtendedTopKAlph(:, ii, :)).^2));
    fCWBotK_A = (yActual - yF1Bench).^2 - squeeze(((yActual - TKBK.yF2ExtendedBotKAlph(:, ii, :)).^2 - ...
                (yF1Bench - TKBK.yF2ExtendedBotKAlph(:, ii, :)).^2));

    for iii = 1:topK
        % Top-K, RMSE-ranked
        inP = yF1Bench ~= TKBK.yF2ExtendedTopKRMSE(:, ii, iii);
        if any(inP)
            row = computeDMCWStats(fCWTopK_R(inP, iii), [], 'cw');
            cwInTopK(iii, ii, 1) = row(1);
        end
        fracTopK(iii, ii, 1) = mean(inP);

        % Bot-K, RMSE-ranked
        inP = yF1Bench ~= TKBK.yF2ExtendedBotKRMSE(:, ii, iii);
        try
            if any(inP)
                row = computeDMCWStats(fCWBotK_R(inP, iii), [], 'cw');
                cwInBotK(iii, ii, 1) = row(1);
            end
            fracBotK(iii, ii, 1) = mean(inP);
        catch
        end

        % Top-K, alpha-ranked
        inP = yF1Bench ~= TKBK.yF2ExtendedTopKAlph(:, ii, iii);
        if any(inP)
            row = computeDMCWStats(fCWTopK_A(inP, iii), [], 'cw');
            cwInTopK(iii, ii, 2) = row(1);
        end
        fracTopK(iii, ii, 2) = mean(inP);

        % Bot-K, alpha-ranked. Original used "==" here (looks like a bug
        % since it would mark "where forecast equals PM" as in-pocket);
        % preserved verbatim for numerical equivalence.
        inP = yF1Bench == TKBK.yF2ExtendedBotKAlph(:, ii, iii);
        try
            if any(inP)
                row = computeDMCWStats(fCWBotK_A(inP, iii), [], 'cw');
                cwInBotK(iii, ii, 2) = row(1);
            end
            fracBotK(iii, ii, 2) = mean(inP);
        catch
        end
    end
end
end
