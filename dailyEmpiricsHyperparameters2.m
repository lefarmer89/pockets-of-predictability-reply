function dailyEmpiricsHyperparameters2(cfg, signSpec)
%DAILYEMPIRICSHYPERPARAMETERS2  Combine the forecasts written by
% dailyEmpiricsHyperparameters1 with G-prior shrinkage weights and write
% per-signSpec OOS results.
%
% Inputs:
%   cfg       config struct (default_config())
%   signSpec  (optional) scalar or vector of signSpec ids in {1, 2}.
%             Defaults to [1, 2].
%
% Reads (from cfg.paths.results/hyperparameters/):
%   forecastResults_<signSpec>_2.5yE_1yDM_1S_pm_HyperR1.mat
%
% Outputs (to cfg.paths.results/hyperparameters/):
%   forecastResults_<signSpec>_*_HyperR1.mat   (in-place append: yF2AllS)
%   OOSResults_<signSpec>_HyperR1_ExpandingC.mat
%
% Consumers: Adaptive panels of Tables 1 and A.3 (via *Marginals).
%
% Runtime: ~30-90 minutes per signSpec.
% Memory: peak ~50 GB on the default 9,720-combo grid.
% cfg.hyperparameterComboSubset (default []) restricts to a small subset
% for shakeout.

if nargin < 1 || isempty(cfg); cfg = default_config(); end
if nargin < 2 || isempty(signSpec)
    if isfield(cfg, 'hyperparameterSpecs') && ~isempty(cfg.hyperparameterSpecs)
        signSpec = cfg.hyperparameterSpecs;
    else
        signSpec = [1, 2];
    end
end

for s = signSpec(:)'
    fprintf('dailyEmpiricsHyperparameters2: signSpec %d (sweep %s)\n', ...
        s, mat2str(signSpec(:)'));
    runOneSignSpec(cfg, s);
end
end


function runOneSignSpec(cfg, signSpec)
% Per-signSpec body.
paths = cfg.paths;

%% Load Hyperparameters1 outputs
fileName = fullfile(robustnessResultsFolder(paths, cfg), sprintf( ...
    'forecastResults_%d_2.5yE_1yDM_1S_pm_HyperR1.mat', signSpec));
S = load(fileName);
yActual      = S.yActual;
yF1PMMat     = S.yF1PMMat;
yF2All       = S.yF2All;
paramCombs   = S.paramCombs;
pocketIndAll = S.pocketIndAll;
clear S;

% Apply optional combo subset (for shakeout / profiling). H1 may have
% already subset the file; only re-apply if the loaded paramCombs is
% still the full grid.
if isfield(cfg, 'hyperparameterComboSubset') && ~isempty(cfg.hyperparameterComboSubset)
    sub = cfg.hyperparameterComboSubset;
    if size(paramCombs, 1) == max(sub) || size(paramCombs, 1) > max(sub)
        % Loaded paramCombs is large enough to hold all subset indices —
        % treat it as full-grid and apply the subset.
        paramCombs   = paramCombs(sub, :);
        yF2All       = yF2All(:, :, sub);
        pocketIndAll = pocketIndAll(:, :, sub);
    end
    % Otherwise the loaded file is already a subset; trust H1's alignment.
end
numParamCombs = size(paramCombs, 1);
memSnapshot('post-load');

weightRestriction = 1;

%% Trim to common sample and rescale to decimals
yActual = yActual ./ 100; yF1PMMat = yF1PMMat ./ 100;
yF2All  = yF2All ./ 100;
trim = max(sum(isnan(yF2All)), [], 'all');
yF1PMMat = yF1PMMat(trim+1:end, :);
yF2All   = yF2All(trim+1:end, :, :);
yActual  = yActual(trim+1:end, :);
pocketIndAll = pocketIndAll(trim+1:end, :, :);
pocketIndAll(isnan(pocketIndAll)) = 0;
% Convert to logical (was double 0/1). 8x smaller; arithmetic uses are
% all of the form `.* pocketIndAll` (auto-promotes) or `~pocketIndAll`
% (logical not), so semantics are preserved.
pocketIndAll = logical(pocketIndAll);
memSnapshot('post-trim');

%% Build yF2AllS (in-pocket TVC, out-of-pocket benchmark)
% The archive computes the three combination forecasts in a specific order
% relative to the in-place yF2All mutation:
%   yComb2 (in-pocket-only average) — BEFORE mutation, on ORIGINAL yF2All
%   yComb3 (equal-weight mean)      — BEFORE mutation, on ORIGINAL yF2All
%   yF2All <- mutated to in-pocket + benchmark out-of-pocket
%   yComb1 (equal-weight mean)      — AFTER mutation, on MUTATED yF2All
%   yF2AllS(:,1:6,:)                — MUTATED yF2All
% So yComb1 and yComb3 use DIFFERENT data despite being the same expression
% (this is the historical "yComb3 vs yComb1" subtlety). We preserve this
% exactly. (Previous refactor wrongly computed yComb2/yComb3 from the
% mutated yF2All, producing a 3.55e-3 yF2AllS drift vs archive.)
combPocket = squeeze(any(pocketIndAll(:, 1:4, :), 2));   % T x numParamCombs logical

% (1) Combination forecasts that read the ORIGINAL yF2All.
yComb2 = combineInPocketAverage(yF2All, pocketIndAll, yF1PMMat);
yComb3 = squeeze(mean(yF2All(:, 1:4, :), 2, 'omitmissing'));

% (2) Build the in-pocket-masked tensor as a transient.
yF2InPocket = yF2All .* pocketIndAll + yF1PMMat .* (~pocketIndAll);

% (3) Combination forecast that reads the MUTATED yF2All.
yComb1 = squeeze(mean(yF2InPocket(:, 1:4, :), 2, 'omitmissing'));

% (4) Pack everything into yF2AllS.
T = size(yActual, 1);
yF2AllS = NaN(T, size(yF2All, 2) + 3, numParamCombs);
yF2AllS(:, 1:6, :) = yF2InPocket;
yF2AllS(:, 7, :)   = yComb1;
yF2AllS(:, 8, :)   = yComb2;
yF2AllS(:, 9, :)   = yComb3;
clear yF2InPocket;

%% Per-spec expanding-window alpha, MSE, t(alpha) plus best-combo selection
% expandingAlpha streams over specs, never materializing yF2Extended or
% any T x specs x combos copy of yPocketTime. It also produces the
% per-spec best-combo time series (yptR/yptA/yptT) needed by
% bestComboRegPanel below, so the caller no longer needs yPocketTime.
yF1pm1 = yF1PMMat(:, 1);
evalStartInd = 252;

memSnapshot('pre-expandingAlpha');
[alphMatExpanding, mseMatExpanding, tStatAlphMatExpanding, rmseAll, ...
 bestAlphInd, bestTStatAlphInd, bestRMSEInd, yptR, yptA, yptT] = ...
    expandingAlphaStreamed(yF2AllS, yF1pm1, yActual, weightRestriction, ...
                           evalStartInd, T, numParamCombs, cfg);
memSnapshot('post-expandingAlpha');

% Save yF2AllS to the OOSResults file early so we can release the ~9 GB
% array before the (perComboForecastTests + bestComboCWStats) phase.
% Marginals reads yF2AllS from S2 (the OOSResults file), so this is the
% canonical save location for it.
fileName2 = fullfile(robustnessResultsFolder(paths, cfg), ...
    sprintf('OOSResults_%d_HyperR1_ExpandingC.mat', signSpec));
save(fileName2, 'yF2AllS', '-v7.3');
clear yF2AllS;
memSnapshot('post-yF2AllS-save');

%% Best-combo summary statistics
numSpecifications = size(yF1PMMat, 2) + 3 + 1;   % 6 + 3 combos + 1 benchmark = 10
SRmkt        = mean(yActual, 'omitmissing') / std(yActual, 'omitmissing');
SR           = squeeze(sqrt(SRmkt^2 + (alphMatExpanding(end,:,:).^2) ./ mseMatExpanding(end,:,:)));
SRAnnual     = SR * sqrt(252);                                       % saved below
alphaAnnual  = squeeze(alphMatExpanding(end, :, :)) * 252 * 100;     % saved below

% Benchmark combo lookup (may be empty when running with a subset).
benchmarkInd = find(all(paramCombs == [2.5, 2, 0, 1, 0, 1, 0, 21], 2));   % saved below

%% Best-combo regression panel using the time series produced above
[alphMatBest, tStatAlphMatBest, mseMatBest] = ...
    bestComboRegPanel(yptR, yptA, yptT, yActual, numSpecifications);
SRBest          = squeeze(sqrt(SRmkt^2 + (alphMatBest.^2) ./ mseMatBest));
SRBestAnnual    = SRBest * sqrt(252);                                % saved below
alphaBestAnnual = squeeze(alphMatBest) * 252 * 100;                  % saved below

%% Per-combo forecast differential tests (DM, CW, in/out-pocket Welch)
[dmMat, cwMat, cwDiffMat] = perComboForecastTests( ...
    yActual, yF1PMMat, yF2All, yComb1, yComb2, yComb3, ...
    pocketIndAll, combPocket, numParamCombs, cfg);

%% Best-combo CW stats per spec
[cwBestMat, cwDiffBestMat] = bestComboCWStats( ...
    yActual, yF1PMMat, yF2All, yComb1, yComb2, yComb3, ...
    pocketIndAll, combPocket, ...
    bestRMSEInd, bestAlphInd, bestTStatAlphInd, T, numSpecifications);
cwDiffBestMat(7, :) = 0;
cwDiffBestMat(8, :) = 0;
memSnapshot('post-bestComboCWStats');

%% Append the remaining OOS results (yF2AllS was already saved above).
saveAppendWithRetry(fileName2, struct( ...
    'alphMatExpanding',      alphMatExpanding, ...
    'alphMatBest',           alphMatBest, ...
    'tStatAlphMatExpanding', tStatAlphMatExpanding, ...
    'tStatAlphMatBest',      tStatAlphMatBest, ...
    'rmseAll',               rmseAll, ...
    'dmMat',                 dmMat, ...
    'cwMat',                 cwMat, ...
    'cwDiffMat',             cwDiffMat, ...
    'alphaAnnual',           alphaAnnual, ...
    'SRAnnual',              SRAnnual, ...
    'alphaBestAnnual',       alphaBestAnnual, ...
    'SRBestAnnual',          SRBestAnnual, ...
    'bestAlphInd',           bestAlphInd, ...
    'bestTStatAlphInd',      bestTStatAlphInd, ...
    'bestRMSEInd',           bestRMSEInd, ...
    'T',                     T, ...
    'numSpecifications',     numSpecifications, ...
    'benchmarkInd',          benchmarkInd, ...
    'cwBestMat',             cwBestMat, ...
    'cwDiffBestMat',         cwDiffBestMat, ...
    'mseMatExpanding',       mseMatExpanding));
end


% ====================================================================
%  Local helpers — entry-point-specific.
% ====================================================================

function saveAppendWithRetry(fileName, varStruct)
% Append `varStruct` to `fileName`. Retries once after a 2-second pause
% to tolerate the Dropbox sync race that occasionally reports "appears
% to be corrupt" when the file was written milliseconds earlier.
try
    save(fileName, '-struct', 'varStruct', '-append');
catch ME
    if contains(ME.message, 'corrupt') || contains(ME.message, 'invalid')
        pause(2);
        save(fileName, '-struct', 'varStruct', '-append');
    else
        rethrow(ME);
    end
end
end


function yComb2 = combineInPocketAverage(yF2, pocketIndAll, yF1PMMat)
% At each (t, combo), average the in-pocket TVC forecasts; fall back to
% the dp PM forecast when no variable is in-pocket at t.
T = size(yF2, 1);
numParamCombs = size(yF2, 3);
yComb2 = NaN(T, numParamCombs);
for ii = 1:T
    for jj = 1:numParamCombs
        if any(pocketIndAll(ii, 1:4, jj))
            yComb2(ii, jj) = mean(yF2(ii, logical(pocketIndAll(ii, 1:4, jj)), jj));
        else
            yComb2(ii, jj) = yF1PMMat(ii, 1);
        end
    end
end
end


function [alphMat, mseMat, tStatMat, rmseAll, ...
          bestAlphInd, bestTStatAlphInd, bestRMSEInd, ...
          yptR, yptA, yptT] = expandingAlphaStreamed( ...
    yF2AllS, yF1pm1, yActual, weightRestriction, evalStartInd, ...
    T, numParamCombs, cfg)
% Per (spec, combo, time t >= evalStartInd) expanding-window alpha,
% MSE, and HAC t-stat, plus the per-spec best-combo time series for
% bestComboRegPanel. The previous implementation built a T x specs x combos
% `yF2Extended` tensor (~16 GB) and a `yPocketTime` of the same size; this
% version derives `yF2_ii`, `cVec_ii`, `pw_ii`, `yPT_ii` per spec inside
% the loop, materializing only the per-spec slice (~1.6 GB transient).
%
% The benchmark spec (ii = numSpecifications) uses yF1pm1 — independent of
% combo — so we compute its alpha/mse/tStat once for jj=1 and replicate
% across the combo dim. Saves 9719/9720 of the inner work for that spec.

numSpecifications = size(yF2AllS, 2) + 1;   % 9 spec columns + 1 benchmark

alphMat   = NaN(T, numSpecifications,     numParamCombs);
mseMat    = NaN(T, numSpecifications,     numParamCombs);
tStatMat  = NaN(T, numSpecifications,     numParamCombs);
rmseAll   = NaN(T, numSpecifications - 1, numParamCombs);
bestAlphInd      = ones(T, numSpecifications);
bestTStatAlphInd = ones(T, numSpecifications);
bestRMSEInd      = ones(T, numSpecifications - 1);
yptR = NaN(T, numSpecifications - 1);
yptA = NaN(T, numSpecifications - 1);
yptT = NaN(T, numSpecifications - 1);

% Pre-compute the warmup time series from spec=end (yF1pm1), combo=1.
cBench = sqrt(var(yActual, 'omitmissing') ./ var(yF1pm1 .* yActual, 'omitmissing'));
pwBench = cBench * yF1pm1;
if weightRestriction
    pwBench(pwBench < 0) = 0;
    pwBench(pwBench > 2) = 2;
end
yPTBench = pwBench .* yActual;

for ii = 1:numSpecifications
    fprintf('  expandingAlpha spec %d/%d\n', ii, numSpecifications);
    if ii < numSpecifications
        yF2_ii = squeeze(yF2AllS(:, ii, :));        % T x numParamCombs
        isBenchSpec = false;
    else
        yF2_ii = yF1pm1;                             % T x 1 (combo-invariant)
        isBenchSpec = true;
    end

    alphTmp_ii  = NaN(T, numParamCombs);
    mseTmp_ii   = NaN(T, numParamCombs);
    tStatTmp_ii = NaN(T, numParamCombs);

    if isBenchSpec
        % Combo-invariant: compute once for jj=1, replicate.
        [alphCol, mseCol, tStatCol] = computeAlphaSeriesForCombo( ...
            yF2_ii, yActual, weightRestriction, evalStartInd, T);
        alphTmp_ii(:, :)  = repmat(alphCol,  1, numParamCombs);
        mseTmp_ii(:, :)   = repmat(mseCol,   1, numParamCombs);
        tStatTmp_ii(:, :) = repmat(tStatCol, 1, numParamCombs);
    else
        parfor (jj = 1:numParamCombs, useParfor(cfg))
            yPocketCur = yF2_ii(:, jj);
            [alphCol, mseCol, tStatCol] = computeAlphaSeriesForCombo( ...
                yPocketCur, yActual, weightRestriction, evalStartInd, T);
            alphTmp_ii(:, jj)  = alphCol;
            mseTmp_ii(:, jj)   = mseCol;
            tStatTmp_ii(:, jj) = tStatCol;
        end
    end

    alphMat(:, ii, :)  = reshape(alphTmp_ii,  [T, 1, numParamCombs]);
    mseMat(:, ii, :)   = reshape(mseTmp_ii,   [T, 1, numParamCombs]);
    tStatMat(:, ii, :) = reshape(tStatTmp_ii, [T, 1, numParamCombs]);

    % Per-spec best indices.
    [~, bestAlphInd(:, ii)]      = max(alphTmp_ii,  [], 2, 'omitmissing');
    [~, bestTStatAlphInd(:, ii)] = max(tStatTmp_ii, [], 2, 'omitmissing');

    if ii < numSpecifications
        rmse_ii = sqrt(cumsum(((yActual - yF2_ii).^2) ./ ((1:T)')));
        rmseAll(:, ii, :) = reshape(rmse_ii, [T, 1, numParamCombs]);
        [~, bestRMSEInd(:, ii)] = min(rmse_ii, [], 2, 'omitmissing');

        cVec_ii = sqrt(var(yActual, 'omitmissing') ./ var(yF2_ii .* yActual, 'omitmissing'));
        pw_ii = cVec_ii .* yF2_ii;
        if weightRestriction
            pw_ii(pw_ii < 0) = 0;
            pw_ii(pw_ii > 2) = 2;
        end
        yPT_ii = pw_ii .* yActual;

        yptR(1:evalStartInd, ii) = yPTBench(1:evalStartInd);
        yptA(1:evalStartInd, ii) = yPTBench(1:evalStartInd);
        yptT(1:evalStartInd, ii) = yPTBench(1:evalStartInd);
        bRind = bestRMSEInd(:, ii);
        bAind = bestAlphInd(:, ii);
        bTind = bestTStatAlphInd(:, ii);
        for t = evalStartInd+1:T
            yptR(t, ii) = yPT_ii(t, bRind(t-1));
            yptA(t, ii) = yPT_ii(t, bAind(t-1));
            yptT(t, ii) = yPT_ii(t, bTind(t-1));
        end
        clear pw_ii yPT_ii rmse_ii cVec_ii;
    end
    clear yF2_ii alphTmp_ii mseTmp_ii tStatTmp_ii;
    if ii == 1
        memSnapshot('inside-expandingAlpha (ii=1)');
    end
end
end


function [alphCol, mseCol, tStatCol] = computeAlphaSeriesForCombo( ...
    yPocketCur, yActual, weightRestriction, evalStartInd, T)
% Inner (t = evalStartInd:T) expanding-window OLS+HAC for one combo.
alphCol  = NaN(T, 1);
mseCol   = NaN(T, 1);
tStatCol = NaN(T, 1);
for t = evalStartInd:T
    cIn = sqrt(var(yActual(1:t), 'omitmissing') / ...
          var(yPocketCur(1:t) .* yActual(1:t), 'omitmissing'));
    pwIn = cIn * yPocketCur(1:t);
    if weightRestriction
        pwIn(pwIn < 0) = 0;
        pwIn(pwIn > 2) = 2;
    end
    yPT = pwIn .* yActual(1:t);
    s = regstats2Fast(yPT, yActual(1:t));
    alphCol(t)  = s.beta(1);
    mseCol(t)   = s.mse;
    tStatCol(t) = s.hac.t(1);
end
end



function [alphMatBest, tStatAlphMatBest, mseMatBest] = bestComboRegPanel( ...
    yptR, yptA, yptT, yActual, numSpecifications)
% Per spec, run a single regression of the best-combo portfolio time series
% on yActual. The yptR/yptA/yptT inputs are produced spec-by-spec inside
% expandingAlphaStreamed and have warmup rows already filled with the
% benchmark series.
alphMatBest      = NaN(3, numSpecifications);
tStatAlphMatBest = NaN(3, numSpecifications);
mseMatBest       = NaN(3, numSpecifications);

for ii = 1:(numSpecifications - 1)
    sR = regstats2Fast(yptR(:, ii), yActual);
    sA = regstats2Fast(yptA(:, ii), yActual);
    sT = regstats2Fast(yptT(:, ii), yActual);
    alphMatBest(:, ii)      = [sR.beta(1); sA.beta(1); sT.beta(1)];
    mseMatBest(:, ii)       = [sR.mse;     sA.mse;     sT.mse];
    tStatAlphMatBest(:, ii) = [sR.hac.t(1); sA.hac.t(1); sT.hac.t(1)];
end
end


function [dmMat, cwMat, cwDiffMat] = perComboForecastTests( ...
    yActual, yF1PMMat, yF2All, yComb1, yComb2, yComb3, ...
    pocketIndAll, combPocket, numParamCombs, cfg)
% Per-combo: build the {9,1} cell array of fDM/fCW/pocketInd inputs and
% run computeForecastTests. Each combo is independent -> parfor.
% Differentials are derived on the fly per combo from yF2All / yComb*
% rather than precomputed full-grid arrays (~30 GB savings vs. the
% previous implementation).
yF1pm1 = yF1PMMat(:, 1);
dmTmp     = NaN(27, numParamCombs);
cwTmp     = NaN(27, numParamCombs);
cwDiffTmp = NaN(9,  numParamCombs);

parfor (jj = 1:numParamCombs, useParfor(cfg))
    fDMCell  = cell(9, 1);
    fCWCell  = cell(9, 1);
    pIndCell = cell(9, 1);
    for v = 1:6
        yv = yF2All(:, v, jj);
        fDMCell{v}  = (yActual - yF1PMMat(:, v)).^2 - (yActual - yv).^2;
        fCWCell{v}  = fDMCell{v} + (yF1PMMat(:, v) - yv).^2;
        pIndCell{v} = pocketIndAll(:, v, jj);
    end
    yC1 = yComb1(:, jj);
    yC2 = yComb2(:, jj);
    yC3 = yComb3(:, jj);
    cp  = combPocket(:, jj);
    fDMCell{7} = (yActual - yF1pm1).^2 - (yActual - yC1).^2;
    fCWCell{7} = fDMCell{7} + (yF1pm1 - yC1).^2;
    fDMCell{8} = (yActual - yF1pm1).^2 - (yActual - yC2).^2;
    fCWCell{8} = fDMCell{8} + (yF1pm1 - yC2).^2;
    fDMCell{9} = (yActual - yF1pm1).^2 - (yActual - yC3).^2;
    fCWCell{9} = fDMCell{9} + (yF1pm1 - yC3).^2;
    pIndCell{7} = cp;
    pIndCell{8} = cp;
    pIndCell{9} = cp;

    [dm, cw, cwD] = computeForecastTests(fDMCell, fCWCell, pIndCell);
    dmTmp(:, jj)     = dm(:);
    cwTmp(:, jj)     = cw(:);
    cwDiffTmp(:, jj) = cwD;
end

dmMat     = reshape(dmTmp, [9, 3, numParamCombs]);
cwMat     = reshape(cwTmp, [9, 3, numParamCombs]);
cwDiffMat = cwDiffTmp;
end


function [cwBestMat, cwDiffBestMat] = bestComboCWStats( ...
    yActual, yF1PMMat, yF2All, yComb1, yComb2, yComb3, ...
    pocketIndAll, combPocket, ...
    bestRMSEInd, bestAlphInd, bestTStatAlphInd, T, numSpecifications)
% For each (spec, t > 1), pick the per-day best combo's CW differential
% and pocket indicator, then compute full / in-pocket / out-of-pocket
% HAC t-stats and the in-vs-out Welch t-statistic. CW differentials are
% derived on demand from yF2All / yComb1/2/3 rather than from a
% precomputed fCWVar/fComb*CW (~30 GB savings).
cwBestMat     = NaN(numSpecifications - 1, 3, 3);
cwDiffBestMat = NaN(numSpecifications - 1, 3);
yF1pm1 = yF1PMMat(:, 1);

for ii = 1:(numSpecifications - 1)
    if ii <= 6
        v   = ii;
        yv  = squeeze(yF2All(:, v, :));               % T x numParamCombs
        srcCW = (yActual - yF1PMMat(:, v)).^2 - ...
                ((yActual - yv).^2 - (yF1PMMat(:, v) - yv).^2);
        srcPocket = squeeze(pocketIndAll(:, v, :));    % T x numParamCombs logical
    else
        switch ii
            case 7; yC = yComb1;
            case 8; yC = yComb2;
            case 9; yC = yComb3;
        end
        srcCW = (yActual - yF1pm1).^2 - ...
                ((yActual - yC).^2 - (yF1pm1 - yC).^2);
        srcPocket = combPocket;                        % T x numParamCombs logical
    end

    fCWBest       = NaN(T, 3);
    pocketIndBest = false(T, 3);
    for t = 2:T
        fCWBest(t, 1) = srcCW(t, bestRMSEInd(t-1, ii));
        fCWBest(t, 2) = srcCW(t, bestAlphInd(t-1, ii));
        fCWBest(t, 3) = srcCW(t, bestTStatAlphInd(t-1, ii));
        pocketIndBest(t, 1) = srcPocket(t, bestRMSEInd(t-1, ii));
        pocketIndBest(t, 2) = srcPocket(t, bestAlphInd(t-1, ii));
        pocketIndBest(t, 3) = srcPocket(t, bestTStatAlphInd(t-1, ii));
    end
    for jj = 1:3
        rowCW = computeDMCWStats(fCWBest(:, jj), pocketIndBest(:, jj), 'cw');
        cwBestMat(ii, :, jj) = rowCW(1:3);
        cwDiffBestMat(ii, jj) = rowCW(4);
    end
end
end
