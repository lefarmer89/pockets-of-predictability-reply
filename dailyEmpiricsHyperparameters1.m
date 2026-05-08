function dailyEmpiricsHyperparameters1(cfg)
%DAILYEMPIRICSHYPERPARAMETERS1  Generate OOS forecasts for every
% combination in the 9,720-row hyperparameter grid (winsorization,
% g-prior, weight bounds, window, time-trend, min-pocket; both shrinkage
% targets). Run dailyEmpiricsHyperparameters2.m next to consume the
% output.
%
% Inputs (from cfg.paths.data):
%   Daily_Predictors.xlsx, pcPredictor.mat
%
% Outputs (to cfg.paths.results/hyperparameters/):
%   forecastResults_*_HyperR1.mat   (one per signSpec)
%
% Consumers: Figures 1, A.1, A.2, A.3 and the Adaptive panels of Tables 1
% and A.3 (via dailyEmpiricsHyperparameters2 and *Marginals).
%
% Runtime: ~4-8 hours on an 8-core machine for the full 9,720-combo
% grid. Smoke (10 combos, 1 signSpec) runs in ~35 s.

if nargin < 1 || isempty(cfg); cfg = default_config(); end
paths = cfg.paths;
resultsFolder = robustnessResultsFolder(paths, cfg);

%% Spec grid (signRestriction × coefRestriction; window/SED hard-coded)
combs             = [0,0,2.7,1; 1,0,2.7,1; 1,1,2.7,1];
varNames          = {'dp','tbl','tsp','rvar','mv','pc'};
nVars             = numel(varNames);
benchmarkFlag     = 'pm';
kernelMethod      = '1S';
[paramCombs, idx] = constructHyperparameterGrid();
if isfield(cfg, 'hyperparameterComboSubset') && ~isempty(cfg.hyperparameterComboSubset)
    paramCombs = paramCombs(cfg.hyperparameterComboSubset, :);
end
numParamCombs     = size(paramCombs, 1);

%% Data load (shared across signSpecs)
data = readtable(fullfile(paths.data, 'Daily_Predictors.xlsx'));
loadedPC       = load(fullfile(paths.data, 'pcPredictor.mat'));
pcPredictorAll = loadedPC.predictors;
clear loadedPC;

%% Spec loop
% Default: iSpec in {1, 2} (no sign restriction; sign-restricted). Set
% cfg.hyperparameterSpecs to a subset (e.g., 2) to restrict the sweep.
if isfield(cfg, 'hyperparameterSpecs') && ~isempty(cfg.hyperparameterSpecs)
    iSpecRange = cfg.hyperparameterSpecs;
else
    iSpecRange = [1, 2];
end
for iSpec = iSpecRange
    fprintf('Spec %d\n', iSpec);
    spec = combs(iSpec, :);
    signRestriction   = spec(1);
    coefRestriction   = spec(2);
    windowLengthYears = spec(3);
    dmLengthYears     = spec(4);
    if coefRestriction
        restrictMat = [1, -1, 1, 1, 0, 0];
    else
        restrictMat = [0, 0, 0, 0, 0, 0];
    end
    [eLabel, dmLabel, rLabel] = buildSpecLabels(spec);
    fileName = fullfile(resultsFolder, sprintf( ...
        'forecastResults%s%s%s_%s_%s_HyperR1.mat', rLabel, eLabel, ...
        dmLabel, kernelMethod, benchmarkFlag));

    %% Per-spec result preallocation
    [hPre, ~] = computeKernelBandwidth(windowLengthYears, ...
        numel(data.dp) - 1, 'daily');
    bufferInd = 2 * floor(hPre * (numel(data.dp) - 1));
    TMax      = size(data, 1) - bufferInd;

    yF1PMMat     = NaN(TMax, nVars);
    yF2All       = NaN(TMax, nVars, numParamCombs);
    pocketIndAll = NaN(TMax, nVars, numParamCombs);
    fDMHatAll    = NaN(TMax, nVars, numParamCombs);

    predictorMat = [data.dp, data.tbl, data.tsp, data.rvar];

    %% Per-variable estimation
    for varNum = 1:nVars
        varName = varNames{varNum};

        % Predictor selection. 'pc' loads the pre-computed pcPredictor.mat
        % from cfg.paths.data; 'mv' uses the running normalized stack.
        if strcmp(varName, 'pc')
            predictorMat = (predictorMat - mean(predictorMat, 'omitnan')) ./ ...
                std(predictorMat, 'omitnan');
            predictors = pcPredictorAll;   % already filtered by sampleMask
        elseif strcmp(varName, 'mv')
            predictors = predictorMat;
        else
            predictors = data.(varName);
        end

        S = prepareDataForEstimation(varName, predictors, data.exret, ...
            data.rf, data.Date, predictorMat);

        [h, out] = computeKernelBandwidth(windowLengthYears, S.T, 'daily');
        [Ke1, ~] = oneSidedKernel();

        % --- Pre-compute the building blocks the inner combo loop reuses.
        [~, ~, r2LocalPM, yForecastPM, ~, ~] = lps1_v2(S.y, ones(S.T,1), ...
            h, 0, Ke1, [], 0); %#ok<ASGLU>
        if signRestriction
            yForecastPM(yForecastPM < 0) = 0;
        end
        [~, ~, r2Local, yForecastOG, ~, ~] = lps1_v2(S.y, ...
            [ones(S.T,1), S.X], h, 0, Ke1, [], restrictMat(varNum));
        r2Local = [0; r2Local(1:end-1)];

        % --- Hierarchical caching + parfor.
        % Forecasts depend on paramCombs cols 1-6: 972 unique forecasts
        % among 9,720 combos. SEDs depend on cols 1-7: 1,944 unique.
        % Pocket detection depends on col 8 and is per-combo. Cell-array
        % caches accommodate variable-length output from evaluateCombo
        % (lps1_v2 trims internally by its own warmup of ~252 obs).

        % --- Group identification.
        [~, ~, grp6] = unique(paramCombs(:, 1:6), 'rows', 'stable');
        nGrp6 = max(grp6);
        [~, ~, grp7] = unique(paramCombs(:, 1:7), 'rows', 'stable');
        nGrp7 = max(grp7);
        grp6Rep = zeros(nGrp6, 1);
        for g = 1:nGrp6; grp6Rep(g) = find(grp6 == g, 1); end
        grp7Rep = zeros(nGrp7, 1);
        for g = 1:nGrp7; grp7Rep(g) = find(grp7 == g, 1); end

        % --- Compute jj==1 result up-front for metadata save.
        fprintf('  varNum=%d: computing Rfirst...\n', varNum);
        tStart = tic;
        Rfirst = evaluateCombo(1, paramCombs, idx, S, out, ...
            yForecastOG, yForecastPM, r2Local, signRestriction, ...
            dmLengthYears, Ke1, [], cfg);
        fprintf('  varNum=%d: Rfirst done in %.1fs\n', varNum, toc(tStart));
        if varNum == 1
            yActual  = Rfirst.yActual;
            riskFree = Rfirst.riskFree;
            dateVec  = data.Date(end-TMax+1:end);
            save(fileName, 'yActual','riskFree','dateVec','paramCombs', ...
                'numParamCombs','idx','-v7.3');
            fprintf('  varNum=%d: metadata save done\n', varNum);
        end

        % --- Phase A: forecast cache (parfor over unique 6-col groups).
        fprintf('  varNum=%d: Phase A (%d 6-col groups)...\n', varNum, nGrp6);
        tStart = tic;
        forecastCache = cell(1, nGrp6);
        parfor (g = 1:nGrp6, useParfor(cfg))
            R = evaluateCombo(grp6Rep(g), paramCombs, idx, S, out, ...
                yForecastOG, yForecastPM, r2Local, signRestriction, ...
                dmLengthYears, Ke1, [], cfg); %#ok<PFBNS>
            forecastCache{g} = R.yF2;
        end
        fprintf('  varNum=%d: Phase A done in %.1fs\n', varNum, toc(tStart));

        % --- Phase B: SED cache (parfor over unique 7-col groups).
        fprintf('  varNum=%d: Phase B (%d 7-col groups)...\n', varNum, nGrp7);
        tStart = tic;
        sedCache = cell(1, nGrp7);
        parfor (g = 1:nGrp7, useParfor(cfg))
            R = evaluateCombo(grp7Rep(g), paramCombs, idx, S, out, ...
                yForecastOG, yForecastPM, r2Local, signRestriction, ...
                dmLengthYears, Ke1, [], cfg); %#ok<PFBNS>
            sedCache{g} = R.fDMHat;
        end
        fprintf('  varNum=%d: Phase B done in %.1fs\n', varNum, toc(tStart));

        % --- Phase C: per-combo pocket detection using cached SED.
        fprintf('  varNum=%d: Phase C (%d combos)...\n', varNum, numParamCombs);
        tStart = tic;
        nF2 = size(yF2All, 1);
        nPI = size(pocketIndAll, 1);
        nFD = size(fDMHatAll, 1);
        yF2All_var       = NaN(nF2, numParamCombs);
        pocketIndAll_var = NaN(nPI, numParamCombs);
        fDMHatAll_var    = NaN(nFD, numParamCombs);

        parfor (jj = 1:numParamCombs, useParfor(cfg))
            yF2    = forecastCache{grp6(jj)}; %#ok<PFBNS>
            fDMHat = sedCache{grp7(jj)};       %#ok<PFBNS>

            minRT = paramCombs(jj, idx.minPocket);
            pocketIndices = pocketDetectorWithShift(fDMHat, r2Local, minRT);

            yF2col = NaN(nF2, 1);
            yF2col(end-numel(yF2)+1:end) = yF2;
            yF2All_var(:, jj) = yF2col;

            pIcol = NaN(nPI, 1);
            pIcol(end-numel(pocketIndices)+1:end) = pocketIndices;
            pocketIndAll_var(:, jj) = pIcol;

            fDcol = NaN(nFD, 1);
            fDcol(end-numel(fDMHat)+1:end) = fDMHat;
            fDMHatAll_var(:, jj) = fDcol;
        end
        fprintf('  varNum=%d: Phase C done in %.1fs\n', varNum, toc(tStart));

        fprintf('  varNum=%d: reshape...\n', varNum);
        tStart = tic;
        yF2All(:, varNum, :)       = reshape(yF2All_var,       [size(yF2All,1),       1, numParamCombs]);
        pocketIndAll(:, varNum, :) = reshape(pocketIndAll_var, [size(pocketIndAll,1), 1, numParamCombs]);
        fDMHatAll(:, varNum, :)    = reshape(fDMHatAll_var,    [size(fDMHatAll,1),    1, numParamCombs]);

        yF1PMMat(end-size(Rfirst.yF1PM,1)+1:end, varNum) = Rfirst.yF1PM;
        fprintf('  varNum=%d: reshape done in %.1fs; saving...\n', varNum, toc(tStart));
        tStart = tic;
        save(fileName, 'yF1PMMat','yF2All','pocketIndAll','fDMHatAll','-append');
        fprintf('  varNum=%d: save done in %.1fs\n', varNum, toc(tStart));
    end
end
end


% ====================================================================
%  Local helpers — entry-point-specific.
% ====================================================================

function [eLabel, dmLabel, rLabel] = buildSpecLabels(spec)
windowLengthYears = spec(3);
dmLengthYears     = spec(4);
if windowLengthYears >= 1
    eLabel = sprintf('_%gyE', windowLengthYears);
else
    eLabel = sprintf('_%gmE', 12 * windowLengthYears);
end
if windowLengthYears == 2.7
    eLabel = '_2.5yE';
end
if dmLengthYears >= 1
    dmLabel = sprintf('_%gyDM', dmLengthYears);
else
    dmLabel = sprintf('_%gmDM', dmLengthYears * 12);
end
if ~spec(1) && ~spec(2)
    rLabel = '_1';
elseif spec(1) && ~spec(2)
    rLabel = '_2';
else
    rLabel = '_3';
end
end


function R = evaluateCombo(jj, paramCombs, idx, S, out, ...
    yForecastOG, yForecastPM, r2Local, signRestriction, ...
    dmLengthYears, Ke1, prevCacheSlice, cfg)
% Per-combo body. prevCacheSlice carries the cached fDMHat for the
% combo's parent group; an empty value means recompute.

% --- Combo-specific winsorization (qLB drives qUB = 100 - qLB).
qLB = paramCombs(jj, idx.winsorization);
yForecast = winsorizeRecursiveFast(yForecastOG, qLB, 100 - qLB, out + 1);

% --- Combo-specific G-prior shrinkage.
g            = paramCombs(jj, idx.gPrior);
weightWindow = paramCombs(jj, idx.windowLength) * 252;
minW         = paramCombs(jj, idx.minWeight);
maxW         = paramCombs(jj, idx.maxWeight);
if paramCombs(jj, idx.benchmark)
    benchmarkForecast = yForecastPM;
else
    benchmarkForecast = computeBenchmarkForecast(yForecast, S.naiveForecast, cfg);
end

[yCombinedG, ~, ~] = gPriorWeightsFast( ...
    yForecast, S.y, S.naiveForecast, benchmarkForecast, weightWindow, out, ...
    g, minW, maxW);

% --- Align to OOS sample boundary; sign-restrict.
yF1PM   = S.naiveForecast(2*out:end);
yF1     = yForecastPM(out + 1:end);
yF2     = yCombinedG(out + 1:end);
if signRestriction
    yF2(yF2 < 0) = 0;
end
yActual  = S.y(2*out:end);
riskFree = S.riskFree(2*out:end);

% Forecast differential against the configured benchmark.
fDM = (yActual - yF1PM).^2 - (yActual - yF2).^2;

% --- SED smoothing for pocket detection (cached when first 7 paramCombs
%     columns match a prior combo).
if isempty(prevCacheSlice)
    if paramCombs(jj, idx.timeTrend)
        sedDesign = [ones(size(fDM)), (1:numel(fDM))' ./ numel(fDM)];
    else
        sedDesign = ones(size(fDM));
    end
    [~, ~, ~, ~, ~, fDMHat] = lps1_v2(fDM, sedDesign, ...
        dmLengthYears * 252 / numel(fDM), 0, Ke1, [], 0);
else
    fDMHat = prevCacheSlice;
end

% --- Pocket detection with combo-specific minRT shift.
minRT = paramCombs(jj, idx.minPocket);
pocketIndices = pocketDetectorWithShift(fDMHat, r2Local, minRT);

R.yF2           = yF2;
R.yF1PM         = yF1PM;
R.yF1           = yF1;
R.pocketIndices = pocketIndices;
R.fDMHat        = fDMHat;
R.yActual       = yActual;
R.riskFree      = riskFree;
end


function pocketIndices = pocketDetectorWithShift(fDMHat, r2Local, minRT)
% Pocket detector for the hyperparameter sweep. Differs from
% subroutines/diagnostics/detectPocketsFromSED in two ways: (a) drops
% pockets with length > minRT (NOT >=), and (b) shifts pocket starts
% by +minRT. detectPocketsFromSED matches this when minRT > 0; the
% Hyperparameters1 sweep also wants minRT == 0 to be a no-op.
preStart = numel(r2Local) - numel(fDMHat) + 1;
hits     = find([zeros(preStart, 1); fDMHat(1:end-1) > 0]);
pocketIndices = false(numel(r2Local), 1);
if isempty(hits); return; end

breaks    = [true; diff(hits) ~= 1];
runStarts = hits(breaks);
runEnds   = hits([find(breaks(2:end)); numel(hits)]);
periods   = [runStarts, runEnds];

keep    = (periods(:,2) - periods(:,1) + 1) > minRT;
periods = periods(keep, :);
periods(:,1) = periods(:,1) + minRT;

for k = 1:size(periods, 1)
    pocketIndices(periods(k,1) : periods(k,2)) = true;
end
end
