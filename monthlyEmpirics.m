function monthlyEmpirics(cfg)
% MONTHLYEMPIRICS  Monthly-frequency out-of-sample forecasting analysis
% with one-sided kernel and G-prior Bayesian shrinkage.
%
% Inputs (read from cfg.paths.data):
%   Monthly_Predictors.xlsx
%
% Outputs (written to cfg.paths.results/monthly/):
%   forecastResultsMonthly_*.mat   per spec; per-variable forecasts,
%                                  pocket indicators, durations, integral R^2.
%   OOSResultsMonthly_*.mat        per spec; combination forecasts and the
%                                  Diebold-Mariano / Clark-West / portfolio
%                                  panels.
%
% Tables/figures consumed by:
%   Tables A.5, A.6
%
% Runtime: ~1-3 minutes for the full 3-spec sweep.
%
% R1 hook: cfg.shrinkageTarget ('zero' / 'prevailingMean').

if nargin < 1 || isempty(cfg); cfg = default_config(); end
paths = cfg.paths;
resultsFolder = baselineResultsFolder(paths, cfg, 'monthly');

%% Spec grid
% Each row = [signRestriction, coefRestriction]. Window length is fixed
% at 2.5y and SED window at 1y for the published monthly results.
combs             = [0 0; 1 0; 1 1];
windowLengthYears = 2.5;
dmLengthYears     = 1;
varNames          = {'dp','tbl','tsp','rvar','mv','pc'};
nVars             = numel(varNames);
benchmarkFlag     = 'pm';
kernelMethod      = '1S';

%% Spec loop
specsToRun = 1:size(combs, 1);
if isfield(cfg, 'monthlySpecsSubset') && ~isempty(cfg.monthlySpecsSubset)
    specsToRun = cfg.monthlySpecsSubset;
end
for iSpec = specsToRun
    fprintf('Spec %d of %d\n', iSpec, size(combs, 1));
    signRestriction = combs(iSpec, 1);
    coefRestriction = combs(iSpec, 2);
    if coefRestriction
        restrictMat = [1, -1, 1, 1, 0, 0];
    else
        restrictMat = [0, 0, 0, 0, 0, 0];
    end
    [eLabel, dmLabel, rLabel] = buildSpecLabels(signRestriction, ...
        coefRestriction, windowLengthYears, dmLengthYears);
    fileName  = fullfile(resultsFolder, sprintf( ...
        'forecastResultsMonthly%s%s%s_%s_%s.mat', rLabel, eLabel, dmLabel, ...
        kernelMethod, benchmarkFlag));
    fileName2 = fullfile(resultsFolder, sprintf( ...
        'OOSResultsMonthly%s%s%s.mat', rLabel, eLabel, dmLabel));

    %% Data load and result preallocation
    data = readtable(fullfile(paths.data, 'Monthly_Predictors.xlsx'));
    data = data(2:end, :);            % match the original drop-first-row convention

    bufferInd = 2 * windowLengthYears * 12;
    TMax      = size(data, 1) - bufferInd;
    R         = preallocateResults(TMax, nVars, 100);

    pocketIndicesLPM = [];

    %% Per-variable estimation
    for varNum = 1:nVars
        varName = varNames{varNum};

        % Predictor selection. Monthly's 'pc' computes recursive PC
        % inline (vs daily which uses a pre-computed pcPredictor.mat).
        if strcmp(varName, 'mv')
            predictorMat = [data.dp, data.tbl, data.tsp, data.rvar];
            predictors   = predictorMat;
        elseif strcmp(varName, 'pc')
            predictorMat = (predictorMat - mean(predictorMat, 'omitnan')) ./ ...
                std(predictorMat, 'omitnan');
            predictors   = computeRecursivePC(predictorMat);
        else
            predictors = selectPredictor(varName, data, paths);
            predictorMat = [data.dp, data.tbl, data.tsp, data.rvar];  % for pc trim sizing
        end

        % Trim, stagger, and build the prevailing-mean naive forecast.
        S = prepareDataForEstimation(varName, predictors, data.exret, ...
            data.rf, data.Date, predictorMat);

        % Kernel + estimation hyperparameters.
        [h, out] = computeKernelBandwidth(windowLengthYears, S.T, 'monthly');
        [Ke1, ~] = oneSidedKernel();
        hSpec = struct('h', h, 'out', out, 'order', 0, 'Ke1', Ke1, ...
            'dmLengthYears', dmLengthYears, 'obsPerYear', 12, ...
            'minRT', 1, ...                                               % monthly: 1 month (cfg.minPocketDays = 21 is for daily)
            'weightWindow', 12, ...                                       % preserve original hard-coded value
            'sedTimeTrend', false);

        % Run kernel forecasting + G-prior + pocket detection.
        Rv = estimateVariableForecasts(S, hSpec, restrictMat(varNum), ...
            signRestriction, benchmarkFlag, cfg);

        % Snapshot dp-aligned series and compute the LPM pocket indicator
        % (used by the forecast-combination and figure rendering).
        if varNum == 1
            fDM2 = (Rv.yActual - Rv.yF1PM).^2 - (Rv.yActual - Rv.yF1).^2;
            [~, ~, ~, ~, ~, fDMHat2] = lps1_v2(fDM2, ones(size(fDM2)), ...
                dmLengthYears * 12 / numel(fDM2), 0, Ke1, [], 0);
            pocketIndicesLPM = detectPocketsFromSED(fDMHat2, ...
                numel(Rv.r2LocalPM), 1);   % monthly: 1 month minimum
            R.yActual    = Rv.yActual;
            R.riskFree   = Rv.riskFree;
        end

        R = packVariableResults(R, Rv, varNum, out);
    end

    R.dateVec           = data.Date(end-TMax+1:end);
    R.windowLengthYears = windowLengthYears;
    R.signRestriction   = signRestriction;
    R.coefRestriction   = coefRestriction;
    R.pocketIndicesLPM  = pocketIndicesLPM;

    % forecastResults persists the un-trimmed, un-scaled forecast
    % matrices (column-wise leading NaNs aligned to each variable's
    % predictor history).
    saveForecastResults(fileName, R);

    %% Forecast combination
    R = trimAndAlignResults(R, 100, ...
        {'yActual','yF1Mat','yF1PMMat','yF2Mat','riskFree'});
    R.pocketIndicesLPM = R.pocketIndicesLPM(end-size(R.yF1Mat,1)+1:end);

    combPocket = any(R.pocketIndMat(:, 1:4), 2);
    yComb2 = combineInPocketAverage(R);
    yComb3 = mean(R.yF2Mat(:, 1:4), 2, 'omitnan');

    % Univariate + LPM forecast differentials BEFORE the yF2Mat mutation.
    [fDMCell, fCWCell, pIndCell] = buildUnivariateAndLpmTestCells(R);

    % Mutate yF2Mat / yF1Mat: in-pocket TVC, out-of-pocket benchmark.
    R.yF2Mat = R.yF2Mat .* R.pocketIndMat   + R.yF1PMMat .* (~R.pocketIndMat);
    R.yF1Mat = R.yF1Mat .* R.pocketIndicesLPM + R.yF1PMMat .* (~R.pocketIndicesLPM);

    % yComb1 uses the MUTATED yF2Mat.
    yComb1 = mean(R.yF2Mat(:, 1:4), 2, 'omitnan');

    [fDMCell, fCWCell, pIndCell] = appendCombinationTestCells( ...
        fDMCell, fCWCell, pIndCell, R, yComb1, yComb2, yComb3, combPocket);

    %% Portfolio metrics + performance regressions
    yForecastMat = [R.yF2Mat(:,1:6), yComb1, yComb2, yComb3, ...
                    R.yF1Mat(:,1), R.yF1PMMat(:,1)];
    F = constructPortfolioFactors(yForecastMat, R.yActual, R.riskFree, 'monthly', true, cfg.adjCostBps);
    P = computePerformanceMetrics(F.yPocketTime, R.yActual, ...
        F.ySigFactor, F.yMomFactor, 'monthly');

    %% Forecast differential tests + final OOS save
    [dmMat, cwMat, cwDiffMat] = computeForecastTests(fDMCell, fCWCell, pIndCell);
    % Match the archive convention: only univariate (1-6) and comb3 (9)
    % carry an in-vs-out difference stat. comb1, comb2, and LPM are NaN
    % in tables A.5 / A.6 of the published reply.
    cwDiffMat([7, 8, 10]) = NaN;
    [cumsumPlotIn, cumsumPlotOut] = buildCumsumPlots(R, ...
        yComb1, yComb2, yComb3, combPocket);
    saveOOSResults(fileName2, dmMat, cwMat, cwDiffMat, P, F, R, ...
        cumsumPlotIn, cumsumPlotOut);
end
end


% ====================================================================
%  Local helpers — entry-point-specific.
% ====================================================================

function [eLabel, dmLabel, rLabel] = buildSpecLabels(sR, cR, windowYrs, dmYrs)
if windowYrs >= 1
    eLabel = sprintf('_%gyE', windowYrs);
else
    eLabel = sprintf('_%gmE', 12 * windowYrs);
end
if dmYrs >= 1
    dmLabel = sprintf('_%gyDM', dmYrs);
else
    dmLabel = sprintf('_%gmDM', dmYrs * 12);
end
if ~sR && ~cR
    rLabel = '_1';
elseif sR && ~cR
    rLabel = '_2';
else
    rLabel = '_3';
end
end


function predictors = computeRecursivePC(predictorMat)
% Monthly-only: recursively expanding-window first principal component.
% At each t, fit PCA on rows firstInd:t and store the first PC.
firstInd = max(sum(isnan(predictorMat))) + 1;
nx       = size(predictorMat, 2);
N        = size(predictorMat, 1);
predictors = NaN(N, N - firstInd - nx + 1);
for t = firstInd + nx : N
    [~, scores] = pca(predictorMat(firstInd:t, :));
    predictors(firstInd:t, t - firstInd - nx + 1) = scores(:, 1);
end
end


function R = preallocateResults(TMax, nVars, maxPockets)
R.yF1Mat        = NaN(TMax, nVars);
R.yF1PMMat      = NaN(TMax, nVars);
R.yF2Mat        = NaN(TMax, nVars);
R.yF2AltMat     = NaN(TMax, nVars);
R.aMat          = NaN(TMax, nVars);
R.gam1Mat       = NaN(TMax, nVars);
R.pocketIndMat  = NaN(TMax, nVars);
R.r2Mat         = NaN(TMax, nVars);
R.r21yMat       = NaN(TMax, nVars);
R.r2PMMat       = NaN(TMax, nVars);
R.r2GamMat      = NaN(TMax, nVars);
R.weightMatG    = NaN(TMax, 2, nVars);
R.durationMat   = NaN(maxPockets, nVars);
R.integralR2Mat = NaN(maxPockets, nVars);
R.fDMMat        = NaN(TMax, nVars);
end


function R = packVariableResults(R, Rv, varNum, out)
% Defensive right-aligned assignment: when a sub-sample's predictor trim
% combined with FP-quantization of the kernel exclusion `out` makes the
% per-variable forecast vector ONE row longer than TMax, drop the oldest
% (warmup-boundary) element so the right-aligned copy still fits. The
% full-sample case never overflows so behavior is unchanged there.
R.yF1Mat     = assignRightAligned(R.yF1Mat,     Rv.yF1,    varNum);
R.yF1PMMat   = assignRightAligned(R.yF1PMMat,   Rv.yF1PM,  varNum);
R.yF2Mat     = assignRightAligned(R.yF2Mat,     Rv.yF2,    varNum);
R.yF2AltMat  = assignRightAligned(R.yF2AltMat,  Rv.yF2Alt, varNum);
aTrim        = Rv.a(out+1:end);
R.aMat       = assignRightAligned(R.aMat,       aTrim,     varNum);
R.gam1Mat    = assignRightAligned(R.gam1Mat,    Rv.gamDM(:,1),    varNum);
R.pocketIndMat = assignRightAligned(R.pocketIndMat, Rv.pocketIndices, varNum);
R.r2Mat      = assignRightAligned(R.r2Mat,      Rv.r2Local,       varNum);
R.r21yMat    = assignRightAligned(R.r21yMat,    Rv.r2Local1y(end-numel(Rv.r2Local)+1:end), varNum);
R.r2PMMat    = assignRightAligned(R.r2PMMat,    Rv.r2LocalPM,     varNum);
R.r2GamMat   = assignRightAligned(R.r2GamMat,   Rv.r2Gam,         varNum);
R.durationMat(1:numel(Rv.pocketLengths), varNum) = Rv.pocketLengths;
R.integralR2Mat(1:numel(Rv.integralR2Pocket), varNum) = Rv.integralR2Pocket;
R.fDMMat     = assignRightAligned(R.fDMMat,     [NaN; Rv.fDMHat(1:end-1)], varNum);
% weightMatG is 3-D (TMax x 2 x nVars); apply the same truncation rule.
nW = size(Rv.weightG, 1);
nRowsW = size(R.weightMatG, 1);
if nW > nRowsW
    Rv.weightG = Rv.weightG(end-nRowsW+1:end, :);
    nW = nRowsW;
end
R.weightMatG(end-nW+1:end, :, varNum) = Rv.weightG;
end


function M = assignRightAligned(M, vec, col)
% Right-align `vec` into column `col` of matrix `M`. If vec is longer than
% size(M, 1), drop the oldest elements so the copy still fits.
nRows = size(M, 1);
n = numel(vec);
if n > nRows
    vec = vec(end-nRows+1:end);
    n = nRows;
end
M(end-n+1:end, col) = vec;
end


function yComb2 = combineInPocketAverage(R)
T = size(R.yF1Mat, 1);
yComb2 = NaN(T, 1);
for ii = 1:T
    if any(R.pocketIndMat(ii, 1:4))
        yComb2(ii) = mean(R.yF2Mat(ii, logical(R.pocketIndMat(ii, 1:4))));
    else
        yComb2(ii) = R.yF1PMMat(ii, 1);
    end
end
end


function [fDMCell, fCWCell, pIndCell] = buildUnivariateAndLpmTestCells(R)
% Univariate (1..6) + LPM (10) cells. Combination cells (7..9) are
% appended after yComb1 is computed. Total = 10 specs (no erL in monthly).
fDMCell  = cell(10, 1);
fCWCell  = cell(10, 1);
pIndCell = cell(10, 1);
for ii = 1:6
    [fCWCell{ii}, fDMCell{ii}] = forecastDifferentials( ...
        R.yActual, R.yF1PMMat(:,ii), R.yF2Mat(:,ii));
    pIndCell{ii} = logical(R.pocketIndMat(:,ii));
end
[fCWCell{10}, fDMCell{10}] = forecastDifferentials( ...
    R.yActual, R.yF1PMMat(:,1), R.yF1Mat(:,1));
pIndCell{10} = logical(R.pocketIndicesLPM);
end


function [fDMCell, fCWCell, pIndCell] = appendCombinationTestCells( ...
    fDMCell, fCWCell, pIndCell, R, yComb1, yComb2, yComb3, combPocket)
yF1pm1 = R.yF1PMMat(:, 1);
combos = [yComb1, yComb2, yComb3];
for cc = 1:3
    yC = combos(:, cc);
    fDMCell{6+cc} = (R.yActual - yF1pm1).^2 - (R.yActual - yC).^2;
    fCWCell{6+cc} = fDMCell{6+cc} + (yF1pm1 - yC).^2;
    pIndCell{6+cc} = combPocket;
end
end


function [cumsumPlotIn, cumsumPlotOut] = buildCumsumPlots(R, yComb1, yComb2, yComb3, combPocket)
yF2MatIn  = R.yF2Mat .* R.pocketIndMat   + R.yF1PMMat .* (~R.pocketIndMat);
yF2MatOut = R.yF2Mat .* (~R.pocketIndMat) + R.yF1PMMat .* R.pocketIndMat;
fDMyF2In  = (R.yActual - R.yF1PMMat).^2 - (R.yActual - yF2MatIn).^2;
fDMyF2Out = (R.yActual - R.yF1PMMat).^2 - (R.yActual - yF2MatOut).^2;

yCombIn  = [yComb1.*combPocket    + R.yF1PMMat(:,1).*(~combPocket), ...
            yComb2.*combPocket    + R.yF1PMMat(:,1).*(~combPocket), ...
            yComb3.*combPocket    + R.yF1PMMat(:,1).*(~combPocket)];
yCombOut = [yComb1.*(~combPocket) + R.yF1PMMat(:,1).*combPocket, ...
            yComb2.*(~combPocket) + R.yF1PMMat(:,1).*combPocket, ...
            yComb3.*(~combPocket) + R.yF1PMMat(:,1).*combPocket];
fDMCombIn  = (R.yActual - R.yF1PMMat(:,1:3)).^2 - (R.yActual - yCombIn).^2;
fDMCombOut = (R.yActual - R.yF1PMMat(:,1:3)).^2 - (R.yActual - yCombOut).^2;

cumsumPlotIn  = cumsum([fDMyF2In,  fDMCombIn]);
cumsumPlotOut = cumsum([fDMyF2Out, fDMCombOut]);
end


function saveForecastResults(fileName, R)
% The unpacked locals look unused to the linter (NASGU) but `save` picks
% them up by name; the suppressions silence those false positives.
yF1Mat       = R.yF1Mat;
yF2Mat       = R.yF2Mat;
aMat         = R.aMat;
pocketIndMat = R.pocketIndMat;
durationMat  = R.durationMat;
r2Mat        = R.r2Mat;
integralR2Mat = R.integralR2Mat;
dateVec      = R.dateVec;
fDMMat       = R.fDMMat;
gam1Mat      = R.gam1Mat;
yF1PMMat     = R.yF1PMMat;
r21yMat      = R.r21yMat;
r2PMMat      = R.r2PMMat;
r2GamMat     = R.r2GamMat;
yF2AltMat    = R.yF2AltMat;
weightMatG   = R.weightMatG;
yActual      = R.yActual;
riskFree     = R.riskFree;
windowLengthYears = R.windowLengthYears;
signRestriction   = R.signRestriction;
coefRestriction   = R.coefRestriction;
pocketIndicesLPM  = R.pocketIndicesLPM;
save(fileName, 'yF1Mat','yF2Mat','aMat','pocketIndMat','durationMat', ...
    'r2Mat','integralR2Mat','dateVec','fDMMat','gam1Mat','yF1PMMat', ...
    'r21yMat','r2PMMat','r2GamMat','yF2AltMat','weightMatG', ...
    'yActual','riskFree','windowLengthYears','signRestriction', ...
    'coefRestriction','pocketIndicesLPM');
end


function saveOOSResults(fileName, dmMat, cwMat, cwDiffMat, P, F, R, ...
    cumsumPlotIn, cumsumPlotOut)
% Pack and save the variable layout consumed by tableA6.
% Several inputs are consumed by `save` by name below.
statMat       = cwMat;
statDiffMat   = cwDiffMat;
econMat       = [P.alphaAnnual', P.tStatMat(1,:)', P.SRAnnual'];
coefMatVol    = P.coefMatVol;    tStatMatVol = P.tStatMatVol;
coefMatMom    = P.coefMatMom;    tStatMatMom = P.tStatMatMom;
coefMat3Fac   = P.coefMat3Fac;   tStatMat3Fac = P.tStatMat3Fac;
portfolioWeight = F.portfolioWeight;
dateVec = R.dateVec;
save(fileName, 'statMat','econMat','coefMatVol','tStatMatVol', ...
    'coefMatMom','tStatMatMom','coefMat3Fac','tStatMat3Fac', ...
    'cumsumPlotIn','cumsumPlotOut','dateVec', ...
    'portfolioWeight','statDiffMat','dmMat');
end
