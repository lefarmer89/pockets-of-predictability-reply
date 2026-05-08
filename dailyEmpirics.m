function dailyEmpirics(cfg)
% DAILYEMPIRICS  Daily out-of-sample forecasting analysis with one-sided
% kernel and G-prior Bayesian shrinkage. Drives the Fixed-shrinkage panels
% of every main daily table.
%
% Inputs (read from cfg.paths.data):
%   Daily_Predictors.xlsx, pcPredictor.mat
%
% Outputs (written to cfg.paths.results/daily/, or to
% cfg.paths.results/<cfg.resultsSubdir>/ when set by
% dailyEmpiricsGPriorRobustness):
%   forecastResults_*.mat   one per spec; per-variable forecasts, pocket
%                           indicators, durations, integral R^2, plus the
%                           yF2MatS / pocketIndMatS / dateVecS appended
%                           panels consumed by figure2.
%   OOSResults_*.mat        one per spec; combination forecasts and the
%                           Diebold-Mariano / Clark-West / portfolio panels.
%
% Tables/figures consumed by:
%   Tables 1, 2, 3, A.1, A.2, A.3 (Fixed panels), figure 2.
%
% Runtime: ~5-15 minutes for the full 15-spec sweep (smoke run at
% one spec is ~15 s).
%
% R1 hook: cfg.shrinkageTarget ('zero' / 'prevailingMean').

if nargin < 1 || isempty(cfg); cfg = default_config(); end
paths = cfg.paths;
resultsFolder = baselineResultsFolder(paths, cfg, 'daily');

%% Spec grid
% Each row = [signRestriction, coefRestriction, windowLengthYears, dmLengthYears].
% Index 7 (sign-restricted, 2.5y window, 0.5y SED) is the published baseline.
combs = [
    0 0 2.7 1;   0 0 2.5 0.5; 0 0 2.5 1.5; 0 0 2 1; 0 0 3 1; ...
    1 0 2.7 1;   1 0 2.5 0.5; 1 0 2.5 1.5; 1 0 2 1; 1 0 3 1; ...
    1 1 2.7 1;   1 1 2.5 0.5; 1 1 2.5 1.5; 1 1 2 1; 1 1 3 1];
varNames     = {'dp','tbl','tsp','rvar','mv','pc','erL'};
nVars        = numel(varNames);
benchmarkFlag = 'pm';
kernelMethod  = '1S';

%% Spec loop
specsToRun = 1:size(combs, 1);
if isfield(cfg, 'dailySpecsSubset') && ~isempty(cfg.dailySpecsSubset)
    specsToRun = cfg.dailySpecsSubset;
end
for iSpec = specsToRun
    fprintf('Spec %d of %d\n', iSpec, size(combs, 1));
    spec = combs(iSpec, :);
    signRestriction   = spec(1);
    coefRestriction   = spec(2);
    windowLengthYears = spec(3);
    dmLengthYears     = spec(4);
    if coefRestriction
        restrictMat = [1, -1, 1, 1, 0, 0, 0];
    else
        restrictMat = [0, 0, 0, 0, 0, 0, 0];
    end
    [eLabel, dmLabel, rLabel] = buildSpecLabels(spec);
    fileName  = fullfile(resultsFolder, sprintf( ...
        'forecastResults%s%s%s_%s_%s.mat', rLabel, eLabel, dmLabel, ...
        kernelMethod, benchmarkFlag));
    fileName2 = fullfile(resultsFolder, sprintf( ...
        'OOSResults%s%s%s.mat', rLabel, eLabel, dmLabel));

    %% Data load and result preallocation
    data = readtable(fullfile(paths.data, 'Daily_Predictors.xlsx'));
    loadedPC = load(fullfile(paths.data, 'pcPredictor.mat'));
    pcPredictorAll = loadedPC.predictors;
    [hPre, ~] = computeKernelBandwidth(windowLengthYears, numel(data.dp) - 1, 'daily');
    bufferInd = 2 * floor(hPre * (numel(data.dp) - 1));
    TMax      = size(data, 1) - bufferInd;
    R         = preallocateResults(TMax, nVars, 500);
    R.econFactors = [data.recession(end-TMax+1:end), data.bwindex(end-TMax+1:end), ...
                     data.lfactor(end-TMax+1:end),  data.cgcoef(end-TMax+1:end)];

    predictorMat     = [data.dp, data.tbl, data.tsp, data.rvar];
    pocketIndicesLPM = [];

    %% Per-variable estimation
    for varNum = 1:nVars
        varName = varNames{varNum};

        % Predictor selection (pc/mv branch is inline because both
        % depend on the running normalization of predictorMat).
        if strcmp(varName, 'pc')
            predictorMat = (predictorMat - mean(predictorMat, 'omitnan')) ./ ...
                std(predictorMat, 'omitnan');
            predictors = pcPredictorAll;   % already filtered by sampleMask
        elseif strcmp(varName, 'mv')
            predictors = predictorMat;
        else
            predictors = selectPredictor(varName, data, paths);
        end

        % Trim, stagger, and build the prevailing-mean naive forecast.
        S = prepareDataForEstimation(varName, predictors, data.exret, ...
            data.rf, data.Date, predictorMat);

        % Kernel + estimation hyperparameters.
        [h, out] = computeKernelBandwidth(windowLengthYears, S.T, 'daily');
        [Ke1, ~] = oneSidedKernel();
        hSpec = struct('h', h, 'out', out, 'order', 0, 'Ke1', Ke1, ...
            'dmLengthYears', dmLengthYears, 'obsPerYear', 252, ...
            'minRT', cfg.minPocketDays, ...
            'weightWindow', round(cfg.weightWindowYears * 252), ...
            'sedTimeTrend', false);

        % Run kernel forecasting + G-prior shrinkage + pocket detection.
        Rv = estimateVariableForecasts(S, hSpec, restrictMat(varNum), ...
            signRestriction, benchmarkFlag, cfg);

        % LPM pocket detection (only for dp; consumed by combo + figure).
        % Snapshot dp-aligned series here; mirrors the original "save when
        % varName == 'dp'" / reload pattern.
        if varNum == 1
            fDM2 = (Rv.yActual - Rv.yF1PM).^2 - (Rv.yActual - Rv.yF1).^2;
            [~, ~, ~, ~, ~, fDMHat2] = lps1_v2(fDM2, ones(size(fDM2)), ...
                dmLengthYears * 252 / numel(fDM2), 0, Ke1, [], 0);
            pocketIndicesLPM = detectPocketsFromSED(fDMHat2, ...
                numel(Rv.r2LocalPM), cfg.minPocketDays);
            yAR1Full     = computeAR1Forecast(data.exret, 252);
            R.yActual    = Rv.yActual;
            R.riskFree   = Rv.riskFree;
            R.dependent  = data.exret(S.trim+1:end);
            R.yAR1       = yAR1Full(end-numel(Rv.yActual)+1:end);
        end

        R = packVariableResults(R, Rv, varNum, out);
    end

    R.dateVec           = data.Date(end-TMax+1:end);
    R.windowLengthYears = windowLengthYears;
    R.signRestriction   = signRestriction;
    R.coefRestriction   = coefRestriction;
    R.pocketIndicesLPM  = pocketIndicesLPM;

    % forecastResults persists the un-trimmed, un-scaled forecast matrices.
    saveForecastResults(fileName, R);

    %% Forecast combination (build the three combination panels)
    R = trimAndAlignResults(R, 100);
    R.pocketIndicesLPM = R.pocketIndicesLPM(end-size(R.yF1Mat,1)+1:end);

    % combPocket / yComb2 / yComb3 use the UNMUTATED yF2Mat (TVC-only).
    combPocket = any(R.pocketIndMat(:, 1:4), 2);
    yComb2 = combineInPocketAverage(R);
    yComb3 = mean(R.yF2Mat(:, 1:4), 2, 'omitnan');

    % Univariate + LPM forecast differentials are also computed BEFORE
    % the yF2Mat mutation so the test stats use raw TVC forecasts.
    [fDMCell, fCWCell, pIndCell] = buildUnivariateAndLpmTestCells(R);

    % Mutate yF2Mat / yF1Mat: in-pocket TVC, out-of-pocket benchmark.
    R.yF2Mat = R.yF2Mat .* R.pocketIndMat   + R.yF1PMMat .* (~R.pocketIndMat);
    R.yF1Mat = R.yF1Mat .* R.pocketIndicesLPM + R.yF1PMMat .* (~R.pocketIndicesLPM);

    % yComb1 uses the MUTATED yF2Mat.
    yComb1 = mean(R.yF2Mat(:, 1:4), 2, 'omitnan');

    % Combination test cells (must come AFTER yComb1).
    [fDMCell, fCWCell, pIndCell] = appendCombinationTestCells( ...
        fDMCell, fCWCell, pIndCell, R, yComb1, yComb2, yComb3, combPocket);

    % Append the *S (trimmed/scaled) panels used by figure2. The -append
    % to a freshly written file occasionally races Dropbox's incremental
    % sync and lands a "file appears to be corrupt" error; retry once
    % with a short pause to let Dropbox release the file handle.
    yF2MatS       = [R.yF2Mat, yComb1, yComb2, yComb3];
    pocketIndMatS = [R.pocketIndMat, combPocket, combPocket, combPocket];
    dateVecS      = R.dateVec;
    try
        save(fileName, 'yF2MatS', 'pocketIndMatS', 'dateVecS', '-append');
    catch ME
        if contains(ME.message, 'corrupt')
            pause(2);
            save(fileName, 'yF2MatS', 'pocketIndMatS', 'dateVecS', '-append');
        else
            rethrow(ME);
        end
    end

    %% Portfolio metrics + performance regressions
    yForecastMat = [R.yF2Mat(:,1:7), yComb1, yComb2, yComb3, ...
                    R.yF1Mat(:,1), R.yF1PMMat(:,1), R.yAR1];
    F = constructPortfolioFactors(yForecastMat, R.yActual, R.riskFree, 'daily', true, cfg.adjCostBps);
    P = computePerformanceMetrics(F.yPocketTime, R.yActual, ...
        F.ySigFactor, F.yMomFactor, 'daily');

    %% Forecast differential tests
    [dmMat, cwMat, cwDiffMat] = computeForecastTests(fDMCell, fCWCell, pIndCell);
    % Only univariate (rows 1-7) and comb3 (row 10) carry an in-vs-out
    % difference stat. comb1, comb2, and LPM are set NaN.
    % in table 3 (and Robustness panels) of the published reply.
    cwDiffMat([8, 9, 11]) = NaN;

    %% Cumsum diagnostics + final OOS save
    [cumsumPlotIn, cumsumPlotOut] = buildCumsumPlots(R, ...
        yComb1, yComb2, yComb3, combPocket);
    pocketIndExtend = [R.pocketIndMat, combPocket, combPocket, combPocket];

    saveOOSResults(fileName2, dmMat, cwMat, cwDiffMat, P, F, R, ...
        cumsumPlotIn, cumsumPlotOut, pocketIndExtend);
end
end


% ====================================================================
%  Local helpers — entry-point-specific. Reusable helpers live in
%  Replication Package/subroutines/.
% ====================================================================

function [eLabel, dmLabel, rLabel] = buildSpecLabels(spec)
% Filename labels for the spec row. Original convention preserved:
% the 2.7y inflated-bandwidth spec is labeled '_2.5yE' in output names.
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
% Slot per-variable outputs into the spec-level result matrices.
nF1   = numel(Rv.yF1);    R.yF1Mat(end-nF1+1:end, varNum)         = Rv.yF1;
nF1PM = numel(Rv.yF1PM);  R.yF1PMMat(end-nF1PM+1:end, varNum)     = Rv.yF1PM;
nF2   = numel(Rv.yF2);    R.yF2Mat(end-nF2+1:end, varNum)         = Rv.yF2;
nA    = numel(Rv.yF2Alt); R.yF2AltMat(end-nA+1:end, varNum)       = Rv.yF2Alt;
aTrim = Rv.a(out+1:end);
R.aMat(end-numel(aTrim)+1:end, varNum)        = aTrim;
R.gam1Mat(end-size(Rv.gamDM,1)+1:end, varNum) = Rv.gamDM(:,1);
R.pocketIndMat(end-numel(Rv.pocketIndices)+1:end, varNum) = Rv.pocketIndices;
nR2 = numel(Rv.r2Local);
R.r2Mat(end-nR2+1:end, varNum)    = Rv.r2Local;
R.r21yMat(end-nR2+1:end, varNum)  = Rv.r2Local1y(end-nR2+1:end);
R.r2PMMat(end-numel(Rv.r2LocalPM)+1:end, varNum) = Rv.r2LocalPM;
R.r2GamMat(end-numel(Rv.r2Gam)+1:end, varNum)    = Rv.r2Gam;
R.durationMat(1:numel(Rv.pocketLengths), varNum) = Rv.pocketLengths;
R.integralR2Mat(1:numel(Rv.integralR2Pocket), varNum) = Rv.integralR2Pocket;
nFD = numel(Rv.fDMHat);
R.fDMMat(end-nFD+1:end, varNum) = [NaN; Rv.fDMHat(1:end-1)];
R.weightMatG(end-size(Rv.weightG,1)+1:end, :, varNum) = Rv.weightG;
end


function yComb2 = combineInPocketAverage(R)
% yComb2: equal-weighted average of the in-pocket TVC forecasts at each
% time t. Falls back to the prevailing-mean benchmark when no variable is
% in-pocket at t.
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
% Univariate (1..7) + LPM (11) cells. Combination cells (8..10) are
% appended in a separate call after yComb1 is computed.
fDMCell  = cell(11, 1);
fCWCell  = cell(11, 1);
pIndCell = cell(11, 1);
for ii = 1:7
    [fCWCell{ii}, fDMCell{ii}] = forecastDifferentials( ...
        R.yActual, R.yF1PMMat(:,ii), R.yF2Mat(:,ii));
    pIndCell{ii} = logical(R.pocketIndMat(:,ii));
end
[fCWCell{11}, fDMCell{11}] = forecastDifferentials( ...
    R.yActual, R.yF1PMMat(:,1), R.yF1Mat(:,1));
pIndCell{11} = logical(R.pocketIndicesLPM);
end


function [fDMCell, fCWCell, pIndCell] = appendCombinationTestCells( ...
    fDMCell, fCWCell, pIndCell, R, yComb1, yComb2, yComb3, combPocket)
% Combination forecast tests (rows 8..10 of fDM/fCWCell). All three
% combos use combPocket as the in-pocket indicator.
yF1pm1 = R.yF1PMMat(:, 1);
combos = [yComb1, yComb2, yComb3];
for cc = 1:3
    yC = combos(:, cc);
    fDMCell{7+cc} = (R.yActual - yF1pm1).^2 - (R.yActual - yC).^2;
    fCWCell{7+cc} = fDMCell{7+cc} + (yF1pm1 - yC).^2;
    pIndCell{7+cc} = combPocket;
end
end


function [cumsumPlotIn, cumsumPlotOut] = buildCumsumPlots(R, yComb1, yComb2, yComb3, combPocket)
% Cumulative-sum-of-DM panels for figure rendering. Match the original
% benchmarkFlag2 == 1 branch (prevailing mean as the per-spec benchmark).
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
% Persist forecast results in the layout that table*.m / figure*.m expect.
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
econFactors       = R.econFactors;
pocketIndicesLPM  = R.pocketIndicesLPM;
dependent         = R.dependent;
yAR1              = R.yAR1;
save(fileName, 'yF1Mat','yF2Mat','aMat','pocketIndMat','durationMat', ...
    'r2Mat','integralR2Mat','dateVec','fDMMat','gam1Mat','yF1PMMat', ...
    'r21yMat','r2PMMat','r2GamMat','yF2AltMat','weightMatG', ...
    'yActual','riskFree','windowLengthYears','signRestriction', ...
    'coefRestriction','econFactors','pocketIndicesLPM','dependent','yAR1');
end


function saveOOSResults(fileName, dmMat, cwMat, cwDiffMat, P, F, R, ...
    cumsumPlotIn, cumsumPlotOut, pocketIndExtend)
% Pack and save the variable layout consumed by table3 / tableA.* / figures.
% Several inputs are consumed by `save` by name below.
statMat       = cwMat;
statDiffMat   = cwDiffMat;
econMat       = [P.alphaAnnual', P.tStatMat(1,:)', P.SRAnnual'];
coefMatVol    = P.coefMatVol;    tStatMatVol = P.tStatMatVol;
coefMatMom    = P.coefMatMom;    tStatMatMom = P.tStatMatMom;
coefMat3Fac   = P.coefMat3Fac;   tStatMat3Fac = P.tStatMat3Fac;
portfolioWeight = F.portfolioWeight;
dateVec = R.dateVec;
% dmMat is computed but not consumed by current downstream (statMat = cwMat).
% Include in the save anyway so it is available for ad-hoc queries.
save(fileName, 'statMat','econMat','coefMatVol','tStatMatVol', ...
    'coefMatMom','tStatMatMom','coefMat3Fac','tStatMat3Fac', ...
    'cumsumPlotIn','cumsumPlotOut','pocketIndExtend','dateVec', ...
    'portfolioWeight','statDiffMat','dmMat');
end
