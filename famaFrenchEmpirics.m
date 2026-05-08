function famaFrenchEmpirics(cfg, ffFactor)
% FAMAFRENCHEMPIRICS  Fama-French SMB and/or HML factor portfolio analysis.
% Default produces both panels of Table A.4 in a single call.
%
% Inputs:
%   cfg       config struct (default_config())
%   ffFactor  (optional) 'SMB', 'HML', or a cell array such as
%             {'SMB','HML'}. Defaults to {'SMB','HML'} (full Table A.4).
%
% Reads from cfg.paths.data:
%   Daily_Predictors.xlsx, pcPredictor.mat
%
% Outputs (written to cfg.paths.results/famaFrench/):
%   forecastResultsSMB_*.mat  /  forecastResultsHML_*.mat
%   OOSResultsSMB_*.mat       /  OOSResultsHML_*.mat
%
% Tables/figures consumed by:
%   Table A.4
%
% Runtime: ~1-2 minutes for both factors.
%
% Convention notes:
%   - Single specification: signRestriction=1, coefRestriction=0,
%     windowLengthYears=2.5 (with the 2.7y inflated bandwidth in the
%     kernel call), dmLengthYears=1.
%   - SED smoothing uses a k=2 design with a time trend (sedTimeTrend=true).
%     gam1Mat captures the time-trend coefficient (gamDM(:,2)).
%   - Daily windows for the portfolio factors but the additive adjTerm
%     formula (freq='famaFrench' in constructPortfolioFactors).

if nargin < 1 || isempty(cfg);      cfg = default_config();      end
if nargin < 2 || isempty(ffFactor); ffFactor = {'SMB', 'HML'};   end
if ischar(ffFactor) || isstring(ffFactor)
    ffFactor = {char(ffFactor)};
end

for kFactor = 1:numel(ffFactor)
    fprintf('famaFrenchEmpirics: factor %d/%d (%s)\n', ...
        kFactor, numel(ffFactor), ffFactor{kFactor});
    runFactor(cfg, ffFactor{kFactor});
end
end


function runFactor(cfg, ffFactor)
% Per-factor body (the original famaFrenchEmpirics implementation).
paths = cfg.paths;
resultsFolder = baselineResultsFolder(paths, cfg, 'famaFrench');

%% Setup
varNames           = {'dp','tbl','tsp','rvar','mv','pc'};
nVars              = numel(varNames);
signRestriction    = 1;
coefRestriction    = 0;
windowLengthYears  = 2.5;          % nominal label
windowLengthYearsK = 2.7;          % inflated bandwidth used in the kernel
dmLengthYears      = 1;
benchmarkFlag      = 'pm';
kernelMethod       = '1S';
restrictMat        = [0, 0, 0, 0, 0, 0];

[eLabel, dmLabel, rLabel] = buildSpecLabels( ...
    signRestriction, coefRestriction, windowLengthYears, dmLengthYears);
fileName  = fullfile(resultsFolder, sprintf( ...
    'forecastResults%s%s%s%s_%s_%s.mat', ffFactor, rLabel, eLabel, ...
    dmLabel, kernelMethod, benchmarkFlag));
fileName2 = fullfile(resultsFolder, sprintf( ...
    'OOSResults%s%s%s%s.mat', ffFactor, rLabel, eLabel, dmLabel));

%% Data load and result preallocation
data = readtable(fullfile(paths.data, 'Daily_Predictors.xlsx'));
[hPre, ~] = computeKernelBandwidth(windowLengthYearsK, numel(data.dp) - 1, 'daily');
bufferInd = 2 * floor(hPre * (numel(data.dp) - 1));
TMax      = size(data, 1) - bufferInd;
R         = preallocateResults(TMax, nVars, 100);
R.econFactors = [data.recession(end-TMax+1:end), data.bwindex(end-TMax+1:end), ...
                 data.lfactor(end-TMax+1:end),  data.cgcoef(end-TMax+1:end)];

predictorMat     = [data.dp, data.tbl, data.tsp, data.rvar];
pocketIndicesLPM = [];

switch ffFactor
    case 'SMB'; dependentFull = data.smb;
    case 'HML'; dependentFull = data.hml;
    otherwise
        error('famaFrenchEmpirics:badFactor', ...
            'ffFactor must be ''SMB'' or ''HML''');
end

%% Per-variable estimation
for varNum = 1:nVars
    varName = varNames{varNum};

    % Predictor selection.
    if strcmp(varName, 'pc')
        predictorMat = (predictorMat - mean(predictorMat, 'omitnan')) ./ ...
            std(predictorMat, 'omitnan');
        predictors = selectPredictor('pc', data, paths);
    elseif strcmp(varName, 'mv')
        predictors = predictorMat;
    else
        predictors = selectPredictor(varName, data, paths);
    end

    % Trim, stagger, and build the prevailing-mean naive forecast.
    S = prepareDataForEstimation(varName, predictors, dependentFull, ...
        data.rf, data.Date, predictorMat);

    % Kernel + estimation hyperparameters. famaFrench uses an inflated
    % bandwidth (2.7 vs 2.5) to compensate for the Epanechnikov kernel's
    % effective sample-size reduction.
    [h, out] = computeKernelBandwidth(windowLengthYearsK, S.T, 'daily');
    [Ke1, ~] = oneSidedKernel();
    hSpec = struct('h', h, 'out', out, 'order', 0, 'Ke1', Ke1, ...
        'dmLengthYears', dmLengthYears, 'obsPerYear', 252, ...
        'minRT', cfg.minPocketDays, ...
        'weightWindow', round(cfg.weightWindowYears * 252), ...
        'sedTimeTrend', true);                              % famaFrench-specific

    % Run kernel forecasting + G-prior + pocket detection.
    Rv = estimateVariableForecasts(S, hSpec, restrictMat(varNum), ...
        signRestriction, benchmarkFlag, cfg);

    % LPM pocket detection (only for dp; consumed by combo + figure).
    if varNum == 1
        fDM2 = (Rv.yActual - Rv.yF1PM).^2 - (Rv.yActual - Rv.yF1).^2;
        [~, ~, ~, ~, ~, fDMHat2] = lps1_v2(fDM2, ...
            [ones(size(fDM2)), (1:numel(fDM2))'], ...
            dmLengthYears * 252 / numel(fDM2), 0, Ke1, [], 0);
        pocketIndicesLPM = detectPocketsFromSED(fDMHat2, ...
            numel(Rv.r2LocalPM), cfg.minPocketDays);
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

% forecastResults persists the un-trimmed, un-scaled forecast matrices.
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

% Append the *S (trimmed/scaled) panels for figure consumers.
yF2MatS       = [R.yF2Mat, yComb1, yComb2, yComb3];
pocketIndMatS = [R.pocketIndMat, combPocket, combPocket, combPocket];
dateVecS      = R.dateVec;
save(fileName, 'yF2MatS', 'pocketIndMatS', 'dateVecS', '-append');

%% Portfolio metrics + performance regressions
yForecastMat = [R.yF2Mat(:,1:6), yComb1, yComb2, yComb3, ...
                R.yF1Mat(:,1), R.yF1PMMat(:,1)];
F = constructPortfolioFactors(yForecastMat, R.yActual, R.riskFree, ...
    'famaFrench', true, cfg.adjCostBps);
P = computePerformanceMetrics(F.yPocketTime, R.yActual, ...
    F.ySigFactor, F.yMomFactor, 'famaFrench');

%% Forecast differential tests
[dmMat, cwMat, cwDiffMat] = computeForecastTests(fDMCell, fCWCell, pIndCell);
% Match the archive convention: only univariate (1-6) and comb3 (9)
% carry an in-vs-out difference stat. comb1, comb2, and LPM are NaN
% in table A.4 of the published reply.
cwDiffMat([7, 8, 10]) = NaN;

%% Final OOS save
saveOOSResults(fileName2, dmMat, cwMat, cwDiffMat, P);
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
R.durationMat   = NaN(maxPockets, nVars);
R.integralR2Mat = NaN(maxPockets, nVars);
R.fDMMat        = NaN(TMax, nVars);
end


function R = packVariableResults(R, Rv, varNum, out)
% NOTE: famaFrench stores gamDM(:,2) (time-trend coefficient) in gam1Mat,
% in contrast to daily/monthly which store gamDM(:,1).
nF1   = numel(Rv.yF1);    R.yF1Mat(end-nF1+1:end, varNum)         = Rv.yF1;
nF1PM = numel(Rv.yF1PM);  R.yF1PMMat(end-nF1PM+1:end, varNum)     = Rv.yF1PM;
nF2   = numel(Rv.yF2);    R.yF2Mat(end-nF2+1:end, varNum)         = Rv.yF2;
nA    = numel(Rv.yF2Alt); R.yF2AltMat(end-nA+1:end, varNum)       = Rv.yF2Alt;
aTrim = Rv.a(out+1:end);
R.aMat(end-numel(aTrim)+1:end, varNum)        = aTrim;
R.gam1Mat(end-size(Rv.gamDM,1)+1:end, varNum) = Rv.gamDM(:,2);
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
% appended after yComb1 is computed. Total = 10 specs (no erL).
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


function renderPerformanceTable(P)
tableValues = [P.alphaAnnual', P.coefMat(2,:)', P.appRatioAnnual', ...
               P.SRAnnual', NaN(numel(P.alphaAnnual), 1)];
rowLabels   = {'dp','tbl','tsp','rvar','mv','pc', ...
               'comb1','comb2','comb3','lpm','pm'};
colLabels   = {'\textbf{Variables}','$\hat{\alpha}$ \textbf{(annualized)}', ...
    '$\hat{\beta}$','\textbf{Appraisal Ratio (annualized)}', ...
    '\textbf{Sharpe Ratio (annualized)}', ...
    '$\hat{\Delta}$ \textbf{with} $\gamma=3$'};
opts.title   = 'Out-of-sample measures of forecasting performance';
opts.caption = '';
LatexTableFull(tableValues, colLabels, rowLabels, '9.2f', [], 0, opts);
end


function renderDMCWTable(dmMat, cwMat)
tableValues = [dmMat; cwMat];
rowLabels   = repmat({'dp','tbl','tsp','rvar','mv','pc', ...
    'comb1','comb2','comb3','lpm'}, 1, 2);
colLabels   = {'\textbf{Variables}','\textbf{Full sample}', ...
    'In-pocket (real time)','Out-of-pocket (real time)'};
opts.title   = 'Out-of-sample measures of forecasting performance';
opts.caption = '';
LatexTableFull(tableValues, colLabels, rowLabels, '9.2f', [], 0, opts);
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
yActual      = R.yActual;
riskFree     = R.riskFree;
windowLengthYears = R.windowLengthYears;
signRestriction   = R.signRestriction;
coefRestriction   = R.coefRestriction;
econFactors       = R.econFactors;
pocketIndicesLPM  = R.pocketIndicesLPM;
save(fileName, 'yF1Mat','yF2Mat','aMat','pocketIndMat','durationMat', ...
    'r2Mat','integralR2Mat','dateVec','fDMMat','gam1Mat','yF1PMMat', ...
    'r21yMat','r2PMMat','r2GamMat','yF2AltMat', ...
    'yActual','riskFree','windowLengthYears','signRestriction', ...
    'coefRestriction','econFactors','pocketIndicesLPM');
end


function saveOOSResults(fileName, dmMat, cwMat, cwDiffMat, P) %#ok<INUSL>
% Layout matches the original: statMat, econMat, statDiffMat, plus the
% multi-factor regression panels. NO cumsumPlotIn/Out or portfolioWeight
% in the famaFrench OOS save (those are daily/monthly only). The dmMat
% input is consumed by `save` by name; INUSL suppression matches.
statMat       = cwMat;
statDiffMat   = cwDiffMat;
econMat       = [P.alphaAnnual', P.tStatMat(1,:)', P.SRAnnual'];
coefMatVol    = P.coefMatVol;    tStatMatVol = P.tStatMatVol;
coefMatMom    = P.coefMatMom;    tStatMatMom = P.tStatMatMom;
coefMat3Fac   = P.coefMat3Fac;   tStatMat3Fac = P.tStatMat3Fac;
save(fileName, 'statMat','econMat','coefMatVol','tStatMatVol', ...
    'coefMatMom','tStatMatMom','coefMat3Fac','tStatMat3Fac', ...
    'statDiffMat','dmMat');
end
