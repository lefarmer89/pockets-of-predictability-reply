function dailyBootstrap(cfg)
% DAILYBOOTSTRAP  Bootstrap replications for the in-pocket and
% out-of-pocket Clark-West / Diebold-Mariano statistics. Drives the
% significance markers in Table 2 of the reply and feeds signifCell.mat.
%
% Inputs (read from cfg.paths.data):
%   Daily_Predictors.xlsx
%
% Outputs (written to cfg.paths.results/bootstrap/):
%   bootstrapResultsOOS_*.mat for each (restriction, spec, variable) tuple
%
% Tables/figures consumed by:
%   Table 2 (significance markers on pocket statistics)
%
% Runtime (B=1000, parallel, 8-core): ~3-6 hours total across the three
% restriction specs. Smoke at B=5 runs in ~90 s.
%
% RNG: rng(cfg.rngSeed) is set at function entry; per-replication seeds
% are pre-generated serially before each parfor so worker count does not
% change results.

if nargin < 1 || isempty(cfg); cfg = default_config(); end
paths = cfg.paths;
rng(cfg.rngSeed);

bootDir = fullfile(paths.results, 'bootstrap');
if ~isfolder(bootDir); mkdir(bootDir); end

%% Setup
data = readtable(fullfile(paths.data, 'Daily_Predictors.xlsx'));

restrictVals = [0,0; 1,0; 1,1];
% Each spec column is [hetType; bootType].
specNames    = {'none','egarch','none'; 'iid','iid','stationary'};
varList      = {'dp','tbl','tsp','rvar','mv'};
if isfield(cfg, 'bootstrapReps') && ~isempty(cfg.bootstrapReps)
    B = cfg.bootstrapReps;
else
    B = 1000;
end
hConst       = 2.7;                              % 2.5y inflated bandwidth
bufferInd    = 2 * 680;
TMin         = sum(~isnan(data.tsp)) - bufferInd;
weightRestriction = 1;
minLength    = 21;

[K, ~] = oneSidedKernel();

%% Restriction loop
for restrict = 1:size(restrictVals, 1)
    signRestriction = restrictVals(restrict, 1);
    coefRestriction = restrictVals(restrict, 2);
    if ~coefRestriction
        restrictMat = [0, 0, 0, 0, 0];
    else
        restrictMat = [1, -1, 1, 1, 0];
    end

    %% Bootstrap-spec loop (iid+none, iid+egarch, stationary+none)
    for spec = 1:size(specNames, 2)
        fprintf('Restriction %d/%d  Bootstrap spec %d/3\n', ...
            restrict, size(restrictVals,1), spec);
        hetType  = specNames{1, spec};
        bootType = specNames{2, spec};

        %% Variable loop
        for varNum = 1:numel(varList)
            varName = varList{varNum};
            predictors = selectVariablePredictors(varName, data);

            % Trim, stagger, prepare X and y for bootstrap fitting.
            trim = max(sum(isnan(predictors)));
            X = predictors(trim+1:end, :);
            y = data.exret(trim+1:end);
            T = size(X, 1) - 1;

            % Pre-fit OLS, AR(1) on each predictor, optional EGARCH on
            % both residual streams, residual ECDFs, and the optimal
            % stationary-bootstrap block length.
            calib = calibrateBootstrap(y, X, T, hetType);

            % Bandwidth and parfor seed setup.
            h   = hConst * 252 / T;
            out = floor(h * T);
            bootSeeds = randi(2^31 - 1, B, 1);

            opts = struct( ...
                'varName',           varName, ...
                'corrSetting',       'default', ...
                'bootType',          bootType, ...
                'hetType',           hetType, ...
                'restrictVal',       restrictMat(varNum), ...
                'signRestriction',   signRestriction, ...
                'weightRestriction', weightRestriction, ...
                'h',                 h, ...
                'out',               out, ...
                'T',                 T, ...
                'TMin',              TMin, ...
                'minLength',         minLength, ...
                'K',                 K);

            %% Bootstrap parfor
            dmMat         = NaN(B, 3);
            cwMat         = NaN(B, 3);
            cwDiffMat     = NaN(B, 1);
            econMat       = NaN(B, 3, 2);
            tStatAlphaMat = NaN(B, 2);
            blockLengths  = NaN(B, 500);
            integralR2    = NaN(B, 500);

            parfor (b = 1:B, useParfor(cfg))
                rng(bootSeeds(b));
                R = runReplication(b, calib, opts, cfg);
                dmMat(b,:)         = R.dmRow;
                cwMat(b,:)         = R.cwRow;
                cwDiffMat(b)       = R.cwDiff;
                econMat(b,:,:)     = R.econRow;
                tStatAlphaMat(b,:) = R.tStatRow;
                blockLengths(b,:)  = R.blockLengths;
                integralR2(b,:)    = R.integralR2;
            end

            outFile = fullfile(bootDir, sprintf( ...
                'bootstrapResultsOOS_signRestriction_%d_coefRestriction_%d_%s_%sbs_%s_25y_0.mat', ...
                signRestriction, coefRestriction, varName, bootType, hetType));
            save(outFile, 'dmMat','cwMat','cwDiffMat','econMat', ...
                'tStatAlphaMat','blockLengths','integralR2');
            fprintf('  wrote %s\n', outFile);
        end
    end
end
end


% ====================================================================
%  Local helpers — entry-point-specific.
% ====================================================================

function predictors = selectVariablePredictors(varName, data)
% dailyBootstrap-specific predictor dispatcher (5-variable subset; mv is
% the 4-column multivariate stack). selectPredictor in subroutines/ is
% empirics-style and includes erL/pc which dailyBootstrap does not use.
switch varName
    case 'dp';   predictors = data.dp;
    case 'tbl';  predictors = data.tbl;
    case 'tsp';  predictors = data.tsp;
    case 'rvar'; predictors = data.rvar;
    case 'mv';   predictors = [data.dp, data.tbl, data.tsp, data.rvar];
    otherwise
        error('dailyBootstrap:badVar', 'unknown varName %s', varName);
end
end


function calib = calibrateBootstrap(y, X, T, hetType)
% Pre-fit the predictive regression, AR(1) on each predictor, optional
% EGARCH on both residual streams, residual ECDFs, and the optimal
% stationary-bootstrap block length. Captures everything the per-
% replication body needs from the actual data.
nx = size(X, 2);

% Predictive regression on full sample.
calib.gamma0 = [ones(T,1), X(1:end-1,:)] \ y(2:end);
eHatR = y(2:end) - [ones(T,1), X(1:end-1,:)] * calib.gamma0;

% Variance model on the predictive residuals.
[calib.uHatR, calib.htR, calib.egR] = fitVariance(eHatR, hetType);
calib.eStdR = std(eHatR);

% AR(1) on each predictor.
calib.rho       = NaN(2, nx);
calib.uHatX     = NaN(T, nx);
calib.htX       = NaN(T, nx);
calib.eStdX     = NaN(1, nx);
calib.egX(nx,1) = struct('constant',NaN,'garch',NaN,'arch',NaN,'leverage',NaN,'meanNorm',NaN);
uXf = NaN(T+1, nx);
uXx = NaN(T+1, nx);

for kk = 1:nx
    calib.rho(:, kk) = [ones(T,1), X(1:end-1,kk)] \ X(2:end,kk);
    eHatX = X(2:end, kk) - [ones(T,1), X(1:end-1,kk)] * calib.rho(:,kk);
    [calib.uHatX(:,kk), calib.htX(:,kk), calib.egX(kk)] = ...
        fitVariance(eHatX, hetType);
    calib.eStdX(kk) = std(eHatX);
    [tempXf, tempXx] = ecdf(calib.uHatX(:,kk));
    uXf(end-numel(tempXf)+1:end, kk) = tempXf;
    uXx(end-numel(tempXx)+1:end, kk) = tempXx;
end

[calib.uRf, calib.uRx] = ecdf(calib.uHatR);
calib.uXf = uXf;
calib.uXx = uXx;

% Optimal block length for stationary bootstrap.
calib.optBlock     = ceil(opt_block_length([calib.uHatR, calib.uHatX]));
calib.optBlock     = calib.optBlock(1, 1);
calib.residExtend  = [calib.uHatR, calib.uHatX; calib.uHatR, calib.uHatX];
calib.X            = X;
calib.T            = T;
calib.nx           = nx;
calib.garchLag     = 1;
calib.archLag      = 1;
end


function [uHat, ht, params] = fitVariance(eHat, hetType)
% Fit either a t-distributed EGARCH(1,1) model or no variance model and
% return the standardized residuals and conditional variances. Captures
% the per-residual-stream calibration the bootstrap evolution needs.
switch hetType
    case 'egarch'
        mdl = egarch(1, 1); mdl.Distribution = 't';
        m   = estimate(mdl, eHat);
        ht  = infer(m, eHat);
        uHat = eHat ./ sqrt(ht);
        params.constant = m.Constant;
        params.garch    = cell2mat(m.GARCH)';
        params.arch     = cell2mat(m.ARCH)';
        params.leverage = cell2mat(m.Leverage)';
        nu = m.Distribution.DoF;
        params.meanNorm = sqrt((nu-2)/pi) * gamma((nu-1)/2) / gamma(nu/2);
    case 'none'
        uHat   = eHat ./ std(eHat);
        ht     = NaN(size(eHat));
        params = struct('constant',NaN,'garch',NaN,'arch',NaN, ...
                        'leverage',NaN,'meanNorm',NaN);
end
end


function R = runReplication(~, calib, opts, cfg)
% One bootstrap replication: generate bootstrap sample, fit TVC model,
% G-prior shrinkage, pocket detection, DM/CW stats, portfolio + Delta.
T = calib.T; out = opts.out; TMin = opts.TMin;

% --- Generate bootstrap sample (yb, Xb)
[yb, Xb] = generateBootstrapSample(calib, opts);

% Trim and stagger for one-step-ahead prediction.
yemb = cumsum(yb(1:end-1)); yemb = yemb ./ ((1:numel(yemb))');
yb   = yb(2:end);
Xb   = Xb(1:end-1, :);

% Optional centering: subtract the predictive-regression yhat.
sFull = regstats2Fast(yb, Xb, 'linear', {'beta','hac','adjrsquare','yhat'});
yb    = yb - sFull.yhat;

% --- TVC kernel forecast + winsorize + G-prior shrinkage
[~, ~, r2Local, yForecast] = lps1_v2(yb, [ones(T,1), Xb], opts.h, ...
    0, opts.K, 1.2, opts.restrictVal);
r2Local   = [0; r2Local(1:end-1)];
yForecast = winsorizeRecursiveFast(yForecast, 2.5, 97.5, out + 1);

weightWindow = 252;
benchmarkForecast = computeBenchmarkForecast(yForecast, yemb, cfg);
[yCombinedG, ~, ~] = gPriorWeightsFast(yForecast, yb, yemb, ...
    benchmarkForecast, weightWindow, out, cfg.gPrior, 0, 1);

% Align series at the OOS sample boundary; rescale to decimals.
yF1 = yemb(2*out:end) ./ 100;
yF2 = yCombinedG(end-numel(yF1)+1:end) ./ 100;
if opts.signRestriction
    yF2(yF2 < 0) = 0;
end
yActual = yb(2*out:end) ./ 100;

% --- Pocket detection (custom min-length filter; same as stickyExp /
%      dailyAssetPricingBootstrap convention).
fDM = (yActual - yF1).^2 - (yActual - yF2).^2;
[~, ~, ~, ~, ~, fDMHat] = lps1_v2(fDM, ...
    [ones(size(fDM)), (1:numel(fDM))'], 252/numel(fDM), 0, opts.K, 1.2, 0);
pocketIndices = stickyExpPocketFilter(fDM, fDMHat, TMin, opts.minLength);
yF1 = yF1(end-TMin+1:end); yF2 = yF2(end-TMin+1:end);
yActual = yActual(end-TMin+1:end);
r2Local = r2Local(end-TMin+1:end);

% --- Block boundaries for diagnostic block-length / integral-R^2 panels.
[blockLengths, integralR2] = computePocketBlocks(pocketIndices, r2Local, 500);

% --- DM / CW test stats.
fCW = (yActual - yF1).^2 - ((yActual - yF2).^2 - (yF1 - yF2).^2);
fDM = (yActual - yF1).^2 - (yActual - yF2).^2;
R.dmRow  = computeDMCWStats(fDM, pocketIndices, 'dm');
cwFull   = computeDMCWStats(fCW, pocketIndices, 'cw');
R.cwRow  = cwFull(1:3);
R.cwDiff = cwFull(4);

% --- Portfolio metrics: replace TVC with PM benchmark out-of-pocket
%      before constructing yPocketTime (matches original line 573).
yF2(~pocketIndices) = yF1(~pocketIndices);
[~, yPocketTime] = computeGomezCramWeights( ...
    [yF2, yF1], yActual, 0, 2, opts.weightRestriction);
[alphaA, srA, deltaA, tStatA] = computeUtilityPanel(yPocketTime, yActual, 3);
R.econRow      = [alphaA; srA; deltaA];
R.tStatRow     = tStatA;
R.blockLengths = blockLengths;
R.integralR2   = integralR2;
end


function [yb, Xb] = generateBootstrapSample(calib, opts)
% Build (yb, Xb) of length T+1 by drawing residuals (iid or stationary
% block) and either rescaling them by std (hetType='none') or evolving
% them through the EGARCH(1,1) recursion (hetType='egarch').
T  = calib.T; nx = calib.nx;
yb = zeros(T + 1, 1);
startInd = randi(T, 1);
Xb = [calib.X(startInd, :); zeros(T, nx)];

% Residual draws.
[ubR, ubX] = drawResiduals(calib, opts.bootType, opts.varName, ...
    opts.corrSetting, T);

switch opts.hetType
    case 'egarch'
        p = calib.garchLag; q = calib.archLag;
        mpq = max(p, q);
        startInd = randi(T, mpq, 1);
        htbR = [calib.htR(startInd);          zeros(T - mpq + 1, 1)];
        htbX = [calib.htX(startInd, :);       zeros(T - mpq + 1, nx)];
        ebR  = [ubR(1:mpq) .* sqrt(htbR(1:mpq)); zeros(T - mpq + 1, 1)];
        ebX  = [ubX(1:mpq, :) .* sqrt(htbX(1:mpq, :)); zeros(T - mpq + 1, nx)];

        for t = mpq:T
            htbR(t+1) = exp(calib.egR.constant + ...
                sum(calib.egR.garch .* log(htbR(t:-1:t-p+1))) + ...
                sum(calib.egR.arch .* (abs(ebR(t:-1:t-q+1)) ./ sqrt(htbR(t:-1:t-q+1)) - calib.egR.meanNorm)) + ...
                sum(calib.egR.leverage .* ebR(t:-1:t-q+1) ./ sqrt(htbR(t:-1:t-q+1))));
            xConst = vertcat(calib.egX.constant);
            xGarch = vertcat(calib.egX.garch);
            xArch  = vertcat(calib.egX.arch);
            xLev   = vertcat(calib.egX.leverage);
            xMnorm = vertcat(calib.egX.meanNorm);
            htbX(t+1, :) = exp(xConst' + ...
                sum(xGarch' .* log(htbX(t:-1:t-p+1, :)), 1) + ...
                sum(xArch' .* (abs(ebX(t:-1:t-q+1, :)) ./ sqrt(htbX(t:-1:t-q+1, :)) - xMnorm'), 1) + ...
                sum(xLev' .* ebX(t:-1:t-q+1, :) ./ sqrt(htbX(t:-1:t-q+1, :)), 1));
            ebR(t+1)    = sqrt(htbR(t+1)) * ubR(t+1);
            ebX(t+1, :) = sqrt(htbX(t+1, :)) .* ubX(t+1, :);

            yb(t+1)     = [1, Xb(t, :)] * calib.gamma0 + ebR(t+1);
            Xb(t+1, :)  = sum([ones(1, nx); Xb(t, :)] .* calib.rho) + ebX(t+1, :);
        end

    case 'none'
        ebR = ubR * calib.eStdR;
        ebX = ubX .* calib.eStdX;
        for t = 1:T
            yb(t+1)    = [1, Xb(t, :)] * calib.gamma0 + ebR(t+1);
            Xb(t+1, :) = sum([ones(1, nx); Xb(t, :)] .* calib.rho) + ebX(t+1, :);
        end
end
end


function [ubR, ubX] = drawResiduals(calib, bootType, varName, corrSetting, T)
% iid or stationary block bootstrap draws; optional Stambaugh correlation
% injection on dp.
useStambaugh = strcmp(varName, 'dp') && ~isequal(corrSetting, 'default');

switch bootType
    case 'iid'
        if useStambaugh
            simNorm = mvnrnd([0, 0], [1, corrSetting; corrSetting, 1], T+1);
            simUnif = normcdf(simNorm);
            posR = sum(simUnif(:, 1)' <= calib.uRf);
            posX = sum(simUnif(:, 2)' <= calib.uXf);
            ubR = calib.uRx(posR);
            ubX = calib.uXx(posX);
        else
            ind = randsample(1:T, T+1, 'true');
            ubR = calib.uHatR(ind);
            ubX = calib.uHatX(ind, :);
        end
    case 'stationary'
        [~, ind] = stationary_bootstrap((1:T+1)', 1, calib.optBlock);
        if useStambaugh
            ubR = calib.residExtend(ind, 1);
            ubX = corrSetting * ubR + sqrt(1 - corrSetting^2) * randn(T+1, 1);
        else
            temp = calib.residExtend(ind, :);
            ubR  = temp(:, 1);
            ubX  = temp(:, 2:end);
        end
end
end


function pocketIndices = stickyExpPocketFilter(fDM, fDMHat, TMin, minLength)
% Custom min-length pocket filter (NOT detectPocketsFromSED). Clips the
% FIRST minLength entries of each pocket but retains a final-period
% pocket fully because the loop does not enter the cleanup branch when
% pocketIndices(TMin) is true. Same convention as stickyExpectationsSim.
pocketIndices = logical([zeros(numel(fDM) - numel(fDMHat) + 1, 1); ...
                         fDMHat(1:end-1) > 0]);
pocketIndices = pocketIndices(end-TMin+1:end);

curLength = pocketIndices(1);
for t = 2:TMin
    if pocketIndices(t)
        curLength = curLength + 1;
    elseif ~pocketIndices(t) || t == TMin
        pocketIndices(t-curLength : t-curLength + min(curLength, minLength) - 1) = false;
        curLength = 0;
    end
end
end


function [blockLengths, integralR2] = computePocketBlocks(pocketIndices, r2Local, maxBlocks)
% Identify pocket runs in pocketIndices and report each pocket's length
% and integral R^2. Pads to maxBlocks with NaN.
pocketStart = find([0; diff(pocketIndices)] == 1);
if pocketIndices(1)
    pocketStart = [1; pocketStart];
end
pocketEnd = find([0; diff(pocketIndices)] == -1);
if pocketIndices(end)
    pocketEnd = [pocketEnd; numel(pocketIndices)];
end

blockLengths = NaN(1, maxBlocks);
integralR2   = NaN(1, maxBlocks);
for np = 1:numel(pocketStart)
    integralR2(np)   = sum(r2Local(pocketStart(np):pocketEnd(np)));
    blockLengths(np) = pocketEnd(np) - pocketStart(np) + 1;
end
end


function [alphaA, srA, deltaA, tStatA] = computeUtilityPanel(yPocketTime, yActual, gam)
% Per-spec annualized alpha, Sharpe, per-day Delta certainty-equivalent
% return (per-period fsolve to preserve original calibration), and the
% alpha t-stat. Spec 2 (PM baseline) uses the unconditional market Sharpe.
S = size(yPocketTime, 2);
alphaA = NaN(1, S); srA = NaN(1, S); deltaA = NaN(1, S); tStatA = NaN(1, S);
SRmkt    = mean(yActual, 'omitnan') / std(yActual, 'omitnan');
SRmktAnn = SRmkt * sqrt(252);

opts = optimoptions(@fsolve, 'Display', 'off');

for jj = 1:S
    s = regstats2Fast(yPocketTime(:, jj), yActual);
    if isfield(s, 'beta'); coef = s.beta; else; coef = [NaN; NaN]; end
    if isfield(s, 'mse');  mseV = s.mse;  else; mseV = NaN;        end
    if isfield(s, 'hac') && isfield(s.hac, 't')
        tStatA(jj) = s.hac.t(1);
    end

    alphaA(jj) = coef(1) * 252 * 100;
    srA(jj)    = sqrt(SRmkt^2 + (coef(1)^2) / mseV) * sqrt(252);
    if jj == S
        srA(jj) = SRmktAnn;
    end

    DeltaDay = NaN(numel(yActual), 1);
    for t = 1:numel(yActual)
        if ~isnan(yPocketTime(t, jj))
            objectiveFun = @(D) ((yPocketTime(t,jj) - D) - ...
                gam/(2*(1+gam)) * (yPocketTime(t,jj) - D).^2) - ...
                (yActual(t) - gam/(2*(1+gam)) * (yActual(t).^2));
            DeltaDay(t) = fsolve(objectiveFun, 0, opts);
        end
    end
    deltaA(jj) = mean(DeltaDay, 'omitnan') * 252 * 10000;
end
end
