function dailyAssetPricingBootstrap(cfg)
% DAILYASSETPRICINGBOOTSTRAP  Bootstrap on simulated returns from five
% asset-pricing models (Bansal-Yaron, Campbell-Cochrane, Duffie-Tomz,
% Gabaix-Postemski, Wachter, plus a Wachter-no-disasters variant). Drives
% the model-comparison statistics in Table 4.
%
% Inputs (read from cfg.paths.csvSims):
%   BY_sim.csv, CC_sim.csv, DT_sim.csv, GP_sim.csv, W_sim.csv
%
% Outputs (written to cfg.paths.results/assetPricing/):
%   <Model>_asset_pricing_sims_*.mat for each of the six model variants
%
% Tables/figures consumed by:
%   Table 4 (asset-pricing model OOS simulations)
%
% Runtime (B=1000, parallel, 8-core): ~2-4 hours total. The bottleneck
% is `readtable` on the multi-GB simulation CSVs (one per model), which
% does not scale with B; bootstrap compute scales with B but is small
% relative to I/O.
%
% RNG: rng(cfg.rngSeed) is set at function entry; per-replication seeds
% are pre-generated serially before each parfor.

if nargin < 1 || isempty(cfg); cfg = default_config(); end
paths = cfg.paths;
rng(cfg.rngSeed);

apOutDir = fullfile(paths.results, 'assetPricing');
if ~isfolder(apOutDir); mkdir(apOutDir); end

%% Setup
if isfield(cfg, 'bootstrapReps') && ~isempty(cfg.bootstrapReps)
    B = cfg.bootstrapReps;
else
    B = 1000;
end
TMax = 23786;
TVec = [TMax, 15859, 23726];
bufferInd = 2 * 1360;
TMin      = min(TVec) - bufferInd;
weightRestriction = 1;
minLength         = cfg.minPocketDays;

% Each row = (modelName, csv path index inside W_sim.csv if filtered).
% readIndices: [exretCol, predictorCol (raw, then -log applied), rfCol].
% erInd: column of the true expected return.
modelSpecs = {
    'BY',   [12, 5, 3], 9,  false;
    'CC',   [6,  3, 2], 7,  false;
    'DT',   [7, 11, 3], 6,  false;
    'GP',   [6,  3, 2], 5,  false;
    'W',    [9,  7, 6], 8,  false;
    'W_nd', [9,  7, 6], 8,  true;       % filter Z==0 from W_sim.csv
};
restrictVals = [0,0; 1,0; 1,1];

[K1, ~] = oneSidedKernel();

%% Model loop
for fileNum = 1:size(modelSpecs, 1)
    modelName   = modelSpecs{fileNum, 1};
    readIndices = modelSpecs{fileNum, 2};
    filterZ     = modelSpecs{fileNum, 4};
    fprintf('Model %d/6 (%s)\n', fileNum, modelName);

    %% Load and reshape simulation data
    if filterZ
        data = readtable(fullfile(paths.csvSims, 'W_sim.csv'));
        data = data(data.Z == 0, :);
    else
        data = readtable(fullfile(paths.csvSims, [modelName, '_sim.csv']));
    end
    data = data(1:(TMax * B), :);

    exretMat = reshape(data{:, readIndices(1)} - data{:, readIndices(3)} ./ 252, ...
        TMax, B) * 100;
    rvarTemp = movsum((data{:, readIndices(1)} - data{:, readIndices(3)} ./ 252).^2, ...
        [59, 0]);
    predictorMat = NaN(TMax, B, 3);
    predictorMat(:, :, 1) = -log(reshape(data{:, readIndices(2)}, TMax, B));
    predictorMat(:, :, 2) = reshape(data{:, readIndices(3)}, TMax, B) * 100;
    predictorMat(:, :, 3) = reshape(rvarTemp, TMax, B);

    clear data

    %% Restriction loop
    for restrict = 1:size(restrictVals, 1)
        signRestriction = restrictVals(restrict, 1);
        coefRestriction = restrictVals(restrict, 2);
        if ~coefRestriction
            restrictMat = [0, 0, 0];
        else
            restrictMat = [1, -1, 1];
        end

        bootSeeds = randi(2^31 - 1, B, 1);

        econMat       = NaN(B, 3, 3, 2);
        dmMat         = NaN(B, 3, 3);
        cwMat         = NaN(B, 3, 3);
        cwDiffMat     = NaN(B, 3);
        tStatAlphaMat = NaN(B, 3, 2);

        parfor (b = 1:B, useParfor(cfg))
            rng(bootSeeds(b));
            [econSlc, dmSlc, cwSlc, cwDiffSlc, tStatSlc] = runReplication( ...
                b, exretMat, predictorMat, TVec, TMin, K1, restrictMat, ...
                signRestriction, weightRestriction, minLength, cfg);
            econMat(b,:,:,:)     = econSlc;
            dmMat(b,:,:)         = dmSlc;
            cwMat(b,:,:)         = cwSlc;
            cwDiffMat(b,:)       = cwDiffSlc;
            tStatAlphaMat(b,:,:) = tStatSlc;
        end

        outFile = fullfile(apOutDir, sprintf( ...
            '%s_asset_pricing_sims_signRestriction_%d_coefRestriction_%d_25ybandwidth_1yDM_OOS.mat', ...
            modelName, signRestriction, coefRestriction));
        save(outFile, 'dmMat', 'cwMat', 'cwDiffMat', 'econMat', 'tStatAlphaMat');
        fprintf('  wrote %s\n', outFile);
    end

    clear exretMat predictorMat
end
end


% ====================================================================
%  Local helpers — entry-point-specific.
% ====================================================================

function [econSlc, dmSlc, cwSlc, cwDiffSlc, tStatSlc] = runReplication( ...
    b, exretMat, predictorMat, TVec, TMin, K1, restrictMat, ...
    signRestriction, weightRestriction, minLength, cfg)
% One bootstrap replication: 3 predictor variants for the same b'th draw.
econSlc   = NaN(3, 3, 2);
dmSlc     = NaN(3, 3);
cwSlc     = NaN(3, 3);
cwDiffSlc = NaN(1, 3);
tStatSlc  = NaN(3, 2);

for ii = 1:3
    T   = TVec(ii);
    h   = 2.7 * 252 / T;            % famaFrench-style inflated bandwidth
    out = floor(h * T);

    % Slice the b'th replication and stagger for one-step-ahead prediction.
    y   = exretMat(:, b);
    yem = cumsum(y) ./ ((1:numel(y))');
    X   = predictorMat(:, b, ii);
    yem = yem(end-T+1:end);
    y   = y(end-T+1:end);
    X   = X(end-T+1:end);
    yem = yem(1:end-1);
    y   = y(2:end);
    X   = X(1:end-1);
    X   = X ./ std(X);
    T   = T - 1;

    % Kernel forecast + recursive winsorization (qLB=2.5, qUB=97.5 hard-coded
    % in the original — preserved).
    [~, ~, ~, yForecast] = lps1_v2(y, [ones(T,1), X], h, 0, K1, 1.2, ...
        restrictMat(ii));
    yForecast = winsorizeRecursiveFast(yForecast, 2.5, 97.5, out + 1);

    % G-prior shrinkage. With cfg.shrinkageTarget='zero' (default reply
    % setting), gPriorWeights is bit-identical to the original inline loop.
    weightWindow = 252;
    benchmarkForecast = computeBenchmarkForecast(yForecast, yem, cfg);
    [yCombinedG, ~, ~] = gPriorWeightsFast(yForecast, y, yem, ...
        benchmarkForecast, weightWindow, out, cfg.gPrior, 0, 1);

    % Align all series to the OOS sample boundary; rescale to decimals.
    yF1 = yem(2*out:end) ./ 100;
    yF2 = yCombinedG(end-numel(yF1)+1:end) ./ 100;
    if signRestriction
        yF2(yF2 < 0) = 0;
    end
    yActual = y(2*out:end) ./ 100;

    % Smooth squared-error differential and detect pockets.
    fDM = (yActual - yF1).^2 - (yActual - yF2).^2;
    [~, ~, ~, ~, ~, fDMHat] = lps1_v2(fDM, ...
        [ones(size(fDM)), (1:numel(fDM))'], 252/numel(fDM), 0, K1, 1.2, 0);

    pocketIndices = stickyExpPocketFilter(fDM, fDMHat, TMin, minLength);

    % Trim everything to the common sample.
    yF1     = yF1(end-TMin+1:end);
    yF2     = yF2(end-TMin+1:end);
    yActual = yActual(end-TMin+1:end);

    fCW = (yActual - yF1).^2 - ((yActual - yF2).^2 - (yF1 - yF2).^2);
    fDM = (yActual - yF1).^2 - (yActual - yF2).^2;

    dmSlc(ii, :)   = computeDMCWStats(fDM, pocketIndices, 'dm');
    cwFull         = computeDMCWStats(fCW, pocketIndices, 'cw');
    cwSlc(ii, :)   = cwFull(1:3);
    cwDiffSlc(ii)  = cwFull(4);

    % Portfolio metrics: 2 specs (TVC vs prevailing-mean baseline).
    [~, yPocketTime] = computeGomezCramWeights( ...
        [yF2, yF1], yActual, 0, 2, weightRestriction);

    [alphaA, srA, deltaA, tStatA] = computeUtilityPanel( ...
        yPocketTime, yActual, 3);
    econSlc(ii, :, :) = [alphaA; srA; deltaA];
    tStatSlc(ii, :)   = tStatA;
end
end


function pocketIndices = stickyExpPocketFilter(fDM, fDMHat, TMin, minLength)
% Custom min-length pocket filter (NOT detectPocketsFromSED). Clips the
% FIRST minLength entries of each pocket but retains a final-period
% pocket fully because the loop does not enter the cleanup branch when
% pocketIndices(TMin) is true. Matches the original convention shared by
% dailyAssetPricingBootstrap and stickyExpectationsSim.
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


function [alphaA, srA, deltaA, tStatA] = computeUtilityPanel(yPocketTime, yActual, gam)
% Per-spec annualized alpha, Sharpe, per-day Delta certainty-equivalent
% return (per-period fsolve to preserve original calibration), and the
% alpha t-stat. Spec 2 (prevailing-mean baseline) overrides Sharpe with
% the unconditional market Sharpe.
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
    if jj == 2
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
