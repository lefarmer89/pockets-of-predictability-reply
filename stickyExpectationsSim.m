function stickyExpectationsSim(cfg)
% STICKYEXPECTATIONSSIM  Run the pocket-detection / G-prior shrinkage
% pipeline on simulated paths from three asset-pricing models (sticky
% expectations, rational expectations, rational recalibrated) and write
% the bootstrapped distributions used by Table 5 of the reply.
%
% Inputs (read from cfg.paths.results/simulatedPaths/):
%   sticky_sims.mat, rational_sims.mat, rational_recalibrated_sims.mat
%
% Outputs (written to cfg.paths.results/assetPricing/):
%   SE_asset_pricing_sims_*.mat
%   RE_asset_pricing_sims_*.mat
%   RE_Recalibrated_asset_pricing_sims_*.mat
%
% Tables/figures consumed by:
%   Table 5 (sticky expectations vs rational benchmark)
%
% Runtime (B=500 replications, 3 models, parallel): ~1-2 hours total
% on an 8-core machine. Smoke at B=5 runs in ~35 s.

if nargin < 1 || isempty(cfg); cfg = default_config(); end
paths = cfg.paths;

apOutDir = fullfile(paths.results, 'assetPricing');
if ~isfolder(apOutDir); mkdir(apOutDir); end

%% Setup
modelFiles = {
    'SE',              'sticky_sims.mat';
    'RE',              'rational_sims.mat';
    'RE_Recalibrated', 'rational_recalibrated_sims.mat';
};
if isfield(cfg, 'bootstrapReps') && ~isempty(cfg.bootstrapReps)
    B = cfg.bootstrapReps;
else
    B = 500;
end
TMax      = 23786;
TVec      = [TMax, 15859, 23726];
bufferInd = 2 * 680;
TMin      = min(TVec) - bufferInd;

% The reply's Table 5 uses only the unrestricted spec
% (signRestriction=0, coefRestriction=0).
restrictMat       = [0, 0, 0];
signRestriction   = 0;
coefRestriction   = 0;
weightRestriction = 1;

[K1, ~] = oneSidedKernel();

%% Model loop
for fileNum = 1:size(modelFiles, 1)
    modelName = modelFiles{fileNum, 1};
    modelMat  = modelFiles{fileNum, 2};
    fprintf('Model %s (file %d/3)\n', modelName, fileNum);

    sims = load(fullfile(paths.results, 'simulatedPaths', modelMat));
    raw  = [sims.rSim, sims.dpSim, sims.rfSim, sims.rvarSim, sims.rpSim];
    raw  = raw(1:(TMax * B), :);

    exretMat     = reshape(raw(:,1) - raw(:,3), TMax, B);
    predictorMat = NaN(TMax, B, 3);
    predictorMat(:,:,1) = reshape(raw(:,2), TMax, B);
    predictorMat(:,:,2) = reshape(raw(:,3), TMax, B);
    predictorMat(:,:,3) = reshape(raw(:,4), TMax, B);

    %% Bootstrap loop
    econMat       = NaN(B, 3, 3, 2);
    dmMat         = NaN(B, 3, 3);
    cwMat         = NaN(B, 3, 3);
    cwDiffMat     = NaN(B, 3);
    tStatAlphaMat = NaN(B, 3, 2);

    parfor (b = 1:B, useParfor(cfg))
        [econSlc, dmSlc, cwSlc, cwDiffSlc, tStatSlc] = runReplication( ...
            b, exretMat, predictorMat, TVec, TMin, K1, restrictMat, ...
            signRestriction, weightRestriction, cfg);
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
end


% ====================================================================
%  Local helpers — entry-point-specific.
% ====================================================================

function [econSlc, dmSlc, cwSlc, cwDiffSlc, tStatSlc] = runReplication( ...
    b, exretMat, predictorMat, TVec, TMin, K1, restrictMat, ...
    signRestriction, weightRestriction, cfg)
% One bootstrap replication: 3 model variants (different sample sizes
% TVec) for the same b'th draw. Returns slices sized to fit the parent
% econMat (3 x 3 x 2), dmMat / cwMat (3 x 3), cwDiffMat (1 x 3),
% tStatAlphaMat (3 x 2).

econSlc   = NaN(3, 3, 2);
dmSlc     = NaN(3, 3);
cwSlc     = NaN(3, 3);
cwDiffSlc = NaN(1, 3);
tStatSlc  = NaN(3, 2);

for ii = 1:3
    T   = TVec(ii);
    % Sticky-expectations sim uses the inflated 2.7y kernel bandwidth
    % (matches dailyAssetPricingBootstrap.m and the FST published
    % baseline), NOT cfg.coefWindowYears = 2.5.
    [h, out] = computeKernelBandwidth(2.7, T, 'daily');

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
    T   = T - 1;

    % Kernel forecast + recursive winsorization.
    [~, ~, ~, yForecast] = lps1_v2(y, [ones(T,1), X], h, 0, K1, 1.2, ...
        restrictMat(ii));
    yForecast = winsorizeRecursiveFast(yForecast, cfg.winsorPct, ...
        100 - cfg.winsorPct, out + 1);

    % G-prior shrinkage with the wider [weightLB, weightUB] clip range
    % stickyExp uses (vs the [0,1] range hard-coded in gPriorWeights).
    yCombinedG = applyGPriorWide(yForecast, y, yem, out, cfg);

    % Align all series to the OOS sample boundary.
    yF1     = yem(2*out:end);
    yF2     = yCombinedG(end-numel(yF1)+1:end);
    if signRestriction
        yF2(yF2 < 0) = 0;
    end
    yActual = y(2*out:end);

    % Smooth squared-error differential and detect pockets.
    fDM = (yActual - yF1).^2 - (yActual - yF2).^2;
    [~, ~, ~, ~, ~, fDMHat] = lps1_v2(fDM, ...
        [ones(size(fDM)), (1:numel(fDM))'], 252/numel(fDM), 0, K1, 1.2, 0);

    pocketIndices = stickyExpPocketFilter(fDM, fDMHat, TMin, ...
        cfg.minPocketDays);

    % Trim everything to the common sample.
    yF1     = yF1(end-TMin+1:end);
    yF2     = yF2(end-TMin+1:end);
    yActual = yActual(end-TMin+1:end);

    fCW = (yActual - yF1).^2 - ((yActual - yF2).^2 - (yF1 - yF2).^2);
    fDM = (yActual - yF1).^2 - (yActual - yF2).^2;

    dmSlc(ii, :)        = computeDMCWStats(fDM, pocketIndices, 'dm');
    cwFull              = computeDMCWStats(fCW, pocketIndices, 'cw');
    cwSlc(ii, :)        = cwFull(1:3);
    cwDiffSlc(ii)       = cwFull(4);

    % Portfolio metrics: 2 specs (TVC vs prevailing-mean baseline).
    [~, yPocketTime] = computeGomezCramWeights( ...
        [yF2, yF1], yActual, cfg.weightLB, cfg.weightUB, weightRestriction);

    [alphaA, srA, deltaA, tStatA] = computeUtilityPanel( ...
        yPocketTime, yActual, 3);
    econSlc(ii, :, :)   = [alphaA; srA; deltaA];
    tStatSlc(ii, :)     = tStatA;
end
end


function yCombinedG = applyGPriorWide(yForecast, y, yem, out, cfg)
% G-prior shrinkage loop. The G-prior weights are clipped to [0, 1]
% (matching subroutines/forecasting/gPriorWeights) — NOT to the wider
% portfolio-weight range (cfg.weightLB / cfg.weightUB), which applies
% only to the Gomez-Cram portfolio weights computed downstream.
yCombinedG    = yForecast;
weightWindow  = round(cfg.weightWindowYears * 252);
g             = cfg.gPrior;
bet0          = [1; 0];
benchmarkForecast = computeBenchmarkForecast(yForecast, yem, cfg);

for t = max(weightWindow, out) + 1 : (numel(yForecast) - 1)
    yReg   = y(t + out - weightWindow : t + out - 1);
    F      = yForecast(t - weightWindow + 1 : t);
    betHat = F \ yReg;
    betHat = min(max(betHat, 0), 1);

    w2 = bet0(2) + (1 / (1 + g)) * (betHat - bet0(2));
    w2 = min(max(w2, 0), 1);

    yCombinedG(t + 1) = (1 - w2) * benchmarkForecast(t + 1) + ...
                        w2 * yForecast(t + 1);
end
end


function pocketIndices = stickyExpPocketFilter(fDM, fDMHat, TMin, minLength)
% Custom min-length pocket filter (not detectPocketsFromSED). Clips the
% first `minLength` entries of each pocket; a final-period pocket is
% retained fully because the loop's cleanup branch does not fire when
% pocketIndices(TMin) is true.
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
% return (via per-period fsolve to preserve the original calibration),
% and the alpha t-stat.
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
