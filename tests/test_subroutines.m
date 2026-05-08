function r = test_subroutines()
% TEST_SUBROUTINES  Exercise shared helper functions on synthetic data.
%
% Verifies that the helper API hasn't drifted: lps1_v2, regstats2,
% opt_block_length, LatexTableFull each accept expected inputs and return
% expected shapes. Synthetic data only; no .mat or .xlsx is loaded.

r.name = 'test_subroutines';
r.pass = false;
r.message = '';

rng(42);
T = 500;
y = randn(T,1);
X = [ones(T,1), randn(T,1)];

checks = {};

% ---- regstats2 ----
try
    s = regstats2(y, X(:,2), 'linear', {'beta','hac','adjrsquare'});
    assert(numel(s.beta) == 2, 'beta length');
    assert(isfield(s.hac, 't'), 'hac.t missing');
    assert(isfield(s.hac, 'covb'), 'hac.covb missing');
    assert(numel(s.hac.t) == 2, 'hac.t length');
    assert(isscalar(s.adjrsquare), 'adjrsquare scalar');
    checks{end+1} = 'regstats2 OK';
catch ME
    r.message = sprintf('regstats2 failed: %s', ME.message);
    return
end

% ---- regstats2Fast equivalence ----
% Fast must match regstats2 to 1e-12 on beta/mse and 1e-10 on hac.t/hac.se
% across these regimes:
%   (a) 'linear' p=2 (univariate predictor + intercept) — H2 hot path
%   (b) 'linear' p=3 (two predictors + intercept) — perfMetrics 2-factor
%   (c) 'linear' p=4 (three predictors + intercept) — perfMetrics 3-factor
%   (d) 'onlydata' p=1 (intercept-only design) — DM/CW HAC of mean
%   (e) NaN-stripping in (a)
%   (f) {'beta','hac','adjrsquare','yhat'} for dailyBootstrap-style call
try
    % (a) p=2 / linear / {beta,hac,mse}
    sFast = regstats2Fast(y, X(:,2));
    sRef  = regstats2(y, X(:,2), 'linear', {'beta','hac','mse'});
    dB2 = max(abs(sFast.beta - sRef.beta));
    dM2 = abs(sFast.mse - sRef.mse);
    dT2 = max(abs(sFast.hac.t  - sRef.hac.t));
    dS2 = max(abs(sFast.hac.se - sRef.hac.se));
    assert(dB2 < 1e-12 && dM2 < 1e-12 && dT2 < 1e-10 && dS2 < 1e-10, ...
        sprintf('p=2 drift dB=%.3g dM=%.3g dT=%.3g dS=%.3g', dB2, dM2, dT2, dS2));

    % (b) p=3 / linear / {beta,hac,mse}
    rng(7);
    X3data = randn(T, 2);
    sFast = regstats2Fast(y, X3data, 'linear', {'beta','hac','mse'});
    sRef  = regstats2(y, X3data, 'linear', {'beta','hac','mse'});
    dB3 = max(abs(sFast.beta - sRef.beta));
    dM3 = abs(sFast.mse - sRef.mse);
    dT3 = max(abs(sFast.hac.t - sRef.hac.t));
    assert(dB3 < 1e-12 && dM3 < 1e-12 && dT3 < 1e-10, ...
        sprintf('p=3 drift dB=%.3g dM=%.3g dT=%.3g', dB3, dM3, dT3));

    % (c) p=4 / linear / {beta,hac,mse}
    X4data = randn(T, 3);
    sFast = regstats2Fast(y, X4data, 'linear', {'beta','hac','mse'});
    sRef  = regstats2(y, X4data, 'linear', {'beta','hac','mse'});
    dB4 = max(abs(sFast.beta - sRef.beta));
    dT4 = max(abs(sFast.hac.t - sRef.hac.t));
    assert(dB4 < 1e-12 && dT4 < 1e-10, ...
        sprintf('p=4 drift dB=%.3g dT=%.3g', dB4, dT4));

    % (d) 'onlydata' p=1 (intercept-only) — DM/CW pattern
    sFast = regstats2Fast(y, ones(T,1), 'onlydata', 'hac');
    sRef  = regstats2(y, ones(T,1), 'onlydata', 'hac');
    dT1 = max(abs(sFast.hac.t - sRef.hac.t));
    dS1 = max(abs(sFast.hac.se - sRef.hac.se));
    assert(dT1 < 1e-10 && dS1 < 1e-10, ...
        sprintf('onlydata p=1 drift dT=%.3g dS=%.3g', dT1, dS1));

    % (e) NaN-stripping equivalence in p=2.
    yN = y; xN = X(:,2);
    yN([5, 100, 250]) = NaN;
    xN([20, 200])     = NaN;
    sFastN = regstats2Fast(yN, xN);
    sRefN  = regstats2(yN, xN, 'linear', {'beta','hac','mse'});
    dBN = max(abs(sFastN.beta - sRefN.beta));
    dTN = max(abs(sFastN.hac.t - sRefN.hac.t));
    assert(dBN < 1e-12 && dTN < 1e-10, ...
        sprintf('NaN-strip drift dB=%.3g dT=%.3g', dBN, dTN));

    % (f) {beta,hac,adjrsquare,yhat} multivariate (dailyBootstrap pattern)
    sFast = regstats2Fast(y, X4data, 'linear', {'beta','hac','adjrsquare','yhat'});
    sRef  = regstats2(y, X4data, 'linear', {'beta','hac','adjrsquare','yhat'});
    dAR = abs(sFast.adjrsquare - sRef.adjrsquare);
    dY  = max(abs(sFast.yhat - sRef.yhat));
    assert(dAR < 1e-12 && dY < 1e-12, ...
        sprintf('adjrsquare/yhat drift dAR=%.3g dY=%.3g', dAR, dY));

    checks{end+1} = sprintf(['regstats2Fast OK ' ...
        '(p=2 dB=%.1e dT=%.1e; p=3 dB=%.1e dT=%.1e; p=4 dB=%.1e dT=%.1e; ' ...
        'onlydata p=1 dT=%.1e; NaN-strip dB=%.1e dT=%.1e; ' ...
        'adjR2/yhat dAR=%.1e dY=%.1e)'], ...
        dB2, dT2, dB3, dT3, dB4, dT4, dT1, dBN, dTN, dAR, dY);
catch ME
    r.message = sprintf('regstats2Fast failed: %s', ME.message);
    return
end

% ---- opt_block_length ----
try
    bl = opt_block_length(y);
    assert(all(bl > 0), 'block length positive');
    assert(numel(bl) >= 1, 'block length non-empty');
    checks{end+1} = sprintf('opt_block_length OK (block=%.1f)', bl(1));
catch ME
    r.message = sprintf('opt_block_length failed: %s', ME.message);
    return
end

% ---- stationary_bootstrap basic sanity ----
% Confirms the restored Politis-Romano helper returns the expected shapes
% and that bsdata == data(indices_modwrapped) (the post-doubling lookup).
try
    rng(2024);
    Tsb = 500;
    [bs, ind] = stationary_bootstrap((1:Tsb)', 5, 20);
    assert(isequal(size(bs),  [Tsb, 5]), 'bsdata shape');
    assert(isequal(size(ind), [Tsb, 5]), 'indices shape');
    indMod = mod(ind - 1, Tsb) + 1;
    assert(all(bs(:) == indMod(:)), 'bsdata == data(mod(indices))');
    checks{end+1} = sprintf('stationary_bootstrap OK (shape=[%d %d], wrap-consistent)', ...
        size(bs, 1), size(bs, 2));
catch ME
    r.message = sprintf('stationary_bootstrap failed: %s', ME.message);
    return
end

% ---- winsorizeRecursiveFast equivalence ----
% Must match winsorizeRecursive at machine precision, including the
% production case where yForecast has NaN in the warm-up region.
try
    rng(123);
    Tw = 2000;
    yF = randn(Tw, 1);
    refOut  = winsorizeRecursive(yF,     2.5, 97.5, 50);
    fastOut = winsorizeRecursiveFast(yF, 2.5, 97.5, 50);
    dW = max(abs(refOut - fastOut));
    assert(dW < 1e-12, sprintf('winsorizeRecursive drift %.3g >= 1e-12', dW));

    refOut2  = winsorizeRecursive(yF,     5, 95, 252);
    fastOut2 = winsorizeRecursiveFast(yF, 5, 95, 252);
    dW2 = max(abs(refOut2 - fastOut2));
    assert(dW2 < 1e-12, sprintf('winsorizeRecursive 5/95 drift %.3g', dW2));

    % NaN warm-up case (matches the lps1_v2 output pattern: first 252 NaN).
    yFnan = yF; yFnan(1:252) = NaN;
    refOutN  = winsorizeRecursive(yFnan,     2.5, 97.5, 253);
    fastOutN = winsorizeRecursiveFast(yFnan, 2.5, 97.5, 253);
    diffN = abs(refOutN - fastOutN);
    diffN(isnan(refOutN) & isnan(fastOutN)) = 0;
    dWnan = max(diffN);
    assert(dWnan < 1e-12, sprintf('NaN-warmup drift %.3g', dWnan));

    checks{end+1} = sprintf('winsorizeRecursiveFast OK (dW=%.1e dW2=%.1e dWnan=%.1e)', ...
        dW, dW2, dWnan);
catch ME
    r.message = sprintf('winsorizeRecursiveFast failed: %s', ME.message);
    return
end

% ---- gPriorWeightsFast equivalence ----
% Must match gPriorWeights at near-machine precision for both rolling
% (W = 252) and expanding (W = NaN) windows.
try
    rng(456);
    Tg = 1500;
    out = 200;
    yF      = randn(Tg, 1) * 0.01;
    yDep    = randn(Tg + out - 1, 1) * 0.01;
    naiveF  = yDep(1:Tg);
    benchF  = yF * 0;
    cfgG    = struct('gPrior', 5);

    % Rolling W = 252.
    [yRef,   wRef,   bRef]   = gPriorWeights(yF, yDep, naiveF, benchF, 252, out, cfgG);
    [yFast,  wFast,  bFast]  = gPriorWeightsFast(yF, yDep, naiveF, benchF, 252, out, cfgG.gPrior, 0, 1);
    dY = max(abs(yRef - yFast), [], 'omitmissing');
    dW = max(abs(wRef(:) - wFast(:)), [], 'omitmissing');
    dB = max(abs(bRef - bFast), [], 'omitmissing');
    assert(dY < 1e-10, sprintf('gPrior rolling y drift %.3g', dY));
    assert(dW < 1e-10, sprintf('gPrior rolling weight drift %.3g', dW));
    assert(dB < 1e-10, sprintf('gPrior rolling beta drift %.3g', dB));

    % Expanding (NaN window).
    [yRefE,  wRefE,  bRefE]  = gPriorWeights(yF, yDep, naiveF, benchF, NaN, out, cfgG);
    [yFastE, wFastE, bFastE] = gPriorWeightsFast(yF, yDep, naiveF, benchF, NaN, out, cfgG.gPrior, 0, 1);
    dYe = max(abs(yRefE - yFastE), [], 'omitmissing');
    dWe = max(abs(wRefE(:) - wFastE(:)), [], 'omitmissing');
    dBe = max(abs(bRefE - bFastE), [], 'omitmissing');
    assert(dYe < 1e-10, sprintf('gPrior expanding y drift %.3g', dYe));
    assert(dWe < 1e-10, sprintf('gPrior expanding weight drift %.3g', dWe));
    assert(dBe < 1e-10, sprintf('gPrior expanding beta drift %.3g', dBe));

    % NaN warm-up in yForecast (production pattern: first `out` NaN).
    yFn = yF; yFn(1:out) = NaN;
    benchFn = benchF; benchFn(1:out) = NaN;
    [yRefN,  ~, bRefN]  = gPriorWeights(yFn, yDep, naiveF, benchFn, 252, out, cfgG);
    [yFastN, ~, bFastN] = gPriorWeightsFast(yFn, yDep, naiveF, benchFn, 252, out, cfgG.gPrior, 0, 1);
    diffYn = abs(yRefN - yFastN);
    diffBn = abs(bRefN - bFastN);
    diffYn(isnan(yRefN) & isnan(yFastN)) = 0;
    diffBn(isnan(bRefN) & isnan(bFastN)) = 0;
    dYn = max(diffYn, [], 'omitmissing');
    dBn = max(diffBn, [], 'omitmissing');
    assert(dYn < 1e-10, sprintf('gPrior NaN-warmup y drift %.3g', dYn));
    assert(dBn < 1e-10, sprintf('gPrior NaN-warmup beta drift %.3g', dBn));

    checks{end+1} = sprintf(['gPriorWeightsFast OK ' ...
        '(rolling dY=%.1e dW=%.1e dB=%.1e; expanding dY=%.1e dW=%.1e dB=%.1e; ' ...
        'NaN-warmup dY=%.1e dB=%.1e)'], ...
        dY, dW, dB, dYe, dWe, dBe, dYn, dBn);
catch ME
    r.message = sprintf('gPriorWeightsFast failed: %s', ME.message);
    return
end

% ---- lps1_v2 ----
% Time-varying coefficient kernel regression on synthetic data, using the
% one-sided Epanechnikov kernel that all production callers use.
try
    Tlps = 300;
    Xlps = [ones(Tlps,1), randn(Tlps,1)];
    ylps = randn(Tlps,1);
    h = 0.3;
    p = 0;
    K = @(z) 1.5*(1-z.^2).*(abs(z) <= 1).*(z <= 0);
    k2 = 1.2;
    [a, sigDefault, r2, yForecast, ~, ~] = lps1_v2(ylps, Xlps, h, p, K, k2, 0);
    assert(size(a,2) == 2, 'a width');
    assert(~isempty(r2), 'r2 non-empty');
    assert(~isempty(yForecast), 'forecast non-empty');
    assert(isempty(sigDefault), 'computeHAC default false: sigmaHatA empty');
    % HAC opt-in: pass computeHAC = true.
    [~, sigOpt, ~, ~, ~, ~] = lps1_v2(ylps, Xlps, h, p, K, k2, 0, true);
    assert(isequal(size(sigOpt), [2, 2, size(a,1)]), 'computeHAC=true: sigma shape');
    % k > 2 path (used by the 'mv' multivariate case in the empirics).
    Xlps5 = [ones(Tlps,1), randn(Tlps, 4)];
    [a5, sig5, r25, yF5, ~, ~] = lps1_v2(ylps, Xlps5, h, p, K, k2, 0);
    assert(size(a5,2) == 5, 'k>2 a width');
    assert(isempty(sig5), 'k>2 sigma empty by default');
    assert(numel(r25) == numel(r2), 'k>2 r2 length matches k=2');
    assert(numel(yF5) == size(a5,1), 'k>2 yForecast length');
    % fullDesign path (k == T+1, used by the 'pc' case in dailyEmpirics).
    Xfd = [ones(Tlps,1), randn(Tlps, Tlps)];
    [aFD, sigFD, r2FD, yFFD, ~, ~] = lps1_v2(ylps, Xfd, h, p, K, k2, 0);
    assert(size(aFD,2) == 2, 'fullDesign a width (effective k=2)');
    assert(isempty(sigFD), 'fullDesign sigma empty by default');
    assert(numel(yFFD) == size(aFD,1), 'fullDesign yForecast length');
    assert(~isempty(r2FD), 'fullDesign r2 non-empty');
    checks{end+1} = sprintf('lps1_v2 OK (a=%dx%d, k>2 + fullDesign + HAC opt-in)', ...
        size(a,1), size(a,2));
catch ME
    r.message = sprintf('lps1_v2 failed: %s', ME.message);
    return
end

% ---- LatexTableFull ----
try
    A = [1.23 4.56; 7.89 0.12];
    pval = 0.01*ones(2);
    rowLabels = {'a','b'};
    colLabels = {'Col','x','y'};
    if exist('LatexTableFull','file') ~= 2
        r.message = 'LatexTableFull not on path';
        return
    end
    LatexTableFull(A, colLabels, rowLabels, '9.2f', pval, 0);
    checks{end+1} = 'LatexTableFull OK';
catch ME
    r.message = sprintf('LatexTableFull failed: %s', ME.message);
    return
end

% ---- computeKernelBandwidth ----
try
    [h, outBuf] = computeKernelBandwidth(2.5, 5040, 'daily');
    assert(abs(h - 2.5*252/5040) < 1e-15, 'h mismatch');
    assert(outBuf == floor(h*5040), 'out mismatch');
    [hM, ~] = computeKernelBandwidth(2.5, 600, 'monthly');
    assert(abs(hM - 2.5*12/600) < 1e-15, 'monthly h mismatch');
    checks{end+1} = sprintf('computeKernelBandwidth OK (h=%.4g)', h);
catch ME
    r.message = sprintf('computeKernelBandwidth failed: %s', ME.message);
    return
end

% ---- computeAR1Forecast ----
try
    rng(7);
    dep = randn(300,1);
    yA = computeAR1Forecast(dep, 252);
    assert(numel(yA) == 300, 'yAR1 length mismatch');
    assert(all(yA(1:251) == 0), 'pre-startT entries should be zero');
    % Spot-check: yA(252) should equal AR(1) on dep(1:251) at lag dep(251).
    bet = [ones(250,1), dep(1:250)] \ dep(2:251);
    expected = [1, dep(251)] * bet;
    assert(abs(yA(252) - expected) < 1e-12, 'AR(1) value mismatch');
    checks{end+1} = 'computeAR1Forecast OK';
catch ME
    r.message = sprintf('computeAR1Forecast failed: %s', ME.message);
    return
end

% ---- prepareDataForEstimation ----
try
    rng(11);
    Traw = 200;
    pred = [NaN(20,1); randn(Traw-20,1)];   % first 20 NaN
    dep  = randn(Traw,1);
    rf   = 0.001*ones(Traw,1);
    dts  = (1:Traw)';
    pmat = randn(Traw,4);
    Sout = prepareDataForEstimation('dp', pred, dep, rf, dts, pmat);
    assert(Sout.trim == 20, 'trim mismatch');
    assert(Sout.T == Traw - 20 - 1, 'T mismatch');
    assert(numel(Sout.y) == Sout.T, 'y length mismatch');
    assert(numel(Sout.naiveForecast) == Sout.T, 'naiveForecast length');
    checks{end+1} = sprintf('prepareDataForEstimation OK (T=%d)', Sout.T);
catch ME
    r.message = sprintf('prepareDataForEstimation failed: %s', ME.message);
    return
end

% ---- computeForecastTests ----
try
    rng(31);
    n = 3;
    fDMCell = cell(n,1); fCWCell = cell(n,1); pIndCell = cell(n,1);
    for i = 1:n
        fDMCell{i} = randn(500,1);
        fCWCell{i} = randn(500,1);
        pIndCell{i} = rand(500,1) > 0.7;
    end
    [dmM, cwM, cwD] = computeForecastTests(fDMCell, fCWCell, pIndCell);
    assert(isequal(size(dmM),  [n 3]), 'dmMat size');
    assert(isequal(size(cwM),  [n 3]), 'cwMat size');
    assert(isequal(size(cwD),  [n 1]), 'cwDiff size');
    % Empty-pocket case must give NaN, not a fall-back stat.
    pIndCell2 = pIndCell; pIndCell2{1} = false(500,1);
    [dmM2, cwM2, ~] = computeForecastTests(fDMCell, fCWCell, pIndCell2);
    assert(isnan(dmM2(1,2)) && isnan(cwM2(1,2)), 'empty in-pocket NaN');
    checks{end+1} = 'computeForecastTests OK';
catch ME
    r.message = sprintf('computeForecastTests failed: %s', ME.message);
    return
end

% ---- computeDMCWStats ----
try
    rng(41);
    diffSer = randn(400,1);
    pIdx    = rand(400,1) > 0.6;
    rowDM = computeDMCWStats(diffSer, pIdx, 'dm');
    rowCW = computeDMCWStats(diffSer, pIdx, 'cw');
    assert(numel(rowDM) == 3, 'dm row width');
    assert(numel(rowCW) == 4, 'cw row width');
    checks{end+1} = 'computeDMCWStats OK';
catch ME
    r.message = sprintf('computeDMCWStats failed: %s', ME.message);
    return
end

% ---- constructPortfolioFactors ----
try
    rng(53);
    Tp = 1000;
    yAct  = randn(Tp,1)*0.01;
    rfp   = 0.0001*ones(Tp,1);
    yFmat = randn(Tp,3)*0.01;
    Fout = constructPortfolioFactors(yFmat, yAct, rfp, 'daily', true);
    assert(all(Fout.portfolioWeight(:) >= 0 & Fout.portfolioWeight(:) <= 2), ...
        'weights not clipped to [0,2]');
    assert(isequal(size(Fout.yPocketTime), [Tp 3]), 'yPocketTime size');
    checks{end+1} = 'constructPortfolioFactors OK';
catch ME
    r.message = sprintf('constructPortfolioFactors failed: %s', ME.message);
    return
end

% ---- computePerformanceMetrics ----
try
    rng(61);
    Tp = 800;
    yAct = randn(Tp,1)*0.01;
    yPT  = repmat(yAct, 1, 4) + randn(Tp,4)*0.005;
    yS   = randn(Tp,1)*0.01;
    yM   = randn(Tp,1)*0.01;
    Pout = computePerformanceMetrics(yPT, yAct, yS, yM, 'daily');
    assert(isequal(size(Pout.coefMat), [2 4]), 'coefMat size');
    assert(numel(Pout.alphaAnnual) == 4, 'alphaAnnual length');
    assert(Pout.SRmktAnnual == Pout.SRmkt * sqrt(252), 'SR annualization');
    checks{end+1} = 'computePerformanceMetrics OK';
catch ME
    r.message = sprintf('computePerformanceMetrics failed: %s', ME.message);
    return
end

% ---- computeGomezCramWeights ----
try
    rng(67);
    Tp = 500;
    yAct = randn(Tp,1)*0.01;
    yF   = randn(Tp,2)*0.05;
    [pw, ypt] = computeGomezCramWeights(yF, yAct, 0, 2, true);
    assert(all(pw(:) >= 0 & pw(:) <= 2), 'clipping');
    assert(isequal(size(ypt), [Tp 2]), 'ypt size');
    [pw0, ~] = computeGomezCramWeights(yF, yAct, 0, 2, false);
    assert(any(pw0(:) < 0) || any(pw0(:) > 2), 'unclipped should have outliers');
    checks{end+1} = 'computeGomezCramWeights OK';
catch ME
    r.message = sprintf('computeGomezCramWeights failed: %s', ME.message);
    return
end

% ---- computeEconomicsMetrics ----
try
    rng(71);
    Tp = 600;
    yAct = randn(Tp,1)*0.01;
    yPT  = yAct + randn(Tp,1)*0.005;
    rowE = computeEconomicsMetrics(yPT, yAct, 'daily');
    assert(numel(rowE) == 3, 'econ row width');
    assert(all(isfinite(rowE)), 'econ row finite');
    checks{end+1} = 'computeEconomicsMetrics OK';
catch ME
    r.message = sprintf('computeEconomicsMetrics failed: %s', ME.message);
    return
end

% ---- trimAndAlignResults ----
try
    Rin = struct();
    Rin.yF2Mat   = [NaN(5,2); randn(95,2)];
    Rin.yActual  = randn(100,1);
    Rin.yF1Mat   = [NaN(5,2); randn(95,2)];
    Rin.yF1PMMat = [NaN(5,2); randn(95,2)];
    Rin.riskFree = randn(100,1);
    Rin.dateVec  = (1:100)';
    Rin.pocketIndMat = [NaN(5,2); rand(95,2) > 0.5];
    Rout = trimAndAlignResults(Rin, 100);
    assert(size(Rout.yF2Mat,1) == 95, 'yF2Mat trim');
    assert(size(Rout.yActual,1) == 95, 'yActual trim');
    assert(all(~isnan(Rout.pocketIndMat(:))), 'pocketIndMat NaN -> 0');
    checks{end+1} = 'trimAndAlignResults OK';
catch ME
    r.message = sprintf('trimAndAlignResults failed: %s', ME.message);
    return
end

% ---- extractAndReshapeData ----
try
    Traw = 50; Bp = 4;
    raw = randn(Traw*Bp, 5);
    Sout2 = extractAndReshapeData(raw, Traw, Bp);
    assert(isequal(size(Sout2.exretMat), [Traw Bp]), 'exret size');
    assert(isequal(size(Sout2.predictorMat), [Traw Bp 3]), 'predictor size');
    checks{end+1} = 'extractAndReshapeData OK';
catch ME
    r.message = sprintf('extractAndReshapeData failed: %s', ME.message);
    return
end

r.pass = true;
r.message = strjoin(checks, '; ');
end
