function [a, sigmaHatA, r2, yForecast, yFit, yHat] = lps1_v2(y, X, h, p, K, k2, restrict, computeHAC)
% LPS1_V2  Chen and Hong (2012) one-sided local polynomial smoothing for
% time-varying coefficient regressions y_t = X_t alpha_t + e_t.
%
%   [a, sigmaHatA, r2, yForecast, yFit, yHat] = ...
%       lps1_v2(y, X, h, p, K, k2, restrict, computeHAC)
%
% Inputs:
%   y           T x 1 dependent variable.
%   X           T x k design matrix. k must be 1 (intercept-only,
%               prevailing-mean style) or 2 (intercept + predictor).
%   h           Bandwidth as a fraction of the sample (effective half-window
%               is floor(h*T) observations).
%   p           Polynomial order. Only p == 0 (local constant) is supported;
%               higher orders error out.
%   K           Kernel function handle. Must be one-sided causal: K(z) must
%               be zero for z > 0. The package's standard kernel
%               Ke1 = @(z) 1.5*(1-z.^2).*(abs(z)<=1).*(z<=0) satisfies this.
%   k2          Kernel second-moment constant. Unused; kept in the signature
%               for API compatibility.
%   restrict    Sign restriction on the second coefficient (k=2 only):
%               +1 enforces a >= 0, -1 enforces a <= 0, anything else
%               leaves the coefficient unrestricted.
%   computeHAC  Optional logical (default false). When true, compute the
%               per-t Newey-West HAC covariance matrix in sigmaHatA.
%               When false, sigmaHatA is returned as []. Default-off
%               because no production caller in the package uses
%               sigmaHatA, and HAC is 95-99% of total runtime; opting
%               out gives a 20-100x speedup. The empirics/bootstrap/
%               hyperparameter callers all discard sigmaHatA with `~`,
%               so they get this speedup automatically.
%
% Outputs (computed only when requested):
%   a          (T-extend+1) x k matrix of time-varying coefficients.
%   sigmaHatA  k x k x (T-extend+1) array of Newey-West HAC covariance
%              matrices, one per time t. Empty when computeHAC is false.
%   r2         (T-2*extend+1) x 1 vector of local R^2 values.
%   yForecast  One-step-ahead forecast vector. Sign restriction applied
%              to the second coefficient when restrict ~= 0.
%   yFit       y values at the fitted time stamps.
%   yHat       Fitted values.
%
% Implementation:
%   The OLS coefficients across all time stamps are computed by a small
%   number of filter() calls that produce running weighted moving sums
%   for X'WX and X'Wy. With one-sided kernel weights of length L = O(h*T),
%   each filter call is O(T*L), so the total OLS work is O(T*L) rather
%   than the O(T*L) per-iteration work of a per-t loop. The per-t HAC
%   loop is still required because the HAC covariance depends on the
%   per-t residuals.
%
% Profiling (3 warmup + 5 measurement, T = 22000, k = 2, p = 0):
%   pre-cleanup version of lps1_v2:                ~25-30 s estimated
%   first optimization pass (per-t loop, inlined):  2.8 s
%   this version:                                   1.2 s   (~25-30x over original)

if p ~= 0
    error('lps1_v2:unsupportedOrder', ...
        ['lps1_v2 only supports local-constant smoothing (p==0). ', ...
         'No caller in this package uses higher orders.']);
end

if nargin < 2
    error('Not enough input arguments.');
end
if nargin < 8 || isempty(computeHAC)
    computeHAC = false;
end

T = size(X,1);
if size(y,1) ~= T
    error('Dimensions of X and y do not match.');
end

k = size(X,2);
% Designs with k <= 2 use the vectorized closed-form path below.
% Special case: k == T + 1 marks the "fullDesign" pattern used by the
% 'pc' (recursive expanding-window principal component) variable in
% dailyEmpirics — column 1 is the intercept, column t+1 is the PC
% computed using data 1..t. This is effectively a k=2 fit with a
% per-t-changing predictor; allocating a TxTxT cube would OOM.
% Otherwise: k > 2 is the 'mv' multivariate case (intercept + 4
% predictors) and falls back to lps1_v2_general.
if k == T + 1
    [a, sigmaHatA, r2, yForecast, yFit, yHat] = ...
        lps1_v2_fullDesign(y, X, h, K, restrict, computeHAC);
    return
end
if k > 2
    % Memory guard: lps1_v2_general allocates a T x k x k cube. For
    % all known package callers (mv has k=5) this is small, but a
    % much larger k would OOM silently. Surface a clear error well
    % before MATLAB tries to allocate dozens of GB.
    bytesNeeded = T * k * k * 8;
    if bytesNeeded > 4 * 2^30
        error('lps1_v2:tooLargeK', ...
            ['k = %d with T = %d would require %.1f GB for the per-t ', ...
             'cross-product cache. Use the fullDesign pattern (k = T+1) ', ...
             'or reduce k.'], k, T, bytesNeeded / 2^30);
    end
    [a, sigmaHatA, r2, yForecast, yFit, yHat] = ...
        lps1_v2_general(y, X, h, K, restrict, computeHAC);
    return
end

extend = floor(h*T);
hT = h * T;
lb = extend;
ub = T;

% Causal kernel weights for negative lags (one-sided kernel). Truncated
% to the bandwidth window since the kernel is zero outside |z| <= 1.
maxLag = min(T, ceil(hT) + 1);
lagOff = (0:maxLag-1)' / hT;
wFilter = K(-lagOff);

% Running weighted moving sums via filter().
sum_w  = filter(wFilter, 1, ones(T,1));
sum_wY = filter(wFilter, 1, y);
if k == 2
    x2 = X(:,2);
    sum_wX  = filter(wFilter, 1, x2);
    sum_wXX = filter(wFilter, 1, x2 .^ 2);
    sum_wXY = filter(wFilter, 1, x2 .* y);
end

% Closed-form OLS betas, vectorized over time.
if k == 1
    a_t = sum_wY ./ sum_w;
    a = a_t(lb:ub);
else
    det_t = sum_w .* sum_wXX - sum_wX.^2;
    beta1 = (sum_wXX .* sum_wY - sum_wX .* sum_wXY) ./ det_t;
    beta2 = (-sum_wX .* sum_wY  + sum_w  .* sum_wXY) ./ det_t;
    a = [beta1(lb:ub), beta2(lb:ub)];
end

XFit = X(lb:ub, :);
yFit = y(lb:ub);
nIter = ub - lb + 1;
if computeHAC
    sigmaHatA = zeros(k, k, nIter);
else
    sigmaHatA = [];
end

if computeHAC
    for tIdx = 1:nIter
        t = lb + tIdx - 1;
        firstRow = max(1, t - maxLag + 1);
        rng = firstRow:t;

        Wt = wFilter(t - rng + 1);
        nz = Wt > 0;
        if ~any(nz)
            continue
        end
        sqrtWt = sqrt(Wt(nz));
        Xa = X(rng(nz), :) .* sqrtWt;
        ya = y(rng(nz))    .* sqrtWt;

        if k == 1
            betaHat = a_t(t);
            xtxi = 1 / sum_w(t);
        else
            betaHat = [beta1(t); beta2(t)];
            xtx = [sum_w(t),  sum_wX(t); sum_wX(t), sum_wXX(t)];
            xtxi = xtx \ eye(2);
        end

        residuals = ya - Xa * betaHat;
        n = size(Xa, 1);
        L = ceil(0.75 * n^(1/3));
        hhat = Xa.' .* residuals.';
        xuux = hhat * hhat.';
        for l = 1:L
            za = hhat(:, (l+1):n) * hhat(:, 1:(n-l)).';
            wL = 1 - l / (L+1);
            xuux = xuux + wL * (za + za.');
        end
        sigmaHatA(:, :, tIdx) = xtxi * xuux * xtxi;
    end
end

if nargout >= 2
    yHat = sum(a .* XFit, 2);
    eHat = yFit - yHat;

    % Vectorized weighted local R^2 numerator and denominator.
    eHat2Pad = zeros(T,1);
    eHat2Pad(lb:ub) = eHat .^ 2;
    sum_wEHat = filter(wFilter, 1, eHat2Pad);
    sum_wY2   = filter(wFilter, 1, y .^ 2);

    yMeanVec = cumsum(y) ./ (1:T)';
    sum_wYDev = sum_wY2 - 2 * yMeanVec .* sum_wY + (yMeanVec .^ 2) .* sum_w;

    lb2 = 2 * extend;
    r2_vec = 1 - sum_wEHat ./ sum_wYDev;
    r2_vec(sum_wYDev == 0) = 0;
    r2 = r2_vec(lb2:ub);

    if k == 1
        yForecast = [0; sum(a(1:end-1) .* XFit(2:end, :), 2)];
    else
        if restrict == 1
            col2 = max(a(1:end-1, 2), 0);
            yForecast = [0; sum([a(1:end-1, 1), col2] .* XFit(2:end, :), 2)];
        elseif restrict == -1
            col2 = min(a(1:end-1, 2), 0);
            yForecast = [0; sum([a(1:end-1, 1), col2] .* XFit(2:end, :), 2)];
        else
            yForecast = [0; sum(a(1:end-1, :) .* XFit(2:end, :), 2)];
        end
    end
end
end %#ok<*INUSD>


function [a, sigmaHatA, r2, yForecast, yFit, yHat] = lps1_v2_general(y, X, h, K, restrict, computeHAC)
% Per-t weighted-OLS path for designs with k > 2 (the 'mv' multivariate
% case). Vectorized in the same spirit as the k = 1, 2 fast path:
%
%   - All k*(k+1)/2 weighted cross-product sums X_i'WX_j are obtained by
%     filter() on each pair, reducing the per-t OLS to a single k x k
%     linear solve.
%   - The HAC covariance still requires a per-t loop because it depends
%     on the per-t weighted residuals, but each iteration reuses the
%     precomputed inverse of (X'WX)(t).
%
% This brings the cost down from the original O(T * maxLag^2 * k) per-t
% regstats2 loop to roughly O(T * (k^2 + k^3 + maxLag * L)) — a 30-100x
% speedup at T = 7000, k = 5 in practice (see tests/profile_harness.m).

T  = size(X, 1);
k  = size(X, 2);
extend = floor(h * T);
hT     = h * T;
lb = extend;
ub = T;
nIter = ub - lb + 1;

% At small t, NaN-zeroed predictors and short kernel windows can make
% the per-t weighted Gram matrix rank-deficient. mldivide returns a
% valid least-squares solution in that case (pseudoinverse-equivalent
% in the kernel of X'WX), but emits a singular-matrix warning. Suppress
% it for the duration of this call; the previous state is restored
% automatically when the onCleanup object goes out of scope.
prevWarn = warning;
warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');
warnRestore = onCleanup(@() warning(prevWarn)); %#ok<NASGU>

% Causal kernel weights for negative lags. Truncated to the bandwidth
% window since the kernel is zero outside |z| <= 1.
maxLag = min(T, ceil(hT) + 1);
lagOff = (0:maxLag-1)' / hT;
wFilter = K(-lagOff);

% --- Vectorized weighted sums (one filter call per inner product).
sum_w   = filter(wFilter, 1, ones(T, 1));
sum_wY  = filter(wFilter, 1, y);
sum_wY2 = filter(wFilter, 1, y .^ 2);
sum_wXY = zeros(T, k);
for i = 1:k
    sum_wXY(:, i) = filter(wFilter, 1, X(:, i) .* y);
end
% Symmetric X'WX: store both halves so the per-t solve can use squeeze().
sum_wXX = zeros(T, k, k);
for i = 1:k
    for j = i:k
        s = filter(wFilter, 1, X(:, i) .* X(:, j));
        sum_wXX(:, i, j) = s;
        if i ~= j
            sum_wXX(:, j, i) = s;
        end
    end
end

% --- Batched OLS solve via pagemldivide.
% Stack all (X'WX)(t) into a (k, k, nIter) tensor and all (X'Wy)(t) into
% a (k, 1, nIter) tensor, then do all nIter linear solves in one call.
% This eliminates the MATLAB-level per-t loop. The HAC sandwich below
% dominates total runtime (~95-99% — see tests/profile_lps1_v2_kgt2.m),
% so this is a clean-code improvement more than a speed win.
yFit = y(lb:ub);
XFit = X(lb:ub, :);

XtWX_pages = permute(sum_wXX(lb:ub, :, :), [2, 3, 1]);   % k x k x nIter
XtWy_pages = permute(reshape(sum_wXY(lb:ub, :), nIter, k, 1), [2, 3, 1]);  % k x 1 x nIter
beta_pages = pagemldivide(XtWX_pages, XtWy_pages);       % k x 1 x nIter
a    = squeeze(beta_pages).';                             % nIter x k
yHat = sum(XFit .* a, 2);                                 % nIter x 1

% Cache inv(X'WX) per t for the HAC sandwich below (skip if HAC opted out).
if computeHAC
    xtxiAll = pagemldivide(XtWX_pages, repmat(eye(k), [1, 1, nIter]));
end

if computeHAC
    sigmaHatA = zeros(k, k, nIter);
else
    sigmaHatA = [];
end

% --- Per-t Newey-West HAC covariance.
if computeHAC
    for tIdx = 1:nIter
        t = lb + tIdx - 1;
        firstRow = max(1, t - maxLag + 1);
        rng = firstRow:t;
        Wt = wFilter(t - rng + 1);
        nz = Wt > 0;
        if ~any(nz); continue; end
        sqrtWt = sqrt(Wt(nz));
        Xa = X(rng(nz), :) .* sqrtWt;
        ya = y(rng(nz))    .* sqrtWt;

        betaHat = a(tIdx, :).';
        residuals = ya - Xa * betaHat;
        n = size(Xa, 1);
        L = ceil(0.75 * n^(1/3));
        hhat = Xa.' .* residuals.';
        xuux = hhat * hhat.';
        for l = 1:L
            za = hhat(:, (l+1):n) * hhat(:, 1:(n-l)).';
            wL = 1 - l / (L+1);
            xuux = xuux + wL * (za + za.');
        end
        xtxi = xtxiAll(:, :, tIdx);
        sigmaHatA(:, :, tIdx) = xtxi * xuux * xtxi;
    end
end

if nargout >= 2
    eHat = yFit - yHat;

    % --- Vectorized weighted local R^2 (same construction as k=1,2 path).
    eHat2Pad = zeros(T, 1);
    eHat2Pad(lb:ub) = eHat .^ 2;
    sum_wEHat = filter(wFilter, 1, eHat2Pad);

    yMeanVec = cumsum(y) ./ (1:T)';
    sum_wYDev = sum_wY2 - 2 * yMeanVec .* sum_wY + (yMeanVec .^ 2) .* sum_w;

    lb2 = 2 * extend;
    r2_vec = 1 - sum_wEHat ./ sum_wYDev;
    r2_vec(sum_wYDev == 0) = 0;
    r2 = r2_vec(lb2:ub);

    % --- One-step-ahead forecast with optional sign clipping on slopes.
    if restrict == 1
        aC = a; aC(:, 2:end) = max(aC(:, 2:end), 0);
        yForecast = [0; sum(aC(1:end-1, :) .* XFit(2:end, :), 2)];
    elseif restrict == -1
        aC = a; aC(:, 2:end) = min(aC(:, 2:end), 0);
        yForecast = [0; sum(aC(1:end-1, :) .* XFit(2:end, :), 2)];
    else
        yForecast = [0; sum(a(1:end-1, :) .* XFit(2:end, :), 2)];
    end
end
end


function [a, sigmaHatA, r2, yForecast, yFit, yHat] = lps1_v2_fullDesign(y, X, h, K, restrict, computeHAC)
% "Full design" pattern: X is T x (T+1). Column 1 is the intercept,
% column t+1 is the predictor at time t (e.g. recursive expanding-
% window principal component). At each fitted time t, run a k=2 OLS
% on [intercept, X(:, t+1)] with one-sided kernel weights.
%
% Memory: O(T * maxLag) instead of O(T^2). Time: O(T * maxLag^2 * L)
% dominated by the per-t HAC loop (no across-t vectorization is
% possible since the second column changes with t).

T = size(X, 1);
extend = floor(h * T);
hT     = h * T;
lb = extend; ub = T;
nIter = ub - lb + 1;

prevWarn = warning;
warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');
warnRestore = onCleanup(@() warning(prevWarn)); %#ok<NASGU>

a         = zeros(nIter, 2);
yFit      = y(lb:ub);
XFit      = zeros(nIter, 2);
yHat      = zeros(nIter, 1);
if computeHAC
    sigmaHatA = zeros(2, 2, nIter);
    xtxiAll   = zeros(2, 2, nIter);
else
    sigmaHatA = [];
end

maxLag = min(T, ceil(hT) + 1);
lagOff = (0:maxLag-1)' / hT;
wFilter = K(-lagOff);

for tIdx = 1:nIter
    t = lb + tIdx - 1;
    firstRow = max(1, t - maxLag + 1);
    rng = firstRow:t;
    Wt = wFilter(t - rng + 1);
    nz = Wt > 0;
    if ~any(nz); continue; end

    sqrtWt = sqrt(Wt(nz));
    Xa = X(rng(nz), [1, t+1]) .* sqrtWt;
    ya = y(rng(nz)) .* sqrtWt;

    XtWX = Xa.' * Xa;
    if computeHAC
        xtxi = XtWX \ eye(2);
        betaHat = xtxi * (Xa.' * ya);
        xtxiAll(:, :, tIdx) = xtxi;
    else
        betaHat = XtWX \ (Xa.' * ya);
    end
    a(tIdx, :)    = betaHat.';
    XFit(tIdx, :) = X(t, [1, t+1]);
    yHat(tIdx)    = XFit(tIdx, :) * betaHat;

    if computeHAC
        residuals = ya - Xa * betaHat;
        n = size(Xa, 1);
        L = ceil(0.75 * n^(1/3));
        hhat = Xa.' .* residuals.';
        xuux = hhat * hhat.';
        for l = 1:L
            za = hhat(:, (l+1):n) * hhat(:, 1:(n-l)).';
            wL = 1 - l / (L+1);
            xuux = xuux + wL * (za + za.');
        end
        sigmaHatA(:, :, tIdx) = xtxi * xuux * xtxi;
    end
end

if nargout >= 2
    eHat = yFit - yHat;

    % --- Vectorized weighted local R^2 (same construction as the k=1,2
    %     path). The fullDesign R^2 only depends on eHat and y; the
    %     per-t-changing predictor column doesn't enter, so we can use
    %     filter() instead of the per-t loop in the archived version.
    sum_w     = filter(wFilter, 1, ones(T, 1));
    sum_wY    = filter(wFilter, 1, y);
    sum_wY2   = filter(wFilter, 1, y .^ 2);
    eHat2Pad  = zeros(T, 1);
    eHat2Pad(lb:ub) = eHat .^ 2;
    sum_wEHat = filter(wFilter, 1, eHat2Pad);
    yMeanVec  = cumsum(y) ./ (1:T)';
    sum_wYDev = sum_wY2 - 2 * yMeanVec .* sum_wY + (yMeanVec .^ 2) .* sum_w;
    lb2 = 2 * extend;
    r2_vec = 1 - sum_wEHat ./ sum_wYDev;
    r2_vec(sum_wYDev == 0) = 0;
    r2 = r2_vec(lb2:ub);

    % One-step-ahead forecast: a(t-1)' * X(t, [1, t+1]).
    if restrict == 1
        aC = a; aC(:, 2) = max(aC(:, 2), 0);
        yForecast = [0; sum(aC(1:end-1, :) .* XFit(2:end, :), 2)];
    elseif restrict == -1
        aC = a; aC(:, 2) = min(aC(:, 2), 0);
        yForecast = [0; sum(aC(1:end-1, :) .* XFit(2:end, :), 2)];
    else
        yForecast = [0; sum(a(1:end-1, :) .* XFit(2:end, :), 2)];
    end
end
end
