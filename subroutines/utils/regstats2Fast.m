function s = regstats2Fast(y, data, model, whichstats)
%REGSTATS2FAST  Fast drop-in for regstats2 covering the call patterns the
% replication package actually uses.
%
%   s = regstats2Fast(y, x)                          % 2-arg fast path
%   s = regstats2Fast(y, data, model, whichstats)    % full
%
% 2-arg fast path: equivalent to
%   regstats2(y, x, 'linear', {'beta','hac','mse'})
% with x univariate. Specialized for p=2; bypasses argument parsing and
% whichstats dispatch entirely. Used by the H2 / Marginals expanding-window
% inner loop where the same call shape is hit ~12k times per (spec, combo).
%
% 4-arg general path supports:
%   model:      'linear' (intercept added) or 'onlydata' (use data as-is)
%   whichstats: any subset of {'beta','mse','r','yhat','rsquare',
%               'adjrsquare','covb','hac'}
%
% NaN handling matches regstats2 in both paths: rows where y or any
% column of the design matrix is NaN are removed before estimation. If
% fewer than p (number of design columns) finite rows remain, returns a
% NaN-filled struct (consumers' isfield checks then propagate NaN).
%
% Numerical equivalence target vs regstats2 with the same arguments:
%   beta / mse: 1e-12; hac.t / hac.se: 1e-10. The HAC math is line-for-line
%   equivalent to regstats2.m:517-535.

% ===== 2-arg specialized fast path =====
if nargin == 2
    y = y(:);
    x = data(:);

    isClean = ~isnan(y) & ~isnan(x);
    if ~all(isClean)
        y = y(isClean);
        x = x(isClean);
    end
    n = numel(y);

    if n < 3
        s.beta   = [NaN; NaN];
        s.mse    = NaN;
        s.hac.se = [NaN; NaN];
        s.hac.t  = [NaN; NaN];
        return
    end

    % Closed-form OLS on [1, x].
    sx  = sum(x);
    sxx = sum(x .* x);
    sy  = sum(y);
    sxy = sum(x .* y);
    d   = n*sxx - sx*sx;
    beta = [sxx*sy - sx*sxy; n*sxy - sx*sy] / d;

    % Residuals and MSE.
    e   = y - beta(1) - beta(2)*x;
    mse = sum(e .* e) / (n - 2);

    % HAC sandwich (vectorized, p=2). hhat_i = [e_i; e_i*x_i]; precompute
    % ex = e .* x once to reuse across lags.
    ex = e .* x;
    xuux11 = sum(e .* e);
    xuux12 = sum(e .* ex);
    xuux22 = sum(ex .* ex);
    L = ceil(0.75 * n^(1/3));
    for l = 1:L
        w  = 1 - l/(L+1);
        eL  = e(l+1:n);    eLag  = e(1:n-l);
        xL  = ex(l+1:n);   xLag  = ex(1:n-l);
        z11 = sum(eL .* eLag);
        z12 = sum(eL .* xLag);
        z21 = sum(xL .* eLag);
        z22 = sum(xL .* xLag);
        xuux11 = xuux11 + 2*w*z11;
        xuux12 = xuux12 + w*(z12 + z21);
        xuux22 = xuux22 + 2*w*z22;
    end

    xtxi = [sxx, -sx; -sx, n] / d;
    covb = xtxi * [xuux11, xuux12; xuux12, xuux22] * xtxi;

    s.beta   = beta;
    s.mse    = mse;
    s.hac.se = sqrt([covb(1,1); covb(2,2)]);
    s.hac.t  = beta ./ s.hac.se;
    return
end

% ===== 4-arg general path =====
if nargin == 3 || isempty(model)
    if nargin < 3; model = 'linear'; end
end
if nargin < 4 || isempty(whichstats)
    whichstats = {'beta','hac','mse'};
elseif ischar(whichstats)
    whichstats = {whichstats};
end

y = y(:);
if isrow(data) && numel(data) == numel(y)
    data = data(:);
end

switch model
    case 'linear'
        X = [ones(size(data,1), 1), data];
    case 'onlydata'
        X = data;
    otherwise
        error('regstats2Fast:badModel', ...
            'model must be ''linear'' or ''onlydata''');
end

wasnan = isnan(y) | any(isnan(X), 2);
if any(wasnan)
    y = y(~wasnan);
    X = X(~wasnan, :);
end

[n, p] = size(X);

% Pre-init requested fields to false (matches regstats2's degenerate-case
% behavior).
varnames = {'beta','mse','r','yhat','rsquare','adjrsquare','covb', ...
            'hac','empty','rankdef'};
hasStat  = false(1, numel(varnames));
for j = 1:numel(whichstats)
    sj = whichstats{j};
    ix = strcmp(sj, varnames);
    if any(ix)
        hasStat = hasStat | ix;
    elseif strcmp(sj, 'all')
        hasStat(:) = true;
    end
end
s = struct();
for j = find(hasStat)
    s.(varnames{j}) = false;
end

if n == 0
    if hasStat(strcmp('empty', varnames))
        s.empty = true;
    end
    return
elseif n < p
    if hasStat(strcmp('rankdef', varnames))
        s.rankdef = true;
    end
    return
end

% ----- OLS -----
if p == 1
    sxx  = sum(X .* X);
    sxy  = sum(X .* y);
    beta = sxy / sxx;
    xtxi = 1 / sxx;
elseif p == 2
    x2  = X(:, 2);
    sx  = sum(x2);
    sxx = sum(x2 .* x2);
    sy  = sum(y);
    sxy = sum(x2 .* y);
    d   = n*sxx - sx*sx;
    beta = [sxx*sy - sx*sxy; n*sxy - sx*sy] / d;
    xtxi = [sxx, -sx; -sx, n] / d;
else
    XtX  = X' * X;
    Xty  = X' * y;
    R    = chol(XtX);
    beta = R \ (R' \ Xty);
    Rinv = R \ eye(p);
    xtxi = Rinv * Rinv';
end

yhat = X * beta;
e    = y - yhat;
sse  = sum(e .* e);
dfe  = n - p;
mse  = sse / dfe;

if hasStat(1); s.beta       = beta; end
if hasStat(2); s.mse        = mse;  end
if hasStat(3); s.r          = e;    end
if hasStat(4); s.yhat       = yhat; end
if hasStat(5)
    sst = sum((y - mean(y)).^2);
    s.rsquare = 1 - sse / sst;
end
if hasStat(6)
    sst = sum((y - mean(y)).^2);
    dft = n - 1;
    s.adjrsquare = 1 - (sse/sst) * (dft/dfe);
end
if hasStat(7)
    s.covb = xtxi * mse;
end

% ----- HAC -----
if hasStat(8)
    L = ceil(0.75 * n^(1/3));

    if p == 1
        ex = e .* X;
        xuux = sum(ex .* ex);
        for l = 1:L
            w  = 1 - l/(L+1);
            za = sum(ex(l+1:n) .* ex(1:n-l));
            xuux = xuux + 2*w*za;
        end
    elseif p == 2
        x2 = X(:, 2);
        ex = e .* x2;
        xuux11 = sum(e .* e);
        xuux12 = sum(e .* ex);
        xuux22 = sum(ex .* ex);
        for l = 1:L
            w  = 1 - l/(L+1);
            eL  = e(l+1:n);    eLag  = e(1:n-l);
            xL  = ex(l+1:n);   xLag  = ex(1:n-l);
            z11 = sum(eL .* eLag);
            z12 = sum(eL .* xLag);
            z21 = sum(xL .* eLag);
            z22 = sum(xL .* xLag);
            xuux11 = xuux11 + 2*w*z11;
            xuux12 = xuux12 + w*(z12 + z21);
            xuux22 = xuux22 + 2*w*z22;
        end
        xuux = [xuux11, xuux12; xuux12, xuux22];
    else
        hhat = X' .* e';
        xuux = hhat * hhat';
        for l = 1:L
            w  = 1 - l/(L+1);
            za = hhat(:, l+1:n) * hhat(:, 1:n-l)';
            xuux = xuux + w*(za + za');
        end
    end

    covbHAC = xtxi * xuux * xtxi;
    if p == 1
        seHAC = sqrt(covbHAC);
    else
        seHAC = sqrt(diag(covbHAC));
    end

    h.beta = beta;
    h.se   = seHAC;
    h.t    = beta ./ seHAC;
    h.pval = 2 * tcdf(-abs(h.t), dfe);
    h.dfe  = dfe;
    h.covb = covbHAC;
    s.hac  = h;
end
end
