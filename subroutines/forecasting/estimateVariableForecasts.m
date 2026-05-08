function R = estimateVariableForecasts(S, hSpec, restrictMatVal, signRestriction, benchmarkFlag, cfg)
% ESTIMATEVARIABLEFORECASTS  Per-variable kernel estimation, G-prior
% shrinkage, forecast-differential construction, and pocket detection.
%
% Wraps the inner variable-loop body shared by dailyEmpirics,
% monthlyEmpirics, and famaFrenchEmpirics.
%
%   R = estimateVariableForecasts(S, hSpec, restrictMatVal, signRestriction, benchmarkFlag, cfg)
%
% Inputs:
%   S                struct from prepareDataForEstimation. Required fields:
%                    X, y, riskFree, naiveForecast, T.
%   hSpec            struct with kernel parameters. Required fields:
%                      h               bandwidth (windowYears*obsPerYear/T)
%                      out             warm-up exclusion (floor(h*T))
%                      order           local-poly order (always 0 in callers)
%                      Ke1             kernel function handle
%                      dmLengthYears   SED smoothing window (years)
%                      obsPerYear      252 daily, 12 monthly
%                      minRT           min pocket length to qualify
%                      weightWindow    G-prior weight estimation window
%                      sedTimeTrend    optional logical (default false).
%                                      Set true for famaFrenchEmpirics
%                                      which uses k=2 SED design with
%                                      time trend.
%   restrictMatVal   scalar coef-restriction code passed to lps1_v2 for
%                    the slope coefficient.
%   signRestriction  logical; clip negative forecasts to zero when true.
%   benchmarkFlag    'pm' (prevailing mean) or 'lpm' (kernel) — selects
%                    the benchmark for the fCW/fDM differentials.
%   cfg              cfg struct (uses .winsorPct, .shrinkageTarget).
%
% Output struct R with fields:
%   yF1, yF1PM, yF2, yF2Alt   aligned forecasts (length T - 2*out + 1)
%   a, weightG                TVC slope (T x 1) and G-prior weights
%   gamDM                     smoothed SED coefficient (T-extend x 1)
%   r2Local, r2Local1y        local R^2 from TVC regression (T x 1)
%   r2LocalPM, r2Gam          PM local R^2 and SED smoothing R^2
%   pocketIndices             logical pocket indicator (T - out x 1)
%   pocketLengths             length of each detected pocket
%   pocketPeriods             [start, end] indices per pocket
%   integralR2Pocket          sum of r2Local within each pocket
%   fDMHat                    smoothed SED
%   yActual, riskFree         aligned realized series
%   fCW, fDM                  forecast differentials at the benchmark
%   yForecastWinsorized       full-length winsorized forecast (for callers
%                             that need pre-alignment access)
%   yForecastPM               full-length prevailing-mean forecast (ditto)

X             = S.X;
y             = S.y;
T             = S.T;
naiveForecast = S.naiveForecast;
riskFree      = S.riskFree;

h             = hSpec.h;
out           = hSpec.out;
order         = hSpec.order;
Ke1           = hSpec.Ke1;
dmLenYears    = hSpec.dmLengthYears;
obsPerYear    = hSpec.obsPerYear;
minRT         = hSpec.minRT;
weightWindow  = hSpec.weightWindow;
if isfield(hSpec, 'sedTimeTrend') && hSpec.sedTimeTrend
    sedTimeTrend = true;
else
    sedTimeTrend = false;
end

% --- Prevailing-mean kernel forecast (intercept-only design).
[~, ~, r2LocalPM, yForecastPM, ~, ~] = lps1_v2(y, ones(T,1), h, order, Ke1, [], 0);
if signRestriction
    yForecastPM(yForecastPM < 0) = 0;
end

% --- Time-varying coefficient kernel regression with the predictor.
[a, ~, r2Local, yForecast, ~, ~] = lps1_v2(y, [ones(T,1), X], h, ...
    order, Ke1, 1.2, restrictMatVal);

% Recursively winsorize the kernel forecast.
yForecast = winsorizeRecursiveFast(yForecast, cfg.winsorPct, 100 - cfg.winsorPct, out + 1);

a       = a(:, end);
r2Local = [0; r2Local(1:end-1)];

% 1-year local R^2 (used by Tables A.1 / A.2).
[~, ~, r2Local1y, ~, ~, ~] = lps1_v2(y, [ones(T,1), X], obsPerYear/T, ...
    order, Ke1, 1.2, restrictMatVal);
r2Local1y = [0; r2Local1y(1:end-1)];

% --- G-prior shrinkage of the kernel forecast toward the configured target.
benchmarkForecast = computeBenchmarkForecast(yForecast, naiveForecast, cfg);
[yCombinedG, weightG, ~] = gPriorWeightsFast(yForecast, y, naiveForecast, ...
    benchmarkForecast, weightWindow, out, cfg.gPrior, 0, 1);

% Align all series at the OOS sample boundary.
yF1PM       = naiveForecast(2*out:end);
yF1         = yForecastPM(out+1:end);
weightG     = weightG(out+1:end, :);
yCombinedG  = yCombinedG(out+1:end);
yF2         = yCombinedG;
if signRestriction
    yF2(yF2 < 0) = 0;
end
yActual  = y(2*out:end);
riskFree = riskFree(2*out:end);

% Forecast differentials against the configured benchmark.
if strcmp(benchmarkFlag, 'pm')
    [fCW, fDM] = forecastDifferentials(yActual, yF1PM, yF2);
else
    [fCW, fDM] = forecastDifferentials(yActual, yF1,   yF2);
end

% --- Smooth squared-error differential and extract pocket runs.
if sedTimeTrend
    sedDesign = [ones(size(fDM)), (1:numel(fDM))'];
else
    sedDesign = ones(size(fDM));
end
[gamDM, ~, r2Gam, ~, ~, fDMHat] = lps1_v2(fDM, sedDesign, ...
    dmLenYears * obsPerYear / numel(fDM), 0, Ke1, [], 0);
% Return the full gamDM matrix; callers slice as needed:
%   daily/monthly use gamDM(:,1) (intercept-only design)
%   famaFrench uses gamDM(:,2) (time-trend column from k=2 design)

[pocketIndices, pocketPeriods, pocketLengths] = ...
    detectPocketsFromSED(fDMHat, numel(r2Local), minRT);

integralR2Pocket = NaN(numel(pocketLengths), 1);
if ~isempty(pocketPeriods)
    for k = 1:size(pocketPeriods, 1)
        integralR2Pocket(k) = sum(r2Local(pocketPeriods(k,1):pocketPeriods(k,2)));
    end
    yF2Alt = yF2 .* pocketIndices + yF1PM .* (~pocketIndices);
else
    yF2Alt = yF1PM;
end

R = struct( ...
    'yF1',                  yF1, ...
    'yF1PM',                yF1PM, ...
    'yF2',                  yF2, ...
    'yF2Alt',               yF2Alt, ...
    'a',                    a, ...
    'weightG',              weightG, ...
    'gamDM',                gamDM, ...
    'r2Local',              r2Local, ...
    'r2Local1y',            r2Local1y, ...
    'r2LocalPM',            r2LocalPM, ...
    'r2Gam',                r2Gam, ...
    'pocketIndices',        pocketIndices, ...
    'pocketLengths',        pocketLengths, ...
    'pocketPeriods',        pocketPeriods, ...
    'integralR2Pocket',     integralR2Pocket, ...
    'fDMHat',               fDMHat, ...
    'yActual',              yActual, ...
    'riskFree',             riskFree, ...
    'fCW',                  fCW, ...
    'fDM',                  fDM, ...
    'yForecastWinsorized',  yForecast, ...
    'yForecastPM',          yForecastPM);
end
