function S = prepareDataForEstimation(varName, predictors, dependent, riskFree, dates, predictorMat)
% PREPAREDATAFORESTIMATION  Trim leading NaNs and stagger X / y for the
% one-step-ahead prediction setup used by every empirics entry point.
%
%   S = prepareDataForEstimation(varName, predictors, dependent, ...
%                                riskFree, dates, predictorMat)
%
% Inputs:
%   varName       'dp','tbl','tsp','rvar','mv','pc','erL'
%   predictors    T x k raw predictor matrix
%   dependent     T x 1 dependent variable (e.g. data.exret)
%   riskFree      T x 1 risk-free rate (e.g. data.rf)
%   dates         T x 1 dates
%   predictorMat  T x 4 matrix used by the 'pc' branch to size the trim;
%                 ignored for other variables but kept in the signature
%                 so callers always pass the same fields.
%
% Output struct S with fields:
%   X              (Tnew x kx) staggered predictor matrix (NaNs zeroed)
%   y              (Tnew x 1)  staggered dependent
%   dates          (Tnew x 1)  staggered dates
%   riskFree       (Tnew x 1)  staggered risk-free
%   naiveForecast  (Tnew x 1)  expanding prevailing-mean forecast
%   yem            (Tnew x 1)  expanding-mean shifted to align with y
%   T              scalar      Tnew, after staggering
%   trim           scalar      number of leading observations dropped

switch varName
    case 'pc'
        trim = sum(isnan(predictors(:, end))) + size(predictorMat, 2);
    otherwise
        trim = max(sum(isnan(predictors)));
end

dates      = dates(trim+1:end);
yem        = cumsum(dependent) ./ (1:size(dependent, 1))';
predictors = predictors(trim+1:end, :);
dependent  = dependent(trim+1:end);
riskFree   = riskFree(trim+1:end);
yem        = yem(trim+1:end);
[T, ~]     = size(predictors);

% Stagger so X(t) predicts y(t+1).
X     = predictors;
y     = dependent;
dates = dates(1:end-1);
switch varName
    case 'pc'
        X = X(1:end-1, 1:end-1);
    otherwise
        X = X(1:end-1, :);
end
y        = y(2:end);
riskFree = riskFree(2:end);

naiveForecast = yem(1:end-1);
yem           = yem(2:end);
T             = T - 1;

X(isnan(X)) = 0;

S = struct( ...
    'X',             X, ...
    'y',             y, ...
    'dates',         dates, ...
    'riskFree',      riskFree, ...
    'naiveForecast', naiveForecast, ...
    'yem',           yem, ...
    'T',             T, ...
    'trim',          trim);
end
