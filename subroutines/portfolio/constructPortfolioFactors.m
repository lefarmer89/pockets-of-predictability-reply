function F = constructPortfolioFactors(yForecastMat, yActual, riskFree, freq, weightRestriction, adjCostBps)
% CONSTRUCTPORTFOLIOFACTORS  Gomez-Cram portfolio construction, vol-signed
% factor, 1y momentum factor, and transaction-cost adjustment term.
%
%   F = constructPortfolioFactors(yForecastMat, yActual, riskFree, freq, weightRestriction, adjCostBps)
%
% Inputs:
%   yForecastMat       T x numSpecifications matrix of OOS forecasts to
%                      convert to portfolio weights.
%   yActual            T x 1 realized excess returns.
%   riskFree           T x 1 risk-free rate (matched units).
%   freq               'daily' or 'monthly'. Sets the moving-window sizes
%                      and the form of the transaction-cost adj term.
%   weightRestriction  logical; if true, clip weights to [0, 2].
%   adjCostBps         per-trade transaction cost in basis points
%                      (default 0 — published baseline). Pass 5 or 10
%                      to populate the Table A.3 transaction-cost columns.
%
% Output struct F with fields:
%   portfolioWeight  T x numSpecifications, post-clip if weightRestriction
%   yPocketTime      portfolioWeight .* yActual minus adjTerm*adjCost
%   ySigFactor       T x 1 vol-signed factor
%   yMomFactor       T x 1 1y momentum-signed factor
%   adjTerm          T x numSpecifications transaction-cost adjustment
%
% Window sizes (short / long): daily 21 / 252, monthly 3 / 12.
% Adj formula:
%   daily       → geometric self-financing (preserves dailyEmpirics behaviour)
%   monthly     → additive return-difference (matches monthlyEmpirics)
%   famaFrench  → daily windows but the additive form, matching the
%                 original famaFrenchEmpirics convention.

switch freq
    case 'daily'
        wShort = 21;  wLong = 252;  useGeoAdj = true;
    case 'monthly'
        wShort = 3;   wLong = 12;   useGeoAdj = false;
    case 'famaFrench'
        wShort = 21;  wLong = 252;  useGeoAdj = false;
    otherwise
        error('constructPortfolioFactors:badFreq', ...
            'freq must be ''daily'', ''monthly'', or ''famaFrench''');
end

% Gomez-Cram c-scaling so signal volatility matches yActual volatility.
cVec = sqrt(var(yActual, 'omitnan') ./ ...
            var(yForecastMat .* yActual, 'omitnan'));
portfolioWeight = cVec .* yForecastMat;

rvar1M = movsum((yActual - movmean(yActual, [wShort, 0])).^2, [wShort, 0]);

portfolioWeightUnadj = portfolioWeight;
if weightRestriction
    portfolioWeight(portfolioWeight < 0) = 0;
    portfolioWeight(portfolioWeight > 2) = 2;
end
yPocketTime = portfolioWeight .* yActual;

% Vol-signed factor (daily/monthly use the same definition with
% frequency-specific window sizes).
ySigFactor = yActual ./ [NaN(wShort+1, 1); rvar1M(wShort+1:end-1)];
ySigFactor = ySigFactor ./ std(ySigFactor, 'omitnan') .* std(yActual, 'omitnan');

% 1-year momentum factor.
cumRet1Y   = movsum(yActual, [wLong-1, 0]);
yMomFactor = [NaN(wLong, 1); sign(cumRet1Y(wLong:end-1))] .* yActual ./ ...
             [NaN(wShort+1, 1); sqrt(rvar1M(wShort+1:end-1))];
yMomFactor = yMomFactor ./ std(yMomFactor, 'omitnan') .* std(yActual, 'omitnan');

numSpecifications = size(yPocketTime, 2);

if useGeoAdj
    % Daily: self-financing geometric weight evolution.
    adjTerm = [NaN(1, numSpecifications); ...
               abs(portfolioWeight(2:end,:) - ...
                   min(max(portfolioWeight(1:end-1,:) .* exp(yActual(1:end-1)) ./ ...
                       (portfolioWeight(1:end-1,:) .* exp(yActual(1:end-1)) + ...
                        (1 - portfolioWeight(1:end-1,:)) .* exp(riskFree(1:end-1))), 0), 2))];
else
    % Monthly: additive return-difference variant on unclipped weights.
    adjTerm = [NaN(1, numSpecifications); ...
               abs(portfolioWeight(2:end,:) - ...
                   min(max(portfolioWeightUnadj(1:end-1,:) .* ...
                       exp(diff(yActual + riskFree - yPocketTime)), 0), 2))];
end

if nargin < 6 || isempty(adjCostBps); adjCostBps = 0; end
adjCost = adjCostBps / 100 / 100;
yPocketTime = yPocketTime - adjTerm * adjCost;

F = struct( ...
    'portfolioWeight', portfolioWeight, ...
    'yPocketTime',     yPocketTime, ...
    'ySigFactor',      ySigFactor, ...
    'yMomFactor',      yMomFactor, ...
    'adjTerm',         adjTerm);
end
