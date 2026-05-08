function [yCombined, weightG, betHatVec] = gPriorWeights( ...
    yForecast, y, naiveForecast, benchmarkForecast, weightWindow, out, cfg)
% GPRIORWEIGHTS  Compute the G-prior shrinkage weights and the combined
% (shrunken) forecast at each time t, mirroring the empirical procedure
% described in Section III of the reply.
%
%   [yCombined, weightG, betHatVec] = gPriorWeights( ...
%       yForecast, y, naiveForecast, benchmarkForecast, weightWindow, out, cfg)
%
% Inputs:
%   yForecast          T x 1 one-sided kernel forecast.
%   y                  T x 1 dependent (excess returns).
%   naiveForecast      T x 1 prevailing-mean forecast (used to seed the
%                      combined series in the warm-up region).
%   benchmarkForecast  T x 1 shrinkage target (zeros by default; see
%                      computeBenchmarkForecast for the cfg.shrinkageTarget
%                      hook).
%   weightWindow       rolling-window length for weight estimation, in
%                      observations. NaN means expanding window.
%   out                burn-in observation count (= floor(h*T)).
%   cfg                config struct. Reads cfg.gPrior.
%
% Outputs:
%   yCombined  T x 1 shrunken forecast.
%   weightG    T x 2 weights on (benchmark, forecast) at each t.
%   betHatVec  T x 1 raw OLS shrinkage weight before applying g-prior.

g = cfg.gPrior;
T = numel(yForecast);
weightG = NaN(T, 2);
yCombined = naiveForecast(out:end);
betHatVec = NaN(T, 1);

bet0 = [1; 0];
minWeight = 0;
maxWeight = 1;

startT = max(weightWindow, out) + 1;
for t = startT:(T - 1)
    if isnan(weightWindow)
        % Expanding window
        yReg = y(out:t+out-1) - benchmarkForecast(1:t);
        F    = yForecast(1:t) - benchmarkForecast(1:t);
    else
        % Rolling window of length weightWindow
        yReg = y(t+out-weightWindow : t+out-1) - benchmarkForecast(t-weightWindow+1 : t);
        F    = yForecast(t-weightWindow+1 : t) - benchmarkForecast(t-weightWindow+1 : t);
    end

    betHat = F \ yReg;
    betHat = min(max(betHat, minWeight), maxWeight);
    betHatVec(t) = betHat;

    % G-prior shrinkage toward bet0(2) = 0
    w2 = bet0(2) + (1 / (1 + g)) * (betHat - bet0(2));
    w2 = min(max(w2, minWeight), maxWeight);
    weightG(t, 2) = w2;
    weightG(t, 1) = 1 - w2;

    yCombined(t + 1) = weightG(t, 1) * benchmarkForecast(t + 1) + ...
                       weightG(t, 2) * yForecast(t + 1);
end
end
