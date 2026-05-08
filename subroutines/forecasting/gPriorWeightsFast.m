function [yCombined, weightG, betHatVec] = gPriorWeightsFast( ...
    yForecast, y, naiveForecast, benchmarkForecast, weightWindow, out, ...
    g, minW, maxW)
%GPRIORWEIGHTSFAST  Fast drop-in for gPriorWeights with O(1)-per-step
% rolling-sum updates of the F'F and F'yReg quantities used in the
% univariate OLS shrinkage weight.
%
%   [yCombined, weightG, betHatVec] = gPriorWeightsFast( ...
%       yForecast, y, naiveForecast, benchmarkForecast, weightWindow, out, ...
%       g, minW, maxW)
%
% Compared to gPriorWeights:
%   - g, minW, maxW are explicit args (not unpacked from cfg). cfg-using
%     callers pass `cfg.gPrior, 0, 1`. The combo-specific gPrior block
%     in dailyEmpiricsHyperparameters1.m's evaluateCombo passes its
%     paramCombs-derived g/minW/maxW directly.
%   - Per-step OLS computed via a sliding-window sum: maintain
%       sumFF(t) = sum_{i in window} F(i)^2
%       sumFy(t) = sum_{i in window} F(i) * yReg(i)
%       nanCount(t) = number of NaN positions in the window
%     and update by add-one/drop-one. Reduces per-call cost from
%     O(T*W) to O(T) (W = weightWindow, default 252).
%
% NaN handling: matches gPriorWeights — when ANY position in the
% current window has NaN in F or yReg, betHat is NaN. Tracked via a
% running nanCount; sums are kept finite by replacing NaN with 0 in
% the working arrays.

if nargin < 8 || isempty(minW); minW = 0; end
if nargin < 9 || isempty(maxW); maxW = 1; end

T = numel(yForecast);
weightG    = NaN(T, 2);
yCombined  = naiveForecast(out:end);
betHatVec  = NaN(T, 1);

% Element-wise products. yReg is aligned to the global index k (1..T):
%   yReg(k) = y(k+out-1) - benchmarkForecast(k)
% so the same array works for both expanding and rolling windows.
F_elem    = yForecast - benchmarkForecast;
yReg_elem = y(out:out+T-1) - benchmarkForecast(1:T);
FF_elem   = F_elem .* F_elem;
Fy_elem   = F_elem .* yReg_elem;

% NaN-safe versions for the running sums (replace NaN with 0; track count).
isNaNF  = isnan(FF_elem) | isnan(Fy_elem);
FF_safe = FF_elem;  FF_safe(isNaNF) = 0;
Fy_safe = Fy_elem;  Fy_safe(isNaNF) = 0;

startT = max(weightWindow, out) + 1;

if isnan(weightWindow)
    % Expanding window: cumulative sums.
    sumFF    = sum(FF_safe(1:startT-1));
    sumFy    = sum(Fy_safe(1:startT-1));
    nanCount = sum(isNaNF(1:startT-1));
    for t = startT:(T-1)
        sumFF    = sumFF    + FF_safe(t);
        sumFy    = sumFy    + Fy_safe(t);
        nanCount = nanCount + isNaNF(t);

        if nanCount > 0
            betHat = NaN;
        else
            betHat = sumFy / sumFF;
        end
        betHat = min(max(betHat, minW), maxW);
        betHatVec(t) = betHat;

        w2 = (1 / (1 + g)) * betHat;          % bet0(2) = 0
        w2 = min(max(w2, minW), maxW);
        weightG(t, 2) = w2;
        weightG(t, 1) = 1 - w2;
        yCombined(t + 1) = (1 - w2) * benchmarkForecast(t + 1) + ...
                            w2     * yForecast(t + 1);
    end
else
    % Rolling window of length W: sliding-window sums + NaN count.
    W = weightWindow;
    sumFF    = sum(FF_safe(startT-W+1 : startT-1));
    sumFy    = sum(Fy_safe(startT-W+1 : startT-1));
    nanCount = sum(isNaNF(startT-W+1 : startT-1));
    for t = startT:(T-1)
        % Add the new (incoming) element at index t.
        sumFF    = sumFF    + FF_safe(t);
        sumFy    = sumFy    + Fy_safe(t);
        nanCount = nanCount + isNaNF(t);

        if nanCount > 0
            betHat = NaN;
        else
            betHat = sumFy / sumFF;
        end
        betHat = min(max(betHat, minW), maxW);
        betHatVec(t) = betHat;

        w2 = (1 / (1 + g)) * betHat;
        w2 = min(max(w2, minW), maxW);
        weightG(t, 2) = w2;
        weightG(t, 1) = 1 - w2;
        yCombined(t + 1) = (1 - w2) * benchmarkForecast(t + 1) + ...
                            w2     * yForecast(t + 1);

        % Drop the outgoing element (window shifts by 1 next iteration).
        sumFF    = sumFF    - FF_safe(t - W + 1);
        sumFy    = sumFy    - Fy_safe(t - W + 1);
        nanCount = nanCount - isNaNF(t - W + 1);
    end
end
end
