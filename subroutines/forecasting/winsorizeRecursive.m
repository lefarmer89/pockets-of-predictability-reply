function yWins = winsorizeRecursive(yForecast, qLB, qUB, startIdx)
% WINSORIZERECURSIVE  Recursively winsorize a forecast series in
% expanding-window fashion: at each t, clip yForecast(t) to the
% [qLB, qUB] empirical percentiles of yForecast(1:t).
%
%   yWins = winsorizeRecursive(yForecast, qLB, qUB, startIdx)
%
% Inputs:
%   yForecast  T x 1 forecast vector
%   qLB        lower percentile (e.g. 2.5)
%   qUB        upper percentile (e.g. 97.5)
%   startIdx   first index to winsorize. Earlier indices are passed
%              through unchanged. Equal to floor(h*T)+1 in callers.
%
% Output:
%   yWins      T x 1 winsorized forecast.

yWins = yForecast;
n = numel(yForecast);
for t = startIdx:n
    curLB = prctile(yForecast(1:t), qLB);
    curUB = prctile(yForecast(1:t), qUB);
    yWins(t) = max(min(yForecast(t), curUB), curLB);
end
end
