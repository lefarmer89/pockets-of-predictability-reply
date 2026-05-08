function yAR1 = computeAR1Forecast(dependent, startT)
% COMPUTEAR1FORECAST  Recursive expanding-window AR(1) forecast.
%
%   yAR1 = computeAR1Forecast(dependent, startT)
%
% Inputs:
%   dependent  T x 1 series to forecast.
%   startT     index from which forecasts are computed (entries below this
%              are zero). Daily callers use 252.
%
% Output:
%   yAR1  T x 1 vector with yAR1(t) = expanding-window AR(1) forecast of
%         dependent(t), regressing dependent(2:t-1) on
%         [1, dependent(1:t-2)]. Entries 1:(startT-1) are zero.

yAR1 = zeros(size(dependent));
for t = startT:size(dependent, 1)
    bet = [ones(t-2, 1), dependent(1:t-2)] \ dependent(2:t-1);
    yAR1(t) = [1, dependent(t-1)] * bet;
end
end
