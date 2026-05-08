function row = computeEconomicsMetrics(yPocketTime, yActual, freq)
% COMPUTEECONOMICSMETRICS  Per-replication portfolio metrics for bootstrap
% or simulation: alpha (annualized), CAPM-style Sharpe (annualized), and
% the per-step alpha t-stat. The Delta certainty-equivalent return
% is caller-specific (depends on the assumed utility) and is computed
% inline by the caller.
%
%   row = computeEconomicsMetrics(yPocketTime, yActual, freq)
%
% Inputs:
%   yPocketTime  T x 1 portfolio return time series (single replication)
%   yActual      T x 1 realized returns
%   freq         'daily' (252) or 'monthly' (12)
%
% Output:
%   row  1 x 3: [alphaAnnual, SRAnnual, tStatAlpha]

switch freq
    case {'daily', 'famaFrench'}; annFactor = 252;
    case 'monthly';               annFactor = 12;
    otherwise
        error('computeEconomicsMetrics:badFreq', ...
            'freq must be ''daily'', ''monthly'', or ''famaFrench''');
end

s = regstats2Fast(yPocketTime, yActual);
alphaAnnual = s.beta(1) * annFactor * 100;
tStatAlpha  = s.hac.t(1);

SRmkt    = mean(yActual,    'omitnan') / std(yActual,    'omitnan');
SR       = sqrt(SRmkt^2 + (s.beta(1)^2) / s.mse);
SRAnnual = SR * sqrt(annFactor);

row = [alphaAnnual, SRAnnual, tStatAlpha];
end
