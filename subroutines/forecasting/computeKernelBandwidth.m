function [h, out] = computeKernelBandwidth(windowLengthYears, T, freq)
% COMPUTEKERNELBANDWIDTH  One-sided kernel bandwidth and warm-up exclusion.
%
%   [h, out] = computeKernelBandwidth(windowLengthYears, T, freq)
%
% Inputs:
%   windowLengthYears  rolling-window length in years (e.g. 2.5)
%   T                  number of observations to fit (post-staggering for
%                      empirics callers, raw T-1 for the pre-loop call in
%                      empirics that sizes the result matrices)
%   freq               'daily' (252 obs/year) or 'monthly' (12 obs/year)
%
% Outputs:
%   h    Bandwidth as a fraction of T: windowLengthYears * obsPerYear / T
%   out  Warm-up exclusion: floor(h * T)
%
% The famaFrench daily variant uses an inflated factor (5.4 / 2.7 in place
% of 5 / 2.5) to compensate for the Epanechnikov kernel's effective
% sample-size loss. Pass that pre-multiplied windowLengthYears (e.g. 2.7
% instead of 2.5) when calling from famaFrenchEmpirics.

switch freq
    case 'daily'
        obsPerYear = 252;
    case 'monthly'
        obsPerYear = 12;
    otherwise
        error('computeKernelBandwidth:badFreq', ...
            'freq must be ''daily'' or ''monthly''');
end

h   = windowLengthYears * obsPerYear / T;
out = floor(h * T);
end
