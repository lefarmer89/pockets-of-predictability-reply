function bf = computeBenchmarkForecast(yForecast, prevailingMean, cfg)
% COMPUTEBENCHMARKFORECAST  Pick the shrinkage target for the G-prior
% Bayesian shrinkage. Defaults to a vector of zeros (the reply's setting);
% an alternative shrinks toward the expanding-window prevailing mean.
%
%   bf = computeBenchmarkForecast(yForecast, prevailingMean, cfg)
%
% Inputs:
%   yForecast       T x 1 vector of one-sided kernel forecasts
%   prevailingMean  T x 1 expanding-window mean of the dependent variable
%                   (length must match yForecast). Pass [] if the caller
%                   only wants the zero target.
%   cfg             struct with at least field shrinkageTarget. If empty,
%                   defaults to 'zero'.
%
% Output:
%   bf              T x 1 benchmark forecast for the shrinkage update.
%
% Cfg hook: cfg.shrinkageTarget = 'prevailingMean' selects the
% prevailing-mean alternative.

if nargin < 3 || isempty(cfg) || ~isfield(cfg, 'shrinkageTarget')
    target = 'zero';
else
    target = cfg.shrinkageTarget;
end

switch target
    case 'zero'
        bf = zeros(size(yForecast));
    case 'prevailingMean'
        if isempty(prevailingMean)
            error('cfg.shrinkageTarget=''prevailingMean'' but prevailingMean argument is empty');
        end
        if numel(prevailingMean) ~= numel(yForecast)
            % If lengths differ, align to the tail of yForecast
            n = numel(yForecast);
            bf = prevailingMean(end-n+1:end);
            bf = bf(:);
        else
            bf = prevailingMean(:);
        end
    otherwise
        error('unknown cfg.shrinkageTarget: %s', target);
end
end
