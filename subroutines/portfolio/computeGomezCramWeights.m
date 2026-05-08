function [portfolioWeight, yPocketTime] = computeGomezCramWeights( ...
    yForecasts, yActual, weightLB, weightUB, weightRestriction)
% COMPUTEGOMEZCRAMWEIGHTS  Bootstrap-flavored Gomez-Cram portfolio weight
% construction (no vol/mom factor side-effects).
%
%   [portfolioWeight, yPocketTime] = computeGomezCramWeights( ...
%       yForecasts, yActual, weightLB, weightUB, weightRestriction)
%
% Inputs:
%   yForecasts         T x S matrix of one-step-ahead forecasts.
%   yActual            T x 1 realized excess returns.
%   weightLB           lower clip (typically 0)
%   weightUB           upper clip (typically 2)
%   weightRestriction  logical: if true, clip weights to [weightLB, weightUB]
%
% Outputs:
%   portfolioWeight  T x S after optional clipping
%   yPocketTime      portfolioWeight .* yActual
%
% Used by dailyBootstrap, dailyAssetPricingBootstrap, and stickyExpectationsSim
% where the leaner per-replication form is enough (no vol/mom factor needed).

cVec = sqrt(var(yActual, 'omitnan') ./ ...
            var(yForecasts .* yActual, 'omitnan'));
portfolioWeight = cVec .* yForecasts;

if weightRestriction
    portfolioWeight(portfolioWeight < weightLB) = weightLB;
    portfolioWeight(portfolioWeight > weightUB) = weightUB;
end

yPocketTime = portfolioWeight .* yActual;
end
