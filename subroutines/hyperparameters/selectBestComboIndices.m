function [bestSharpeInd, bestAlphaInd, bestTAlphaInd, bestRMSEInd] = ...
    selectBestComboIndices(sharpeMat, alphaMat, tAlphaMat, rmseMat, criterion)
% SELECTBESTCOMBOINDICES  Pick the best-combo index per metric. Centralizes
% the small repeating "max / min over combos" calls in Hyperparameters2
% and Hyperparameters_Marginals.
%
%   [bestSharpeInd, bestAlphaInd, bestTAlphaInd, bestRMSEInd] = ...
%       selectBestComboIndices(sharpeMat, alphaMat, tAlphaMat, rmseMat, criterion)
%
% Inputs:
%   sharpeMat / alphaMat / tAlphaMat / rmseMat  vectors or matrices of
%       per-combo statistics. NaN values are ignored.
%   criterion  unused at the moment; kept for API stability so callers
%       can pass cfg.selectionCriterion. The R1 hook for swapping the
%       published Adaptive panel selection is wired in
%       dailyEmpiricsHyperparametersMarginals.
%
% Outputs:
%   *Ind  scalar indices of the best combo by each metric.
%
% Sharpe / alpha / tAlpha favor the maximum; RMSE favors the minimum.
%
% NOTE: the published Adaptive (max alpha) panel selection is built
% in dailyEmpiricsHyperparametersMarginals / adaptivePanelsForSplit
% directly via max(alphMatExpanding); it does NOT route through this
% helper. bestSharpeInd is only kept here as a diagnostic scaffold for
% the R1 selectionCriterion='sharpe' hook (not currently exercised by
% any rendered table or figure).

if nargin < 5; criterion = ''; end %#ok<NASGU>

[~, bestSharpeInd] = max(sharpeMat, [], 'omitnan');
[~, bestAlphaInd]  = max(alphaMat,  [], 'omitnan');
[~, bestTAlphaInd] = max(tAlphaMat, [], 'omitnan');
[~, bestRMSEInd]   = min(rmseMat,   [], 'omitnan');
end
