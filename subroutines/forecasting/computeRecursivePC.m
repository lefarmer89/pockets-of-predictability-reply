function predictors = computeRecursivePC(predictorMat)
% COMPUTERECURSIVEPC  Recursive expanding-window first principal component.
%
%   predictors = computeRecursivePC(predictorMat)
%
% At each t (starting from firstInd+nx where firstInd is the earliest
% all-non-NaN row of predictorMat and nx = size(predictorMat, 2)), fit
% PCA on rows firstInd:t and store the first PC. Returns a (T x K)
% matrix where K = T - firstInd - nx + 1 is the number of recursive PCs.
%
% Used by monthlyEmpirics (always) and dailyEmpirics for non-default
% cfg.sampleSplit (the precomputed pcPredictor.mat is only valid for
% the full-sample range; sample-split sub-samples must recompute).

firstInd = max(sum(isnan(predictorMat))) + 1;
nx       = size(predictorMat, 2);
N        = size(predictorMat, 1);
predictors = NaN(N, N - firstInd - nx + 1);
for t = firstInd + nx : N
    [~, scores] = pca(predictorMat(firstInd:t, :));
    predictors(firstInd:t, t - firstInd - nx + 1) = scores(:, 1);
end
end
