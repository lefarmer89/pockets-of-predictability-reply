function [paramCombs, idx] = constructHyperparameterGrid()
% CONSTRUCTHYPERPARAMETERGRID  9,720-row grid of adaptive hyperparameter
% combinations consumed by dailyEmpiricsHyperparameters1 and the
% downstream Hyperparameters2 / Marginals analyses.
%
%   [paramCombs, idx] = constructHyperparameterGrid()
%
% Outputs:
%   paramCombs   numCombs x 8 matrix. Columns indexed by idx fields.
%                Total = 6 x 3 x 3 x 3 x 2 x 3 x 2 x 5 = 9720
%                (4860 zero-shrinkage + 4860 prevailing-mean-shrinkage).
%   idx          struct with column indices:
%                  .winsorization, .gPrior, .minWeight, .maxWeight,
%                  .benchmark, .windowLength, .timeTrend, .minPocket
%
% Marginals filters by idx.benchmark to either the bench=0 (zero-shrinkage)
% or bench=1 (prevailing-mean-shrinkage) half, leaving 4,860 combos per
% downstream artifact. The published reply uses the bench=0 subset.
%
% The combvec / fliplr ordering is preserved exactly so the prevInd
% cache in Hyperparameters1 hits the same matches it always has.

winsorizationVec = [0, 0.5, 1, 1.5, 2, 2.5];
gPriorVec        = [1, 2, 3];
minWeightVec     = [0, -0.5, -1];
maxWeightVec     = [1, 1.5, 2];
benchmarkVec     = [0, 1];
windowLengthVec  = [0.5, 1, 1.5];
timeTrendVec     = [0, 1];
minPocketVec     = [0, 11, 21, 32, 42];

paramCombs = fliplr(combvec( ...
    minPocketVec, timeTrendVec, windowLengthVec, ...
    benchmarkVec, maxWeightVec, minWeightVec, gPriorVec, ...
    winsorizationVec)');

idx = struct( ...
    'winsorization', 1, ...
    'gPrior',        2, ...
    'minWeight',     3, ...
    'maxWeight',     4, ...
    'benchmark',     5, ...
    'windowLength',  6, ...
    'timeTrend',     7, ...
    'minPocket',     8);
end
