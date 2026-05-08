function dailyEmpiricsGPriorRobustness(cfg)
% DAILYEMPIRICSGPRIORROBUSTNESS  Run the daily empirics at g=1 and g=3
% to populate the G-prior robustness panels of Table A.1 (tab:gRobustness
% in the reply). The published baseline (g=2) is produced by
% dailyEmpirics(cfg) directly.
%
%   dailyEmpiricsGPriorRobustness(cfg)
%
% Inputs (read from cfg.paths.data):
%   Daily_Predictors.xlsx, pcPredictor.mat
%
% Outputs:
%   cfg.paths.results/gPriorRobustness/g1/forecastResults_*.mat, OOSResults_*.mat
%   cfg.paths.results/gPriorRobustness/g3/forecastResults_*.mat, OOSResults_*.mat
%
% Tables/figures consumed by:
%   Table A.1 (G-prior robustness)
%
% Runtime: ~10-30 minutes (two sequential dailyEmpirics calls).

if nargin < 1 || isempty(cfg); cfg = default_config(); end

for g = [1, 3]
    cfgG = cfg;
    cfgG.gPrior        = g;
    cfgG.resultsSubdir = fullfile('gPriorRobustness', sprintf('g%d', g));
    fprintf('dailyEmpiricsGPriorRobustness: g = %d\n', g);
    dailyEmpirics(cfgG);
end
end
