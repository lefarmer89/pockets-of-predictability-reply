function folder = baselineResultsFolder(paths, cfg, kind)
% BASELINERESULTSFOLDER  Resolve the per-cfg subfolder under
% paths.results where forecast / OOS results land for an entry-point
% function.
%
%   folder = baselineResultsFolder(paths, cfg)
%   folder = baselineResultsFolder(paths, cfg, kind)
%
% Routing:
%   When cfg.resultsSubdir is non-empty (set by, e.g.,
%   dailyEmpiricsGPriorRobustness), folder = paths.results/<cfg.resultsSubdir>.
%
%   Otherwise, folder = paths.results/<kind>, where kind is one of:
%       'daily'       (default; dailyEmpirics)
%       'monthly'     (monthlyEmpirics)
%       'famaFrench'  (famaFrenchEmpirics)
%
% The folder is created if it does not exist.

if nargin < 3 || isempty(kind); kind = 'daily'; end

if isstruct(cfg) && isfield(cfg, 'resultsSubdir') && ~isempty(cfg.resultsSubdir)
    folder = fullfile(paths.results, cfg.resultsSubdir);
else
    folder = fullfile(paths.results, kind);
end

if ~isfolder(folder)
    mkdir(folder);
end
end
