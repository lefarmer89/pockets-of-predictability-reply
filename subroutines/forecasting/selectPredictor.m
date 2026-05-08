function predictors = selectPredictor(varName, data, paths)
% SELECTPREDICTOR  Pick the predictor matrix or column for a given
% variable name, replacing the long switch-case duplicated across the
% empirics functions.
%
%   predictors = selectPredictor(varName, data, paths)
%
% varName is one of:
%   'dp', 'tbl', 'tsp', 'rvar' -- a single column from data
%   'mv'                       -- the multivariate matrix [dp, tbl, tsp, rvar]
%   'pc'                       -- the pre-computed first principal
%                                 component, loaded from data/pcPredictor.mat
%   'erL'                      -- excess returns (lagged in the caller)
%
% data is a table with the corresponding columns. paths is the struct
% returned by setup_paths(); paths.data is used to locate pcPredictor.mat.

switch varName
    case 'dp'
        predictors = data.dp;
    case 'tbl'
        predictors = data.tbl;
    case 'tsp'
        predictors = data.tsp;
    case 'rvar'
        predictors = data.rvar;
    case 'mv'
        predictors = [data.dp, data.tbl, data.tsp, data.rvar];
    case 'pc'
        loaded = load(fullfile(paths.data, 'pcPredictor.mat'));
        predictors = loaded.predictors;
    case 'erL'
        predictors = data.exret;
    otherwise
        error('selectPredictor:unknownVariable', ...
            'Unknown variable name "%s"', varName);
end
end
