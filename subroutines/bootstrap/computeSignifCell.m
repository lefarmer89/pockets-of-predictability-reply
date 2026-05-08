function computeSignifCell(cfg)
%COMPUTESIGNIFCELL  Compute per-pocket significance flags for the four
% univariate predictors used in Table 2.
%
%   computeSignifCell(cfg)
%
% A pocket is flagged significant (= 2) when its integral R^2 has a
% bootstrap p-value <= 5% under the IID-bootstrap-with-EGARCH-residuals
% null distribution. Indices 1-4 use each variable's own bootstrap
% distribution; index 5 reuses the tbl bootstrap.
%
% Reads:  results/daily/forecastResults_1_2.5yE_1yDM_1S_pm.mat
%         results/bootstrap/bootstrapResultsOOS_*.mat (per variable)
% Writes: results/aggregates/signifCell.mat
% Consumers: tables/table2.m.

if nargin < 1 || isempty(cfg); cfg = default_config(); end

dailyFolder = fullfile(cfg.paths.results, 'daily');
bootstrapFolder = fullfile(cfg.paths.results, 'bootstrap');
aggDir = fullfile(cfg.paths.results, 'aggregates');
if ~isfolder(aggDir); mkdir(aggDir); end

F = load(fullfile(dailyFolder, 'forecastResults_1_2.5yE_1yDM_1S_pm.mat'), ...
    'integralR2Mat', 'signRestriction', 'coefRestriction');

varNames   = {'dp', 'tbl', 'tsp', 'rvar'};
signifCell = cell(5, 1);

for ii = 1:5
    bsVar = varNames{min(ii, 4)};   % index 5 reuses tbl bootstrap (legacy)
    if ii == 5; bsVar = varNames{2}; end
    bsFile = fullfile(bootstrapFolder, sprintf( ...
        'bootstrapResultsOOS_signRestriction_%d_coefRestriction_%d_%s_iidbs_egarch_25y_0.mat', ...
        F.signRestriction, F.coefRestriction, bsVar));
    BS = load(bsFile, 'integralR2');

    integralR2Cur = F.integralR2Mat(~isnan(F.integralR2Mat(:, ii)), ii);
    pVals = NaN(size(integralR2Cur));
    bsValues  = BS.integralR2(:);
    bsDenom   = sum(~isnan(bsValues));
    parfor (kk = 1:numel(integralR2Cur), useParfor(cfg))
        pVals(kk) = sum(bsValues > integralR2Cur(kk)) / bsDenom;
    end

    plotSignif = ones(size(integralR2Cur));
    plotSignif(pVals <= 0.05) = 2;
    plotSignif(pVals  > 0.05) = 0;
    signifCell{ii} = plotSignif;
end

save(fullfile(aggDir, 'signifCell.mat'), 'signifCell');
fprintf('  wrote %s\n', fullfile(aggDir, 'signifCell.mat'));
end
