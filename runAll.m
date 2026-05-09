function runAll(cfg)
%RUNALL  Master driver. Produces every cached `.mat` artifact and
% rendered table / figure consumed by the reply manuscript.
%
%   runAll()                 % uses default_config()
%   runAll(cfg)
%
% Pipeline order matches the dependency DAG: independent producers come
% first, then the hyperparameter sweep + adaptive selection, then
% bootstraps and asset-pricing sims, transaction-cost variants for Table
% A.3, the per-table aggregates, and finally the render layer plus smoke
% tests.
%
% Runtime estimate (full clean rebuild on an 8-core machine): ~10-20 hours.
% Pause Dropbox sync first; the multi-hour stages write thousands of small
% `.mat` files. Per-stage estimates are in README.md.

if nargin < 1 || isempty(cfg); cfg = default_config(); end
setup_paths();

% --- Daily core analyses (g=2 baseline; g=1 / g=3 robustness panels)
dailyEmpirics(cfg);                              % results/daily/
dailyEmpiricsGPriorRobustness(cfg);              % results/gPriorRobustness/{g1,g3}/

% --- Monthly + Fama-French
monthlyEmpirics(cfg);                            % results/monthly/
famaFrenchEmpirics(cfg, {'SMB', 'HML'});         % results/famaFrench/

% --- Hyperparameter sweep + adaptive selection.
% signSpec = 1 (unrestricted) and 2 (sign-restricted; published baseline).
dailyEmpiricsHyperparameters1(cfg);              % results/hyperparameters/
dailyEmpiricsHyperparameters2(cfg);              % results/hyperparameters/
dailyEmpiricsHyperparametersMarginals(cfg, 1);   % results/aggregates/topKbotK_signSpec1.mat
dailyEmpiricsHyperparametersMarginals(cfg, 2);   % results/aggregates/topKbotK.mat (default)

% --- Bootstraps and asset-pricing sims.
% data/csv_sims/{BY,CC,DT,GP,W}_sim.csv are external simulation outputs
% (Bansal-Yaron, Campbell-Cochrane, Di Tella, Garleanu-Panageas,
% Wachter); they are inputs, not regenerated here.
dailyBootstrap(cfg);                             % results/bootstrap/
generateStickyExpectationsData(cfg);             % results/simulatedPaths/
stickyExpectationsSim(cfg);                      % results/assetPricing/
dailyAssetPricingBootstrap(cfg);                 % results/assetPricing/

% --- Transaction-cost variants for Table A.3 (5 bps and 10 bps).
for bps = [5, 10]
    fprintf('runAll: transaction-cost variant adjCostBps=%d\n', bps);
    cfgBps = cfg;
    cfgBps.adjCostBps    = bps;
    cfgBps.resultsSubdir = sprintf('daily_%dbps', bps);
    dailyEmpirics(cfgBps);
    cfgBpsAdapt = cfgBps;
    cfgBpsAdapt.sampleSplit = 'full';
    [~, alphMatAdaptive, tStatAlphMatAdaptive, SRAdaptive] = ...
        adaptivePanelsForSplit(cfgBpsAdapt);
    aggDir = fullfile(cfg.paths.results, 'aggregates');
    save(fullfile(aggDir, sprintf('cwBestMat_full_%dbps.mat', bps)), ...
        'alphMatAdaptive', 'tStatAlphMatAdaptive', 'SRAdaptive');
end

% --- Per-table aggregates consumed by the table scripts.
computeSignifCell(cfg);                          % results/aggregates/signifCell.mat
computeTab1Results(cfg);                         % results/aggregates/tab1Results.mat
computeTabA3Results(cfg);                        % results/aggregates/tabA3Results.mat

% --- Render and verify
displayResults();
testsDir = fullfile(cfg.paths.root, 'tests');
run(fullfile(testsDir, 'run_smoke_tests.m'));
end
