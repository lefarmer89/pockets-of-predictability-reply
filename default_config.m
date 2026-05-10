function cfg = default_config()
% DEFAULT_CONFIG  Return the parameter struct that reproduces the published
% reply exactly. Every entry-point function accepts cfg = default_config()
% as its baseline argument; alternative settings are for R1 robustness
% investigations.
%
%   cfg = default_config()
%
% Fields:
%   paths              struct returned by setup_paths()
%   coefWindowYears    rolling window for the time-varying coefficient
%                      kernel regression (years)
%   sedWindowYears     rolling window for the SED pocket-detection
%                      regression (years)
%   weightWindowYears  rolling window for the G-prior shrinkage weight
%                      estimation (years)
%   minPocketDays      minimum number of trading days a pocket must persist
%                      to be classified
%   gPrior             G-prior parameter (g)
%   shrinkageTarget    'zero' | 'prevailingMean'  (R1 hook)
%   winsorPct          recursive winsorization percentile (both tails)
%   weightLB           portfolio-weight lower bound
%   weightUB           portfolio-weight upper bound
%   adjCostBps         per-trade transaction cost in basis points
%   sampleSplit        'full' | 'pre1989' | 'post1989'  (only consumed by
%                      computeAdaptiveSampleSplit and adaptivePanelsForSplit
%                      for ex-post date masking; no fresh empirical runs
%                      are produced at the subsample level)
%   selectionCriterion 'maxAlpha' | 'tAlpha' | 'rmse'  (R1 hook for the
%                      Adaptive panel's combo selection at each t.
%                      'maxAlpha' = argmax of recursively-computed alpha
%                                   over the OOS sample (the published
%                                   "Adaptive (max alpha)" column, default).
%                      'tAlpha'  = argmax of HAC t(alpha).
%                      'rmse'    = argmin of cumulative RMSE
%                                  (the "Adaptive (min RMSE)" column).
%                      The active rendering pipeline always emits both
%                      maxAlpha and rmse columns; this hook controls
%                      which one is treated as the "default" alongside
%                      filename suffixes for alternative-criterion runs.)
%   resultsSubdir      override target subfolder under cfg.paths.results
%                      for entry-point .mat saves. Empty (default) routes
%                      to the canonical kind-specific subfolder; set to
%                      e.g. fullfile('gPriorRobustness','g1') for the
%                      g-prior robustness wrapper.
%   hyperparameterComboSubset
%                      [] = use all 9720 combos (default). Set to a
%                      vector of row indices into the hyperparameter
%                      grid for shakeout / profiling on a small subset
%                      (honored by dailyEmpiricsHyperparameters1,
%                      dailyEmpiricsHyperparameters2, and
%                      dailyEmpiricsHyperparametersMarginals).
%   bootstrapReps      [] = use the entry-point default (500 for
%                      stickyExpectationsSim, 1000 for dailyBootstrap and
%                      dailyAssetPricingBootstrap). Set to a smaller
%                      integer (e.g., 10) for shakeout runs.
%   hyperparameterSpecs
%                      [] (default) = run signSpec ∈ {1, 2} when
%                      dailyEmpiricsHyperparameters1/2/Marginals is invoked.
%                      Set to a scalar (e.g. 2) to run only the
%                      published-baseline signSpec.
%   dailySpecsSubset   Default [1, 6, 7, 8, 9, 10, 11] runs only the 7
%                      spec rows whose outputs displayResults consumes:
%                      row 1 (sgn=1 canonical), 6-10 (sgn=2 with all 5
%                      kernel/DM combos for the window-robustness sweep),
%                      and 11 (sgn=3 canonical). Set to a vector of row
%                      indices for shakeout (e.g. [7]) or to 1:15 for the
%                      full sweep.
%   monthlySpecsSubset [] (default) = run all 3 specs in monthlyEmpirics.
%                      Vector of row indices for shakeout.
%   useParallel        [] (auto-detect via useParfor) | true | false
%   marginalsParforWorkers
%                      [] (default = 6) cap on parfor workers for
%                      dailyEmpiricsHyperparametersMarginals only. The
%                      Marginals stage broadcasts ~620 MB per worker
%                      and OOMs at 12 workers on 32 GB machines. Other
%                      stages run at full pool size. The function
%                      shrinks the pool on entry and restores it on
%                      exit via onCleanup.
%   rngSeed            integer seed for any RNG-dependent step
%
% Defaults below match the published reply.

cfg = struct();
cfg.paths = setup_paths();

cfg.coefWindowYears   = 2.5;
cfg.sedWindowYears    = 1;
cfg.weightWindowYears = 1;
cfg.minPocketDays     = 21;

cfg.gPrior            = 2;
cfg.shrinkageTarget   = 'zero';
cfg.winsorPct         = 2.5;

cfg.weightLB          = 0;
cfg.weightUB          = 2;
cfg.adjCostBps        = 0;

cfg.sampleSplit        = 'full';
cfg.selectionCriterion = 'maxAlpha';
cfg.resultsSubdir             = '';
cfg.hyperparameterComboSubset = [];
cfg.hyperparameterSpecs       = [];
cfg.bootstrapReps             = [];
cfg.dailySpecsSubset          = [1, 6, 7, 8, 9, 10, 11];
cfg.monthlySpecsSubset        = [];

cfg.useParallel            = [];
cfg.marginalsParforWorkers = [];
cfg.rngSeed                = 20240101;
end
