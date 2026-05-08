# Replication package: reply to Cakici, Fieberg, Neumaier, Poddig, and Zaremba (2024)

This package reproduces the tables and figures in the authors' reply to
the replication comment on Farmer, Schmidt, and Timmermann (2023,
*Journal of Finance*), "Pockets of Predictability." The reply identifies
a coding error in the original paper (a two-sided kernel that allowed
future information into the pocket-classification step) and re-runs the
analysis with a corrected one-sided kernel and a G-prior Bayesian
shrinkage of the time-varying coefficient forecasts.

## Citation

```bibtex
@article{FarmerSchmidtTimmermann2023,
  title   = {Pockets of Predictability},
  author  = {Farmer, Leland E. and Schmidt, Lawrence and Timmermann, Allan},
  journal = {Journal of Finance},
  year    = {2023},
  volume  = {78},
  number  = {3},
  pages   = {1279--1341}
}
```

## System requirements

- MATLAB R2021a or later (developed and tested on R2026a).
- Statistics and Machine Learning Toolbox (used by `regstats2.m` for HAC
  standard errors).
- Parallel Computing Toolbox (optional). Without a license,
  `subroutines/utils/useParfor.m` returns 0 and every `parfor` falls back
  to a serial `for` without code changes.
- ~32 GB RAM minimum for the empirics and bootstraps; ~64 GB recommended
  for the full hyperparameter pipeline. ~80 GB free disk for cached
  outputs and simulation CSVs.
- Windows: enable long-path support
  (`HKLM\SYSTEM\CurrentControlSet\Control\FileSystem\LongPathsEnabled = 1`)
  to handle deep paths under `results/assetPricing/`.

## Quick start

To render the reply tables and figures without rerunning the pipeline:

1. Download the cached `.mat` archive (`results.zip`) from the URL in
   `MAT_ARCHIVE.md` and extract it at the package root so the contents
   land under `results/`.
2. Verify the archive:

    ```matlab
    cd 'Replication Package'
    verify_artifacts
    ```

3. Render the tables and figures:

    ```matlab
    displayResults
    ```

    Runtime: ~2 minutes. `.tex` files land in `output/tables/`, `.eps`
    files in `output/figures/`.

To run the smoke tests:

```matlab
cd 'Replication Package/tests'
run_smoke_tests
```

To rebuild from scratch:

```matlab
cd 'Replication Package'
runAll(default_config())
```

Runtime: ~10-20 hours on an 8-core machine.

## Folder structure

```
Replication Package/
├── README.md
├── MAT_ARCHIVE.md                  pointer to the cached .mat zip
├── .gitignore
│
├── runAll.m                        master driver: full pipeline + display + smoke
├── displayResults.m                render all tables and figures from cached .mat
├── verify_artifacts.m              confirm every expected .mat is present
├── default_config.m                cfg struct with reply-replicating defaults
├── setup_paths.m                   idempotent path helper
│
├── dailyEmpirics.m                 daily out-of-sample forecasting (Tables 1, 2, 3, A.1, A.2, A.3 Fixed; Figure 2)
├── dailyEmpiricsGPriorRobustness.m wraps dailyEmpirics for g=1 and g=3 (Table A.1)
├── dailyEmpiricsHyperparameters1.m forecasts over the 9,720-combo grid
├── dailyEmpiricsHyperparameters2.m per-signSpec OOS results
├── dailyEmpiricsHyperparametersMarginals.m  adaptive selection (Figures 1, A.1-A.3)
├── monthlyEmpirics.m               monthly variants (Tables A.5, A.6)
├── famaFrenchEmpirics.m            SMB and HML factor portfolios (Table A.4)
├── dailyBootstrap.m                bootstrap for pocket significance (Table 2 markers)
├── dailyAssetPricingBootstrap.m    bootstrap on five asset-pricing models (Table 4)
├── generateStickyExpectationsData.m simulate sticky / rational paths
├── stickyExpectationsSim.m         pocket detection on simulated paths (Table 5)
├── PathSimulator.m                 path-simulation utility
│
├── tables/                         table1..5, tableA1..A6
├── figures/                        figure1, figure2, figureA1..A3
├── subroutines/
│   ├── forecasting/   kernel + G-prior + per-variable estimation
│   ├── diagnostics/   pocket detection + DM/CW test stats
│   ├── portfolio/     portfolio weights, performance regressions
│   ├── bootstrap/     opt_block_length (Andrews 1991) + significance helpers
│   ├── hyperparameters/ hyperparameter grid; adaptive selection helpers
│   ├── io/            path resolution, table rendering, file naming
│   └── utils/         kernel and regression building blocks
│
├── data/                           read-only inputs
│   ├── Daily_Predictors.xlsx       dp, tbl, tsp, rvar, mv, exret, rf, ...
│   ├── Monthly_Predictors.xlsx     monthly counterparts
│   ├── pcPredictor.mat             cached principal-component predictor
│   ├── quarterly_fr_coeffs.csv     Coibion-Gorodnichenko forecast-error coefs
│   ├── quarterly_avg_forecast_revsanderrs.csv
│   ├── rationalCalibration_freevols_r211.mat   sticky-expectations calibration
│   ├── stickyCalibration_freevols_r211.mat
│   └── csv_sims/                   asset-pricing simulation outputs
│
├── results/                        cached pipeline outputs (.mat; .gitignored)
│   ├── daily/                      dailyEmpirics outputs
│   ├── monthly/                    monthlyEmpirics outputs
│   ├── famaFrench/                 famaFrenchEmpirics outputs
│   ├── hyperparameters/            Hyperparameters1/2 outputs
│   ├── bootstrap/                  dailyBootstrap outputs
│   ├── assetPricing/               dailyAssetPricingBootstrap + stickyExpectationsSim outputs
│   ├── simulatedPaths/             generateStickyExpectationsData outputs
│   ├── gPriorRobustness/{g1,g3}/   dailyEmpiricsGPriorRobustness outputs
│   └── aggregates/                 cross-driver aggregate caches
│
├── output/
│   ├── tables/                     rendered .tex
│   ├── figures/                    rendered .eps
│   └── logs/                       diary files from long runs
│
└── tests/
    ├── run_smoke_tests.m
    ├── test_displayResults.m
    ├── test_infrastructure.m
    ├── test_results_validation.m
    ├── test_subroutines.m
    ├── test_yF2AllS_bugfix.m
    └── VERIFICATION.md
```

## Pipeline

| Entry point | Reads | Writes | Produces |
|---|---|---|---|
| `dailyEmpirics(cfg)` | `data/Daily_Predictors.xlsx`, `data/pcPredictor.mat` | `results/daily/`, `results/aggregates/tab1Results.mat`, `results/aggregates/tabA3Results.mat` | Tables 1, 2, 3, A.1, A.2, A.3 (Fixed panels); Figure 2 |
| `dailyEmpiricsGPriorRobustness(cfg)` | same as `dailyEmpirics` | `results/gPriorRobustness/{g1,g3}/` | Table A.1 (g=1 and g=3 columns) |
| `dailyEmpiricsHyperparameters1(cfg)` | same | `results/hyperparameters/forecastResults_*_HyperR1.mat` | OOS forecasts across the 9,720-combo grid |
| `dailyEmpiricsHyperparameters2(cfg)` | `results/hyperparameters/forecastResults_*_HyperR1.mat` | `results/hyperparameters/OOSResults_*_HyperR1_ExpandingC.mat` | Per-signSpec OOS results consumed by Marginals |
| `dailyEmpiricsHyperparametersMarginals(cfg)` | `results/hyperparameters/` | `results/aggregates/topKbotK*.mat`, `oosCAlphas*.mat` | Figures 1, A.1-A.3; Adaptive panels of Tables 1 and A.3 |
| `monthlyEmpirics(cfg)` | `data/Monthly_Predictors.xlsx` | `results/monthly/` | Tables A.5, A.6 |
| `famaFrenchEmpirics(cfg)` | `data/Daily_Predictors.xlsx` | `results/famaFrench/` | Table A.4 |
| `dailyBootstrap(cfg)` | `data/Daily_Predictors.xlsx` | `results/bootstrap/`, `results/aggregates/signifCell.mat` | Significance markers in Table 2 |
| `dailyAssetPricingBootstrap(cfg)` | `data/csv_sims/{BY,CC,DT,GP,W}_sim.csv` | `results/assetPricing/` | Table 4 |
| `generateStickyExpectationsData(cfg)` | calibration `.mat` files | `results/simulatedPaths/` | Inputs to `stickyExpectationsSim` |
| `stickyExpectationsSim(cfg)` | `results/simulatedPaths/` | `results/assetPricing/SE_*.mat`, `RE_*.mat`, `RE_Recalibrated_*.mat` | Table 5 |
| `displayResults` | every cached `.mat` | `output/tables/*.tex`, `output/figures/*.eps` | All reply tables and figures |

## Reproducing specific results

| Reply | Script | Description |
|---|---|---|
| Table 1 | `tables/table1.m` | Out-of-sample forecasting performance (Fixed + Adaptive panels) |
| Table 2 | `tables/table2.m` | Pocket statistics (daily) |
| Table 3 | `tables/table3.m` | Sign-restriction comparison |
| Table 4 | `tables/table4.m` | Asset-pricing model OOS simulations |
| Table 5 | `tables/table5.m` | Sticky-expectations model simulations |
| Table A.1 | `tables/tableA1.m` | G-prior parameter robustness (g ∈ {1, 2, 3}) |
| Table A.2 | `tables/tableA2.m` | Window-length robustness |
| Table A.3 | `tables/tableA3.m` | Transaction-cost panels (0, 5, 10 bps) |
| Table A.4 | `tables/tableA4.m` | Fama-French SMB and HML |
| Table A.5 | `tables/tableA5.m` | Monthly pocket statistics |
| Table A.6 | `tables/tableA6.m` | Monthly out-of-sample performance |
| Figure 1 | `figures/figure1.m` | Best-100 adaptive shrinkage forecasts (t-stats and CW) |
| Figure 2 | `figures/figure2.m` | Coibion-Gorodnichenko forecast-error correlations |
| Figure A.1 | `figures/figureA1.m` | Hyperparameter histograms (α t-stats) |
| Figure A.2 | `figures/figureA2.m` | Hyperparameter histograms (in-pocket CW) |
| Figure A.3 | `figures/figureA3.m` | Hyperparameter scatter (ex-ante vs ex-post c selection) |

## Configuration

`default_config()` returns a `cfg` struct with every tunable parameter
at its reply-replicating default:

| Field | Default | Effect |
|---|---|---|
| `cfg.coefWindowYears` | `2.5` | rolling window for the time-varying coefficient kernel |
| `cfg.sedWindowYears` | `1` | rolling window for the SED pocket-detection regression |
| `cfg.weightWindowYears` | `1` | rolling window for the G-prior shrinkage weight |
| `cfg.minPocketDays` | `21` | minimum pocket duration (trading days) |
| `cfg.gPrior` | `2` | Bayesian shrinkage parameter g |
| `cfg.shrinkageTarget` | `'zero'` | `'zero'` or `'prevailingMean'` |
| `cfg.winsorPct` | `2.5` | recursive winsorization percentile (both tails) |
| `cfg.weightLB`, `cfg.weightUB` | `0`, `2` | portfolio-weight bounds |
| `cfg.adjCostBps` | `0` | per-trade transaction cost (basis points) |
| `cfg.selectionCriterion` | `'maxAlpha'` | adaptive-shrinkage selector: `'maxAlpha'`, `'tAlpha'`, or `'rmse'` |
| `cfg.hyperparameterComboSubset` | `[]` (all 9,720) | row-index vector for shakeout subsets |
| `cfg.hyperparameterSpecs` | `[]` (signSpec ∈ {1, 2}) | restrict signSpec sweep |
| `cfg.bootstrapReps` | `[]` (entry-point default) | smaller value for shakeout |
| `cfg.useParallel` | `[]` (auto) | force `parfor` on/off |
| `cfg.rngSeed` | `20240101` | RNG seed for bootstrap reproducibility |

## Data provenance

`data/Daily_Predictors.xlsx` and `data/Monthly_Predictors.xlsx` follow
the data construction in Welch and Goyal (2008) and in the original
2023 paper. Columns include the dividend-price ratio (`dp`), the
3-month T-bill rate (`tbl`), the term spread (`tsp`), realized variance
(`rvar`), excess returns, the risk-free rate, and standard recession /
business-cycle indicators.

`data/quarterly_fr_coeffs.csv` and
`data/quarterly_avg_forecast_revsanderrs.csv` hold Coibion-Gorodnichenko
(2015) forecast-error coefficients computed from SPF data; they feed
`figures/figure2.m` only.

`data/csv_sims/{BY,CC,DT,GP,W}_sim.csv` are external simulation outputs
from the canonical asset-pricing models (Bansal-Yaron, Campbell-Cochrane,
Duffie-Tomz, Gabaix-Postemski, Wachter). They are inputs to
`dailyAssetPricingBootstrap.m` and **are not regenerated by any script
in this package**. Each CSV is several gigabytes; together they total
~19 GB and are distributed via the cached `.mat` archive (see
`MAT_ARCHIVE.md`).

`data/rationalCalibration_freevols_r211.mat` and
`data/stickyCalibration_freevols_r211.mat` hold the solved-out
parameters for the rational and sticky-expectations versions of the
Bansal-Yaron-style model in Section V of the reply.

## Reproducing from scratch

`runAll(cfg)` walks the dependency DAG and produces all cached `.mat`
artifacts and rendered exhibits. With defaults:

```matlab
cd 'Replication Package'
runAll(default_config())
```

Per-stage runtimes on an 8-core machine, estimated from the most recent
shakeout. Actual times vary with hardware and parallel-worker count.

| # | Stage | Wall time |
|---|---|---|
| 1 | `dailyEmpirics` | 5-15 min |
| 2 | `dailyEmpiricsGPriorRobustness` (g=1, g=3) | 10-30 min |
| 3 | `monthlyEmpirics` | 1-3 min |
| 4 | `famaFrenchEmpirics` (SMB + HML) | 1-2 min |
| 5 | `dailyEmpiricsHyperparameters1` | 4-8 h |
| 6 | `dailyEmpiricsHyperparameters2` | 1-2 h |
| 7 | `dailyEmpiricsHyperparametersMarginals` (signSpec 1, 2) | 1.5-3 h |
| 8 | `dailyBootstrap` | 3-6 h |
| 9 | `generateStickyExpectationsData` + `stickyExpectationsSim` | 1-2 h |
| 10 | `dailyAssetPricingBootstrap` | 2-4 h |
| 11 | Transaction-cost variants (5 bps, 10 bps) | 10-30 min |
| 12 | `displayResults` | 2-5 min |

For a fast pipeline shakeout that exercises every entry point on
subsetted compute (`cfg.dailySpecsSubset = 1`,
`cfg.hyperparameterComboSubset = 1:10`, `cfg.hyperparameterSpecs = 2`,
`cfg.bootstrapReps = 5`) and writes to `results_smoke/` so the cached
`results/` is untouched, run `tests/smoke_pipeline.m`. The full smoke
takes ~2 hours; the bottleneck is the multi-GB CSV reads in
`dailyAssetPricingBootstrap`, which do not scale with `B`.

## Smoke tests

`tests/run_smoke_tests.m` runs five validation tests:

| Test | Checks |
|---|---|
| `test_displayResults` | `displayResults` runs end-to-end and writes every `.tex` and `.eps` |
| `test_infrastructure` | `setup_paths`, `default_config`, `useParfor` return well-formed structs |
| `test_results_validation` | required `.mat` artifacts exist with the expected variables |
| `test_subroutines` | core subroutines (`regstats2Fast`, `lps1_v2`, `gPriorWeightsFast`, ...) match reference outputs |
| `test_yF2AllS_bugfix` | `yF2AllS` construction in H2 reproduces the cached value within 1e-10 |

## Operational notes

- **Long-running stages.** Stages 1, 2, 5, 6, 7, 8, 9, 10, 11 produce
  thousands of `.mat` writes. On a Dropbox-synced filesystem, pause
  sync first or write outputs to a non-Dropbox directory.
- **Windows MAX_PATH.** Some files in `results/assetPricing/` have
  names long enough to push their full paths past 260 characters.
  Enable long-path support in the registry, or use the `\\?\`
  long-path prefix in scripts that walk the tree.
- **Parfor reproducibility.** Reductions inside `parfor` can change
  summation order, so expect drift around 1e-10 in some cells across
  worker counts. Bootstrap loops pre-generate per-replication seeds
  before each parfor, so the bootstrap distribution is unaffected.

## Cached `.mat` archive

A pruned `results/` zip is distributed via Dropbox; download it,
extract it under `Replication Package/`, and call `displayResults` to
skip the multi-hour rebuild. See `MAT_ARCHIVE.md` for the URL and
extraction instructions.

## Contact

Questions about the replication should go to Leland E. Farmer (UVA).
