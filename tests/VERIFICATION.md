# Verifying the refactored package against the original

The refactor preserved every entry-point function's contract: same inputs,
same saved-variable layout, same numerical algorithm. To confirm
numerical equivalence, re-run each entry point and diff its outputs
against the snapshot in
`Archive/Replication Package - Original Backup/`.

## What you need

- MATLAB R2026a (R2024a or newer should work; R2026a is what the smoke
  tests run against).
- Cached predictor data in `Replication Package/data/` (already present).
- The frozen original outputs in
  `Archive/Replication Package - Original Backup/{Baseline_Results,
  Bootstrap_Daily, Asset_Pricing_Results, Robustness_Results}/`.

The verification tool is:
```
Replication Package/tests/verify_against_archive.m
```

It loads every cached `.mat` produced by an entry point, walks each
variable, and reports `max abs diff` against the same file in the
archive. Tolerance defaults to `1e-10`.

## Sanity check first (no compute needed)

Confirm the harness wiring works by running it against the **current
cached results** (which still come from the original code if you have
not re-run anything yet — both should match exactly):

```bash
cd 'Replication Package/tests'
"C:\Program Files\MATLAB\R2026a\bin\matlab.exe" -batch \
  "cd('..'); addpath(pwd); setup_paths(); cd tests; \
   verify_against_archive('monthlyEmpirics');"
```

You should see `87/87 variables match within tol = 1.0e-10`. If you
don't, either (a) `Archive/Replication Package - Original Backup/` is
missing, (b) the current `results/Baseline_Results/` got partially
re-run already, or (c) the harness can't read the .mat files.

## Verification workflow per entry point

For each entry point:

1. **Optional but recommended**: snapshot the current result subdirectory
   so you can roll back if numbers drift.
   ```bash
   cp -r 'Replication Package/results/Baseline_Results' \
         'Replication Package/results/Baseline_Results.bak'
   ```

2. **Re-run** the refactored entry point. This overwrites the per-spec
   `.mat` files in the corresponding `results/<subdir>/`.
   ```bash
   "C:\Program Files\MATLAB\R2026a\bin\matlab.exe" -batch \
     "cd('Replication Package'); addpath(pwd); setup_paths(); \
      monthlyEmpirics(default_config());"
   ```

3. **Verify**:
   ```bash
   "C:\Program Files\MATLAB\R2026a\bin\matlab.exe" -batch \
     "cd('Replication Package'); addpath(pwd); setup_paths(); \
      cd tests; verify_against_archive('monthlyEmpirics');"
   ```

4. **If FAIL**, the failing rows print as
   `FAIL  <fname>  (X/Y vars; max abs diff = …)`. To see per-variable
   detail add `'verbose', true`:
   ```matlab
   r = verify_against_archive('monthlyEmpirics', 'verbose', true);
   ```

5. **If PASS**, optionally delete the `.bak` snapshot.

## Suggested order (cheapest → most expensive)

| # | Entry point                    | Wall-clock | Comment                          |
|---|--------------------------------|-----------|-----------------------------------|
| 1 | `monthlyEmpirics`              | ~30 min × 3 specs | Smallest data, 3 specs.    |
| 2 | `famaFrenchEmpirics('SMB')`    | ~30 min   | Single spec, single factor.       |
| 3 | `famaFrenchEmpirics('HML')`    | ~30 min   | Same.                              |
| 4 | `dailyEmpirics`                | ~2-4 h    | 15 specs over the daily sample.   |
| 5 | `dailyAssetPricingBootstrap`   | ~2-4 h × 6 models | Longest in elapsed time. |
| 6 | `dailyBootstrap`               | ~6-10 h × 3 specs | Longest single-call.     |
| 7 | `dailyEmpiricsHyperparameters1` | ~6-12 h  | 4860-combo sweep; needs the prevInd cache to stay within budget. |
| 8 | `dailyEmpiricsHyperparameters2` | ~2-4 h   | Post-processing of #7 output.     |
| 9 | `dailyEmpiricsHyperparametersMarginals` | ~30-60 min | Loads `topKbotK.mat` (saved by #8) by default. |

You can stop at any tier. The Tier-1 baseline (#1-#4) covers Tables 1,
2, 3, A.1, A.2, A.3, A.4, A.5, A.6 and Figure 2 — the entire
empirical layer of the reply. The bootstrap and hyperparameter tiers
(#5-#9) feed Tables 2 (significance markers), 4, and Figures 1, A.1,
A.2, A.3.

## Quick spot-check (no full re-run)

If you only want a feel for whether anything moved, skip the
overnight runs and just check that **`displayResults` still produces
the published numbers** from the cached files (which were produced by
the original code, so this is purely a smoke check on the rendering
pipeline that table*.m and figure*.m use):

```bash
cd 'Replication Package/tests'
"C:\Program Files\MATLAB\R2026a\bin\matlab.exe" -batch \
  "addpath(pwd); run_smoke_tests"
```

The `test_displayResults` step inside the suite runs the orchestrator
end-to-end and exits non-zero if any table or figure errors out. This
takes ~30-50 s and currently passes 4/4. It does NOT verify the
refactored compute paths — it only confirms the rendering layer still
loads the cached `.mat` files and produces the LaTeX tables and EPS
figures.

## Sentinel values to eyeball

The CLAUDE.md convention is to spot-check Table 1 row 1-2 (the `dp` and
`tbl` α and t(α) values for the published baseline spec). After running
`dailyEmpirics`, render Table 1:

```matlab
folder = fullfile(cfg.paths.results, 'Baseline_Results');
table1     % runs tables/table1.m
```

The `dp` and `tbl` rows should match the values in the published reply
(see `Reply_to_Replication_Study_FST.pdf`, Table 1, panel A — "Fixed
shrinkage" columns).

## What "PASS" really means

`verify_against_archive` reports PASS when `max(abs(current - archive))
<= tol` for every numeric variable in every cached file, AND every NaN
position matches between current and archive, AND every struct/cell
field structure matches.

What it does NOT check:
  - Files in current that are NOT in the archive (skipped — these
    are new outputs; e.g., the refactored `dailyEmpirics` saves an
    extra `dmMat` field in `OOSResults_*.mat` that the original did
    not, so that variable is reported as `new variable (archive lacks
    it)` and counted as PASS by convention).
  - Files in the archive that are NOT in current (reported as 'missing
    in current' — would PASS only if the file was intentionally
    dropped from the refactor).
  - Behavior under non-default `cfg` (the archive was generated with
    `default_config()`; verification of R1 hooks like
    `cfg.shrinkageTarget = 'prevailingMean'` or
    `cfg.sampleSplit = 'pre1989'` requires generating fresh
    archive-equivalents first).
