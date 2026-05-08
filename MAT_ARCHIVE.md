# Cached `.mat` archive

This package's `results/` directory holds the cached MATLAB outputs
that `displayResults` consumes to produce every reply table and figure.
Replicators have two options:

1. **Rebuild from scratch** with `runAll(default_config())` — see the
   main `README.md` for runtime estimates (~10-20 hours on an 8-core
   machine).
2. **Download the cached archive(s)** described below and skip the
   compute, going straight to `displayResults`.

## Two archives

| File | Approx. size | Contents | When you need it |
|---|---|---|---|
| `pop_results_2026-05.zip` | ~0.8 GB | Pruned `results/` — only the `.mat` files actually consumed by `displayResults`. | Always (required for rendering). |
| `pop_data_external_2026-05.zip` | ~10-15 GB | `data/pcPredictor.mat` (2.5 GB) plus `data/csv_sims/` (19 GB raw, ~10-13 GB after zip compression). | Only if you want to rerun the empirical pipeline from scratch. Not needed for rendering tables and figures. |

The `2026-05` suffix is the cache snapshot date; bump it whenever you
publish a new archive so replicators can tell which version they have.

## Download

The archives are hosted at:

> **`[DROPBOX URL FOR pop_results_2026-05.zip TBD]`**
>
> **`[DROPBOX URL FOR pop_data_external_2026-05.zip TBD]`**

Replace these placeholders with the published Dropbox links once the
zips are uploaded.

## Creating the archives (author workflow)

The archives are produced from the package's local `results/` and
`data/` directories after a clean production run.

**Step 0 — Verify the pruned `results/`.**

`results/` should already contain only the cached `.mat` files that
`displayResults` consumes (~30 files, ~810 MB). If the directory is
larger, prune intermediates that `displayResults` does not read before
zipping; the typical drop is from ~98 GB to ~810 MB after removing
hyperparameter-sweep and bootstrap intermediates that the rendering
pipeline never opens.

```matlab
cd 'Replication Package'
verify_artifacts
```

Should print PASS for every expected file.

**Step 1 — Zip the pruned `results/`.**

```powershell
cd "Replication Package"
$dateTag = (Get-Date -Format 'yyyy-MM')
$resultsZip = "..\pop_results_$dateTag.zip"
Compress-Archive -Path .\results -DestinationPath $resultsZip -CompressionLevel Optimal
Get-Item $resultsZip | Format-Table Name, @{N='MB';E={[math]::Round($_.Length/1MB,1)}}
```

**Step 2 — Zip the external data inputs (only if distributing the full
data archive).**

```powershell
cd "Replication Package"
$dataZip = "..\pop_data_external_$dateTag.zip"
Compress-Archive -Path .\data\pcPredictor.mat, .\data\csv_sims `
    -DestinationPath $dataZip -CompressionLevel Optimal
Get-Item $dataZip | Format-Table Name, @{N='MB';E={[math]::Round($_.Length/1MB,1)}}
```

**Step 3 — Upload to Dropbox.**

Move both zips into a Dropbox folder (e.g.,
`Dropbox/PopReplicationArchive/`), right-click each file in the Dropbox
web UI or desktop client, and choose "Copy link". Paste the resulting
URLs into the placeholders at the top of this file and into the
"Cached `.mat` archive" section of `README.md`.

The Dropbox UI will produce links of the form
`https://www.dropbox.com/scl/fi/<id>/<filename>?rlkey=<key>&dl=0`. The
trailing `dl=0` opens the preview page; replace with `dl=1` if you want
links that trigger a direct download.

## Naming convention

- `pop_results_<YYYY-MM>.zip` — pruned `results/` snapshot.
- `pop_data_external_<YYYY-MM>.zip` — external data inputs.
- `pop_full_<YYYY-MM>.zip` (optional) — both of the above merged, if
  you want to distribute a single download.

`pop` = Pockets of Predictability. The date suffix matches the cache
snapshot, not the publication date of the reply.

## Replicator workflow (after download)

```powershell
# From the package root:
cd "Replication Package"

# Extract pop_results_<date>.zip — produces results/
Expand-Archive -Path ..\pop_results_2026-05.zip -DestinationPath .

# (Optional) extract pop_data_external_<date>.zip — produces
# data/pcPredictor.mat and data/csv_sims/
Expand-Archive -Path ..\pop_data_external_2026-05.zip -DestinationPath .\data
```

Then:

```matlab
cd 'Replication Package'
verify_artifacts        % confirm every required .mat is present
displayResults          % render every reply table and figure (~2 min)
```

## Notes

- `pop_data_external_<date>.zip` is only needed for
  `runAll(default_config())`. If you stick to `displayResults`, you
  can skip it.
- If you fork this package and produce your own results, mirror the
  same pruned subdirectory structure and bump the date suffix.
