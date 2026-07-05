# Data: fetch, don't commit

Raw data is **not** committed (GB-scale HDF5; see `.gitignore`). What IS committed:
`data/zenodo_listings/*.json` (Zenodo file ids + sizes + md5s, snapshotted live
2026-06-10) and `data/manifests/*.json` (extraction manifests: what was extracted,
columns, row counts, sha256s, UTC timestamps, the provenance chain for every derived
parquet). Scientific context: the [README](../README.md) and
[the pre-registration](../PRE_REGISTRATION.md).

All records are open-access **CC-BY-4.0**, anonymous download, $0. Every byte below was
verified against the Zenodo REST API on **2026-06-10**.

## Stage-1 acquisition (2026-06-10): what we hold and why

> **Status at the Stage-1 close-out commit (2026-06-10 evening, UTC-stamped in the
> manifests):** a concurrent workload owned most of the machine's disk, so acquisition ran
> against a hard 3 GB budget with live floor checks; transfers were throttled to
> ~2 MB/s/connection. LANDED: all metadata/event lists; O3-era mixture parquet
> (919,033 × 25); first **7/86** O4a PE events (11 parquet subsets, GPS-time order =
> the pre-registered subset rule). RUNNING UNATTENDED (nohup'd, self-guarding,
> idempotent, re-running the same tier verifies/skips finished work): `gwtc3-plp`
> streaming (~2.2/10.07 GB at commit), `selection`'s 1.44 GB O4a injection set
> (waiting for a disk window, then extract-and-delete), and a sequencer that
> restarts `--tier o4a-pe` the moment selection exits. If any of them died early,
> resume with: `.venv/bin/python -m src.fetch_data --tier <gwtc3-plp|selection|o4a-pe>`,
> each is safe to re-run any number of times. Partial O4a acquisition is loudly
> documented in `manifests/o4a_pe_acquisition_log.json` (per-event timestamps).

Executed under a **hard 3 GB new-files budget enforced in code** (`src/fetch_data.py`:
`budget_check()` refuses any download/extraction that would exceed 3.0 GB held or leave
< 1.2 GB free disk). Discipline: download → extract needed columns → **delete raw**.

```bash
.venv/bin/python -m src.fetch_data --plan        # tiers + live budget status
.venv/bin/python -m src.fetch_data --tier metadata   # listings, event lists, docs (~25 MB)
.venv/bin/python -m src.fetch_data --tier selection  # injections -> parquet, raw deleted
.venv/bin/python -m src.fetch_data --tier gwtc3-plp  # STREAM 10.07 GB monolith, keep PLPeak only
.venv/bin/python -m src.fetch_data --tier o4a-pe     # BLIND per-event subsets, one at a time
```

### Integrity firewall (amended 2026-06-10, see PRE_REGISTRATION_DRAFT.md header)

GWTC-4.0 per-event posteriors were acquired **pre-lock** under the amended firewall:
**blind mechanical extraction only** (fixed column subsets → parquet; raw HDF5 deleted;
logs record structure (group names, columns, row counts), never sample values).
**No predictive score, summary statistic, or plot of any GWTC-4/5 posterior sample may
be produced until `PRE_REGISTRATION.md` is locked and committed.** Acquisition log:
`data/manifests/o4a_pe_acquisition_log.json`. GWTC-5.0 posteriors remain hard-blocked
in code (`--tier gwtc5-pe`) until after lock (and exceed the Stage-1 budget anyway).

## Records and what Stage 1 did with each

### (a) Training epoch: the frozen GWTC-3 fit · Zenodo [5655785](https://zenodo.org/records/5655785) (2021-11-05)

The data payload is **monolith-only**: `GWTC-3-population-data.tar.gz`, 10,070,439,649 B,
md5 `3c51561b5d8624210685b179c7d1f6ca`, Zenodo file id
`e5e046ef-be9a-4199-b8b2-1f0de620ad3c`. No per-analysis files exist on Zenodo, and the
paper's DCC page (P2100239) hosts only the PDF, re-checked 2026-06-10. A 10 GB download
breaks the budget, so `--tier gwtc3-plp` **streams** the tarball (HTTP → gunzip → tar,
never stored) and keeps only `analyses/PowerLawPeak/*`, the Power-Law+Peak + iid-spin +
power-law-redshift hyperposterior (`o1o2o3_mass_c_iid_mag_iid_tilt_powerlaw_redshift_result.json`)
and its PPD grids (`*_data.h5`), paths pinned from the release tutorial. Early-stop after
the PowerLawPeak directory ⇒ the full-tarball md5 cannot be re-verified on a partial
stream; instead each kept member is sha256-pinned in
`data/manifests/gwtc3_tarball_stream_manifest.json` (with the member list seen).
Small auxiliaries fetched whole: README, py_requirements.txt, table_data.tar.gz,
the tutorial HTML.

Note, confirmed again 2026-06-10: record **6513631 is the GWTC-2.1 PE release, the wrong
record**; the corrected training record is 5655785 (spec note stands).

### (b) Scoring targets: GWTC-4.0 O4a · Zenodo [16053484](https://zenodo.org/records/16053484) (2025-08-26) + [16911563](https://zenodo.org/records/16911563)

- **Event selection cut (the LVK GWTC-4.0 population paper's exact set):**
  `o4a_event_list.tar` (10,240 B, md5 `35ed4499b4abf8deccd53668a335bbf5`, file id
  `02bdee5e-a5d8-42eb-a762-78f6d752fed5`) → extracted to
  `gwtc4_pop/o4a_event_list/`: **84 BBH events** (`events_list_bbh_only.txt`, FAR < 1/yr)
  with the per-event PE label the LVK analysis used (**41 × `C00:Mixed`,
  43 × `C00:NRSur7dq4`**); full-spectrum list = 76 events.
- **Per-event posteriors:** 86 `IGWN-GWTC4p0-...-combined_PEDataRelease.hdf5` files
  (14.31 GB nominal, 53-575 MB each; per-file md5s + file ids in
  `zenodo_listings/16053484.json`). Acquired **one at a time**: download → extract the
  mass/spin/redshift/distance column subset (float32, zstd parquet) for **both** the
  LVK-used label and `C00:Mixed` → delete the HDF5. Columns wanted:
  `mass_1_source, mass_2_source, chirp_mass_source, total_mass_source, mass_ratio,
  redshift, luminosity_distance, comoving_distance, a_1, a_2, tilt_1, tilt_2,
  cos_tilt_1, cos_tilt_2, chi_eff, chi_p` (presence per event recorded in the log).
  Output: `data/derived/o4a_pe/<event>__<label>.parquet`.
  NOT fetched: `Archived_Skymaps.tar.gz` (268 MB, irrelevant), `PESummaryTable.hdf5`
  (204 KB but contains per-event outcome summaries, deliberately untouched until after
  lock). Fetched docs: release notebook + `md5sums.txt`.
- **GWOSC event lists** (free JSON, no auth): `GWTC-3-confident` (35), `GWTC-4.0` (129;
  **pinned**, GWTC-4.1 now exists with 140 but the population paper used 4.0),
  `GWTC-5.0` (161) → `data/gwosc/`.

### (c) Selection function: O4a sensitivity injections · Zenodo [16740117](https://zenodo.org/records/16740117) + [16740128](https://zenodo.org/records/16740128) (2025-08-04)

Both stored as a single compound dataset `events` (schema in the manifest). Extracted
fields: source masses/redshift/distance, all spin parametrisations, **all `lnpdraw_*`
draw-pdf components (float64)**, `weights`/`semianalytic_weights_*` (float64), all
per-pipeline `*_far` and `*_p_astro` (float32); root attrs (incl. `total_analysis_time`,
`total_generated`, `num_accepted`, `searches`) preserved in the manifest. Raw HDF
deleted after extraction.

| File | Size | md5 | Zenodo file id | Rows |
|---|---|---|---|---|
| `samples-rpo4a_v2_20250503133839UTC-1366933504-23846400.hdf` (O4a-only, rec 16740117) | 1,442,289,200 B | `1cf34f97...` | `8247cd41-9b7f-426f-8090-584374a19815` | 1,499,244 |
| `mixture-semi_o1_o2-real_o3-cartesian_spins_20250503134659UTC.hdf` (O3-era cumulative, rec 16740128) | 257,339,480 B | `e5d860f3...` | `5f995a1e-689e-4f3f-a95e-e7bc81f1f239` | 919,033 |

Output: `data/derived/selection/o4a_injections.parquet`,
`data/derived/selection/o1o2o3_mixture_cartesian.parquet` +
`data/manifests/selection_extraction_manifest.json`. The polar-spin mixture twin (same
size) is the documented alternative if the lock picks polar coordinates.

### (d) Comparison target: LVK GWTC-4.0 population release · Zenodo [16911563](https://zenodo.org/records/16911563)

Fetched now: file listing (`zenodo_listings/16911563.json`), `o4a_event_list.tar` (above),
`README.md`, `figure_scripts.tar` (113 KB). **DEFERRED (documented staged download):**
`analyses_BBH.tar`: 7,465,451,520 B, md5 `ba554251bacbda3979206ba673c644f5`, file id
`e9d4de9d-1643-4026-8c59-66793c02ddbd`. Plan: once the Stage-1 budget lifts (disk
reservation released), fetch via HTTP **range requests** (plain `.tar`, so member
boundaries are seekable): read the tar index headers, pull only the BBH PowerLawPeak /
Default popsummary HDF5 members, never store the full tar; fallback = full download →
selective untar → delete tar (needs ~15 GB transient, deferred to that window). We do NOT fetch
`analyses_BGP.tar` (23.5 GB) or `analyses_AllCBC.tar` for the headline result.

### GWTC-5.0 / O4b epoch (deferred to its own phase)

- **Sensitivity** Zenodo [19500064](https://zenodo.org/records/19500064):
  `samples-rpo4ab-1366933504-55469568-clipped.hdf`, 2,941,792,668 B, md5
  `a1106c27ec6cfd906231613523f7b174`, file id `567694de-3cbe-453c-9857-1f6b8cd07740`.
  Fits the budget alone but not alongside the Stage-1 working set → fetch with the same
  extract-and-delete discipline when the O4b epoch starts (expected parquet ≲ 0.7 GB).
  Release notes (`gwtc-5_o4ab_sensitivity-estimates.md`) already fetched.
- **Population release** Zenodo [20292639](https://zenodo.org/records/20292639): fetched
  `Event_list.tar.gz` (1,595 B → `GWTC5_BBH.txt` = **104 events**, `GWTC4.1_BBH.txt` =
  86) + README. `popsummary_files.tar.gz` (25.03 GB, gzipped = not range-seekable)
  is **out of scope** for the headline (GWTC-4 `analyses_BBH.tar` carries the
  comparison); revisit only if an O4b-epoch comparison figure demands it.
- **Per-event PE** Zenodo 20348005 (87 files, 46.45 GB) + 20348006 (21 files, 8.96 GB):
  **hard-blocked until pre-registration locks**; largest single file (6.59 GB,
  GW240615_113620) exceeds the Stage-1 budget even one-at-a-time → stream in a later
  disk window with ~7 GB transient scratch, same blind extract-and-delete loop as `o4a-pe`.

## Disk accounting

Targets: held new files ≤ 3.0 GB at every instant (checked in code before every
download/extraction), ≥ 1.2 GB free floor. Long-term keep after Stage 1: derived
parquets + extracted PowerLawPeak members + metadata (≲ 1.7 GB) + 0.25 GB venv.
Transients (each deleted before the next begins): one injection HDF (≤ 1.44 GB) or one
PE HDF5 (≤ 0.58 GB).

## Pinned toolchain

Fetch stack (installed, `.venv`, Python 3.14.2): `numpy 2.4.6 · h5py 3.16.0 ·
pandas 3.0.3 · pyarrow 24.0.0`. Science stack pinned at pre-registration lock
(candidates verified on PyPI 2026-06-10: `gwpopulation 1.3.1 · bilby 2.8.0 ·
pesummary 1.6.4 · popsummary 0.1.0`).
