# PRE-REGISTRATION: LOCKED 2026-06-10

> **This document is binding as of the git commit that adds it.** After lock, changes are
> amendments: committed separately, labeled `AMENDMENT`, dated, justified. Anything not
> specified at lock is exploratory and will be labeled as such in the writeup.
>
> **Status at lock, stated plainly:**
> - **No out-of-sample predictive score of any kind has been computed at lock time.** No
>   forward model exists; no scoring code exists (`src/` contains only the fetch/status
>   tooling); no statistic, plot, or summary of any GWTC-4/5 posterior sample value has
>   been produced or viewed. Scoring code must re-verify at runtime that this file is
>   committed before producing any score.
> - GWTC-4.0 per-event posteriors were partially acquired pre-lock (7 of 86 files) under
>   the blind-extraction firewall described in §10, structure-only logging, fixed column
>   subsets, raw files deleted. Disclosed in full in §10.
> - The observed event **counts** of both held-out epochs (84 and 104; §4) are public
>   catalog metadata, were used to pin the comparison sets, and are therefore **known at
>   lock**. The count test (H2) does not pretend otherwise: its integrity rests on the
>   frozen 2021 model being used as released and the count formula being locked here,
>   verbatim from the sensitivity release's own documentation, before it is ever
>   evaluated (§6.2).

# Did We See It Coming? Out-of-Sample Test of Black-Hole Population Models

**Design (two-epoch, single frozen fit):** Freeze the LVK GWTC-3 population fit (O3 data
through 2020-03-27; release 2021-11-05). Forward-model the *detected* BBH population of
two later, fully held-out epochs through the public sensitivity injections, and score the
actual catalogs against that frozen prediction:

- **Epoch 1, O4a / GWTC-4.0 (PRIMARY).** The headline result.
- **Epoch 2, O4b / GWTC-5.0 (SECONDARY).** Same frozen fit, identical statistics and
  thresholds, reported alongside. Declared secondary at lock because the GWTC-5 data
  release (2026-05-26) post-dates the project spec and its PE acquisition is deferred
  (§10); its role is replication, not headline.

All dataset facts below were verified 2026-06-10 against live primary sources (Zenodo
REST API, GWOSC event API, PyPI, arXiv). This document supersedes
`PRE_REGISTRATION_DRAFT.md`, deleted in the lock commit (its full text remains in git
history; the firewall disclosure it introduced is carried into §10).

---

## 1. Research question and hypotheses

Population "recipes" for merging black holes are routinely *fit* to a catalog and shown
to *describe* it. We test whether the field's flagship model, frozen at GWTC-3,
**predicts** the next two catalogs out-of-sample. This is the out-of-sample complement to
Miller et al. (arXiv:2604.06090), whose predictive checks are in-sample.

### H1: distribution shape (primary hypothesis)

The Power-Law+Peak mass model with iid spin magnitude/tilt and power-law redshift
evolution, using the **GWTC-3 hyperposterior only** (§2), correctly predicts the
observable distributions of the detected BBH events of each held-out epoch:

- **Primary observables:** source-frame primary mass `m1` and redshift `z` (the spec's
  named failure modes: peak shifts, new high-mass structure). Mass ratio `q` is reported
  with the same machinery, labeled supporting.
- **Secondary observables (labeled):** spins: `chi_eff`, spin magnitudes `a_1, a_2`,
  tilts `cos_tilt_1, cos_tilt_2`. Miller et al.'s GWTC-4 finding (the Default spin model
  under-predicts high spin magnitudes and over-predicts anti-aligned tilts) is the
  **positive control**: the pipeline must be able to reproduce a known in-sample failure
  signature before its out-of-sample nulls are believed.

**Test statistics and decision rule (unchanged from the project spec / draft):**

1. For each detected event `i` and observable `u`, the **posterior-averaged PIT**
   `PIT_i = (1/K_i) * sum_k F_pred(u_ik)` where `{u_ik}` are the event's PE samples and
   `F_pred` is the CDF of the frozen detected-population predictive distribution for the
   epoch (§6.1). The identical operator is applied to mock events (with measurement-noise
   emulation) when building nulls, so the null distribution of every statistic is exact
   by construction regardless of this operator choice.
2. Aggregate per observable per epoch: **Kolmogorov–Smirnov distance** and
   **Anderson–Darling statistic** of `{PIT_i}` against uniformity.
3. **Null calibration (mandatory gate, §7):** null distributions of every statistic come
   from mock catalogs drawn from the frozen model itself, pushed through the same
   injection-based selection and measurement-noise emulation, **including finite-catalog
   ("catalog variance") fluctuations** per the bootstrap logic of arXiv:2603.00239. An
   apparent failure that sits inside the catalog-variance null is **not** a
   falsification.
4. **H1 is falsified (per epoch, per observable)** iff the observed statistic exceeds
   the **97.5th percentile of the mock-null distribution**. The headline claim uses `m1`
   and `z` in Epoch 1.
5. **Data-level and split/partial predictive checks** in the style of Miller et al. are
   run as labeled supporting diagnostics with the same null-calibration principle (naive
   uncalibrated PPC p-values are not used anywhere, spec reviewer fix).
6. **Per-event surprise index (descriptive, pre-defined):** for each event,
   `s_i = (1/K_i) * sum_k Pr_{x~pred}[ f_pred(x) <= f_pred(x_ik) ]`, the
   highest-density tail probability of the event's `(m1, z)` under the frozen detected-
   population density, averaged over the event's posterior samples, reported with a
   catalog-variance-calibrated reference band. No accept/reject threshold; it feeds the
   per-event table and figures only.

### H2: detected event count (sharper headline)

The frozen GWTC-3 rate + population model predicts the **number of detected FAR < 1/yr
BBH events** per epoch.

- **Predictive distribution:** `P(N | epoch) = (1/S) * sum_s Poisson(N; mu(Lambda_s))`
  over the S frozen hyperposterior samples, with `mu(Lambda_s) = R(Lambda_s) *
  VT_epoch(Lambda_s)` where `VT_epoch` is the injection-estimated sensitive time–volume
  computed **exactly** by the sensitivity release's documented estimator (§5, formula
  locked verbatim from the release notes).
- **H2 is falsified (per epoch)** iff the observed count (Epoch 1: **84**; Epoch 2:
  **104**; §4) lies outside the **central 95%** of the predictive distribution
  (2.5th–97.5th percentiles).
- The numeric predictive distributions are computed by the locked formula and committed
  to `results/` before any H1 unblinding; as disclosed in the header, the observed
  counts are already public and known, so the H2 claim is "the frozen model + locked
  formula yields a count interval that does/does not cover the known N", an
  out-of-sample test of the *model*, not a blinded guess. There is no analyst freedom
  left between this document and the number.

### Null results are wins

"The 2021 model predicts the 2025–26 catalogs" is a publishable validation; "it does
not" is a publishable falsification, provided either claim survives the
catalog-variance calibration. The durable artifact is the **cross-catalog
predictive-scoring framework**, re-runnable on every future catalog.

---

## 2. Frozen model (training epoch: input data, not outcome data)

- **Source:** GWTC-3 population inference data release, Zenodo **5655785** (published
  2021-11-05, CC-BY-4.0). Payload `GWTC-3-population-data.tar.gz`,
  **10,070,439,649 bytes, md5 `3c51561b5d8624210685b179c7d1f6ca`** (Zenodo REST API,
  verified 2026-06-10). This md5 cryptographically fixes the identity of every member
  file. Record 6513631 is the GWTC-2.1 PE release, wrong record; the spec's
  correction stands (re-verified 2026-06-10).
- **Primary frozen artifact (exact member paths, pinned from the release tutorial):**
  - `analyses/PowerLawPeak/o1o2o3_mass_c_iid_mag_iid_tilt_powerlaw_redshift_result.json`,
    the Power-Law+Peak + iid-spin + power-law-redshift hyperposterior. **All** released
    hyperposterior samples are used; no subsampling except seeded computational chunking
    (§9).
  - `analyses/PowerLawPeak/*_data.h5`, the released PPD grids (4 files), used for
    cross-checks of our forward model only, never as the prediction itself.
- **Member checksums:** the streaming extraction (`--tier gwtc3-plp`, in flight at lock;
  raw tarball never stored) records each kept member's sha256 in
  `data/manifests/gwtc3_tarball_stream_manifest.json`. Because the member content is
  already fully determined by the pinned tarball md5 + member paths above, these sha256s
  add provenance, not freedom: when the stream completes they are committed and, if not
  already in this file's lock commit, appended by a labeled **non-substantive provenance
  amendment**.
- **Rate normalization (locked rule, file-determined):** if the released result file
  contains the merger-rate posterior, it is used as-is; if the rate must be
  reconstructed, the reconstruction follows the release's own README/tutorial recipe,
  documented in a labeled amendment **before** any H2 evaluation. No third option.
- **Optional robustness epoch (exploratory):** the GWTC-2-era fit (LIGO DCC
  P2000434-v3). Not part of H1/H2.

## 3. Population model definition

The frozen prediction uses the hyperposterior **as released**, we do not refit.
Forward-modeling uses `gwpopulation`'s implementations of the same models the release
documents for this analysis chain (Power-Law+Peak mass; iid spin magnitude Beta
distribution; iid spin tilt isotropic+aligned mixture; power-law-in-(1+z) redshift
evolution), evaluated at the released hyperparameter samples. Any naming/parameter
convention mismatch between the release file and `gwpopulation 1.3.1` is resolved
against the release README/tutorial, documented in the repo, and is a correctness
question, not a tunable choice.

## 4. Held-out epochs and event-inclusion criteria (locked)

**Inclusion rule (both epochs, identical):** FAR < 1/yr BBH events **as pinned by the
LVK's own released population-analysis event lists**, the comparison sets are byte-pinned
files, not our re-derivation:

| Epoch | Catalog pin | Event list (byte-pinned) | N |
|---|---|---|---|
| 1 (O4a, PRIMARY) | GWTC-**4.0** (GWOSC `GWTC-4.0`, 129 events; **not** 4.1) | `events_list_bbh_only.txt` from `o4a_event_list.tar` (10,240 B, md5 `35ed4499b4abf8deccd53668a335bbf5`, Zenodo 16911563); sha256 `d1edcc80447b628268902ede91ad047c070910c8fb41a3051fa55ed34861b31f` | **84** |
| 2 (O4b, SECONDARY) | GWTC-5.0 (GWOSC `GWTC-5.0`, 161 events) | `GWTC5_BBH.txt` from `Event_list.tar.gz` (1,595 B, md5 `a4d3a105…`, Zenodo 20292639); sha256 `cdf3fd8f30ed842bbe9d2219616366da92eaa207fa83277fb2dfc77be5e9eb5a` | **104** |

Full event-name lists are reproduced verbatim in **Appendix A** (locked).

- **Catalog-version pin:** GWOSC also serves GWTC-4.1 (140 events); the LVK *GWTC-5*
  population paper uses a GWTC-4.**1** O4a list (86 BBH). Our Epoch 1 stays pinned to
  the GWTC-4.0 84-event list so the comparison set is identical to the LVK GWTC-4.0
  population paper's. The 4.0→4.1 delta is a labeled robustness check, not the primary
  test.
- **Posterior samples:** per-event combined PE release files
  (`IGWN-GWTC4p0-*-combined_PEDataRelease.hdf5`, Zenodo 16053484; GWTC-5 analogues in
  Zenodo 20348005 + 20348006). **PE analysis-label convention (locked):** primary = the
  **LVK per-event population-analysis label** exactly as recorded in the pinned event
  list (41 × `C00:Mixed`, 43 × `C00:NRSur7dq4` for Epoch 1), matching the LVK population
  paper's inputs; `C00:Mixed`-everywhere is the labeled robustness variant. Both were
  extracted blind for every acquired event, so the choice required no peeking and no
  re-download.

### Exclusions (locked, applied identically to both epochs)

1. Events with FAR ≥ 1/yr per the pinned lists, excluded.
2. Events the LVK population analysis classifies as non-BBH, excluded (the pinned lists
   carry the LVK classification; e.g. GW230529_181500 is outside the BBH cut).
3. Events lacking the required PE release file/label, excluded and named in the writeup
   (none expected; verified at fetch time).
4. GWTC-4.1-only events, excluded from the primary test (robustness only).
5. No other exclusions. Any event dropped for a technical reason (corrupt file, failed
   extraction) is named, and the analysis is re-run including it where possible.

## 5. Selection function (locked recipe, verbatim from the release documentation)

- **Epoch 1:** O4a injection set
  `samples-rpo4a_v2_20250503133839UTC-1366933504-23846400.hdf` (1,442,289,200 B, md5
  `1cf34f97…`), Zenodo **16740117**.
- **Epoch 2:** joint O4a+O4b injection set
  `samples-rpo4ab-1366933504-55469568-clipped.hdf` (2,941,792,668 B, md5 `a1106c27…`),
  Zenodo **19500064**.
- **Detection definition (locked, from the release notes' own example code):** an
  injection is detected iff `min_over_searches(<search>_far) < 1/yr`, the same FAR
  threshold as the event-inclusion cut, using all per-search `*_far` fields the release
  provides (absent detections default to `inf` per the release format).
- **VT estimator (locked, verbatim from the release notes):**
  `VT(Lambda) = total_analysis_time * sum_{sel} [ weights * exp(lnp(Lambda) - lnpdraw) ] / total_generated`
  with `lnpdraw` the release's summed draw-density fields, `weights` the release's
  mixture weights, and `total_analysis_time`, `total_generated` taken from the file's
  root attributes. Detected-population predictive densities for H1 are the
  correspondingly importance-reweighted injection distributions (normalization cancels).
- **Epoch-2 restriction recipe (locked decision tree):** the O4ab release states O4b =
  GPS `[1396969218, 1422118818]` (15:00 UTC 2024-04-10 → 17:00 UTC 2025-01-28).
  - *Primary:* restrict the joint set to injections with `time_geocenter` in the O4b
    window; the month-wise mixture `weights` partition the joint VT sum by construction,
    so the restricted sum is the O4b share.
  - *Mandatory consistency check (pre-unblinding):* the same construction on the O4a
    window of the joint set must reproduce the O4a-only set's `VT(Lambda)` within 10%
    (median over hyperposterior samples).
  - *Locked fallback if the check fails:* Epoch 2's count test becomes the **joint
    O4a+O4b count** (observed N = 188) using the joint set exactly as released (no
    restriction), and Epoch-2 shape tests use the window-restricted predictive density
    (valid regardless, since shape normalization cancels). The check outcome is reported
    either way.
- **Effective-sample-size diagnostic (locked):** importance-sampling
  `Neff > 4 * N_events` per epoch (the standard rule used in LVK population analyses /
  `gwpopulation`); hyperposterior samples failing it are flagged, and if > 1% of samples
  fail, the failure is reported prominently and the affected statistic labeled unstable.

## 6. PE prior treatment (locked)

1. **PE sampling prior:** assumed to be the LVK standard PE prior (uniform in
   detector-frame component masses, luminosity-distance prior ∝ d_L², isotropic spin
   orientations, uniform spin magnitudes), converted analytically to the source-frame
   (m1, q, z) measure, the same convention the GWTC-3 population analysis and
   `gwpopulation` use. PIT statistics (§1) compare *posterior samples* to the predictive
   distribution; the prior enters the data-level/split checks and the surprise-index
   density and is handled with this fixed convention.
2. **Per-event deviations:** if a release file's stored `priors` group documents a
   non-standard sampling prior for an event, the prior factor for that event is
   recomputed from the stored priors group (re-downloadable from the pinned Zenodo
   files), and the event is named in the writeup.
3. **Waveform systematics:** handled by the locked label convention (§4): primary scoring
   uses the exact per-event samples the LVK population analysis used; `C00:Mixed` is the
   robustness variant. No per-event waveform choice is made by us, ever.

## 7. Mock-calibration gate (blocking; unchanged from spec milestone W3–4)

Before any real-data score: mock catalogs drawn from the frozen model, pushed through
the locked injection selection (§5) and a measurement-noise emulation, must recover
uniform PITs and count coverage (the pipeline scoring its own mocks must pass its own
thresholds at the locked rates). Null envelopes include catalog-variance bootstrap per
arXiv:2603.00239. The calibration report is committed (`results/mock_calibration/`)
**before** any outcome posterior is scored. The measurement-noise emulation method is
internal to null construction and validated solely by this gate; it is documented in the
report. **If the gate fails, the main run does not start** until the pipeline is fixed
and the gate re-passed; any change to locked statistics this forces is a labeled
amendment, permitted only because it precedes all unblinding.

## 8. Analysis sequence (locked order)

1. Complete training-epoch + selection acquisition (in flight at lock; §10). Pin member
   sha256s (§2). *(input data, no firewall implications)*
2. Build the forward model (frozen hyperposterior → detected-population predictive
   density + count distribution per epoch). Commit the H2 numeric predictive
   distributions to `results/` at this step, before step 4.
3. **Mock-calibration gate (§7).** Blocking.
4. **Unblind Epoch 1 (O4a):** finish blind PE acquisition (GPS-time order, §10), then
   score H1/H2 with the locked statistics. First moment any predictive score exists.
5. **Unblind Epoch 2 (O4b):** acquire GWTC-5 PE (`--tier gwtc5-pe`, hard-blocked until
   after lock; ~7 GB transient scratch needed, later disk window), score
   identically.
6. Comparison figures against the LVK's own released GWTC-4.0 hyperposteriors
   (`analyses_BBH.tar`, Zenodo 16911563, staged range-request plan), descriptive
   context, not a locked test.
7. Honesty pass: adversarial attempt to attribute any apparent failure to catalog
   variance, PE-prior mismatch, or selection-model error **before** claiming model
   misspecification; nulls reported prominently.

## 9. Toolchain, seed, reproducibility (locked)

- Python 3.14.2 (`.venv`); installed and pinned: `numpy 2.4.6 · h5py 3.16.0 ·
  pandas 3.0.3 · pyarrow 24.0.0`. Science stack locked to: `gwpopulation 1.3.1 ·
  bilby 2.8.0 · pesummary 1.6.4 · popsummary 0.1.0` (verified on PyPI 2026-06-10;
  installed at step-2 start). Remaining utility pins (scipy/astropy/matplotlib/tqdm) are
  recorded in `requirements.txt` at install time, provenance, not statistical freedom.
- **Seed `20260610`** for all mock-catalog draws, bootstrap resampling, and any
  computational subsampling. Single-seed results are the headline; seed-robustness
  (3 extra seeds) is exploratory.
- All compute CPU-only, laptop-scale, $0 spend; every long computation chunked to fit
  8 GB RAM.

## 10. Data-acquisition state at lock (full disclosure)

- **Firewall (carried from the draft, where it was a pre-lock revision):**
  posterior/hyperposterior files are input data and may be acquired pre-lock, but only
  via blind mechanical extraction (`src/fetch_data.py`): fixed 16-column subsets →
  float32 parquet, raw HDF5 deleted, structure-only logging (group names, columns, row
  counts, never sample values), UTC-stamped in
  `data/manifests/o4a_pe_acquisition_log.json` (committed). **No statistic, plot, or
  summary of any GWTC-4/5 sample value was produced or viewed pre-lock.** GWTC-5 PE
  remains hard-blocked in code until post-lock.
- **Honesty note, restated:** with public LVK releases, *content* blindness was never
  available, the GWTC-4/5 papers are published and read (`docs/related_work.md`). The
  out-of-sample claim rests on (a) the frozen training artifact being LVK's own 2021
  release, (b) scoring rules committed here before any score exists, and (c) full
  reproducibility from public data, not on pretending we cannot know what O4a looks
  like.
- **Acquired at lock:** all metadata/event lists/release notes; O3-era injection mixture
  parquet (919,033 × 25; raw 257 MB deleted); **7 of 86** O4a PE files blind-extracted
  (11 parquet subsets incl. both label variants, 279,894 posterior samples held,
  counted from file metadata only). **Subset rule (pre-registered at acquisition):**
  events are acquired in **GPS-time order**; the 7/86 partial state is a download-budget
  artifact of a disk-constrained machine, not a selection on event properties.
  Acquisition resumes in the same order until 86/86.
- **In flight at lock:** the `gwtc3-plp` tarball stream has reached
  `analyses/PowerLawPeak/` and is extracting members (11 kept at lock commit; the
  primary `*_result.json` not yet reached; floor-waiting on concurrent disk
  usage), so `data/manifests/gwtc3_tarball_stream_manifest.json` does not exist yet, so
  the §2 provenance-amendment provision applies. The O4a injection-set download is
  likewise floor-waiting. Both are input-data acquisitions with no bearing on the
  firewall.
- **Spot-checks performed pre-lock (firewall-compliant):** value-range sanity on the
  *selection mixture* parquet only (input data: masses 1.0–999 M☉, z ≤ 1.90, |spin
  components| < 1, weights positive, no NaNs; noted sentinel −1.0 in
  `o3_pycbc_bbh_p_astro`, irrelevant under the FAR-only detection rule); **structure-only**
  check on one O4a PE parquet (16 float32 columns as specified, 41,712 rows).
- **Spin-coordinate convention (locked):** cartesian spin components for the O3-era
  mixture file (as extracted); the polar-coordinate twin is the documented alternative
  and is *not* used. O4a/O4ab injection sets provide both; the locked population models
  (§3) consume magnitude/tilt, computed deterministically from cartesian components.

## 11. Power / feasibility sketch (non-binding context, not a locked claim)

- N = 84 (Epoch 1), 104 (Epoch 2), 188 combined, vs ~70 BBH in the entire training
  catalog. A KS uniformity test at α = 0.05 with n ≈ 84 has ~80% power against CDF
  distortions D ≈ 0.15: the test resolves gross mis-extrapolation (shifted/absent
  ~35 M☉ peak, missing high-mass tail), not few-percent refinements, consistent with
  the project's question and honestly below in-sample-refit sensitivity.
- The count test inherits the GWTC-3 rate posterior's width (read from the frozen
  release at forward-model time, deliberately **not** quoted from memory here); if that
  width is of order a factor ~2, H2 falsifies only order-50%+ count mispredictions.
  Stated up front in the writeup.

## 12. Prior art engaged (reviewed at lock)

See `docs/related_work.md` (maintained). Anchors: Miller et al. arXiv:2604.06090
(in-sample PPC limitations, the out-of-sample complement); arXiv:2603.00239
(catalog variance, built into the nulls); LVK GWTC-4.0 population paper 2508.18083;
GWTC-5.0 papers 2605.27223–27227; the GWTC-4 rapid-response in-sample wave (10 papers
listed in related_work.md). **Prior-art review run 2026-06-10 at lock** (arXiv API:
astro-ph.HE/gr-qc, GWTC + out-of-sample / posterior-predictive / forecast / held-out /
cross-validation / frozen; all entries since 2026-06-01 reviewed): GWTC-5 has only
the five LVK collaboration papers (2026-05-26); the only new GWTC-mentioning astro-ph.HE
submissions since 06-01 are arXiv:2606.02318 (in-sample, already tracked), 2606.03346,
2606.04810, 2606.08838 (all irrelevant to predictive scoring). As of lock, no frozen-fit
predict-then-score test against GWTC-4 or GWTC-5 exists in the literature.

## 13. Known risks (declared at lock)

1. **Concurrent-publication timeline:** GWTC-5.0 was released 2026-05-26; a
   rapid-response reanalysis wave like GWTC-4's (5+ papers in 2–12 weeks) is expected.
   The pilot gate is budgeted in weeks, not months.
2. **PE-prior/waveform mismatch** is the main technical validity threat; §6 is the
   locked handling; the honesty pass (§8.7) must consider it before any misspecification
   claim.
3. **Catalog variance** can mimic model failure; §7's nulls are the defense, and any
   "failure" claim must survive them.
4. **Selection-effects expertise wall:** §7 is the tripwire; if mocks don't calibrate,
   the main run does not start.
5. **Environment:** disk-constrained shared machine; all acquisition tooling is
   idempotent and budget-guarded; partial states are committed and disclosed, never
   silent.

---

## Appendix A: locked event lists (verbatim from the byte-pinned LVK files)

### Epoch 1: O4a / GWTC-4.0 BBH FAR < 1/yr (N = 84), with the LVK per-event PE label

```
GW230601_224134,C00:NRSur7dq4    GW230605_065343,C00:Mixed
GW230606_004305,C00:NRSur7dq4    GW230608_205047,C00:NRSur7dq4
GW230609_064958,C00:NRSur7dq4    GW230624_113103,C00:Mixed
GW230627_015337,C00:Mixed        GW230628_231200,C00:NRSur7dq4
GW230630_125806,C00:NRSur7dq4    GW230630_234532,C00:Mixed
GW230702_185453,C00:NRSur7dq4    GW230704_021211,C00:NRSur7dq4
GW230704_212616,C00:Mixed        GW230706_104333,C00:Mixed
GW230707_124047,C00:NRSur7dq4    GW230708_053705,C00:NRSur7dq4
GW230708_230935,C00:NRSur7dq4    GW230709_122727,C00:NRSur7dq4
GW230712_090405,C00:Mixed        GW230723_101834,C00:Mixed
GW230726_002940,C00:NRSur7dq4    GW230729_082317,C00:Mixed
GW230731_215307,C00:Mixed        GW230803_033412,C00:Mixed
GW230805_034249,C00:NRSur7dq4    GW230806_204041,C00:NRSur7dq4
GW230811_032116,C00:Mixed        GW230814_061920,C00:NRSur7dq4
GW230814_230901,C00:NRSur7dq4    GW230819_171910,C00:Mixed
GW230820_212515,C00:Mixed        GW230824_033047,C00:NRSur7dq4
GW230825_041334,C00:NRSur7dq4    GW230831_015414,C00:Mixed
GW230904_051013,C00:Mixed        GW230911_195324,C00:Mixed
GW230914_111401,C00:NRSur7dq4    GW230919_215712,C00:Mixed
GW230920_071124,C00:NRSur7dq4    GW230922_020344,C00:NRSur7dq4
GW230922_040658,C00:NRSur7dq4    GW230924_124453,C00:NRSur7dq4
GW230927_043729,C00:NRSur7dq4    GW230927_153832,C00:Mixed
GW230928_215827,C00:NRSur7dq4    GW230930_110730,C00:NRSur7dq4
GW231001_140220,C00:Mixed        GW231004_232346,C00:Mixed
GW231005_021030,C00:Mixed        GW231005_091549,C00:NRSur7dq4
GW231008_142521,C00:NRSur7dq4    GW231014_040532,C00:Mixed
GW231018_233037,C00:Mixed        GW231020_142947,C00:Mixed
GW231028_153006,C00:NRSur7dq4    GW231029_111508,C00:NRSur7dq4
GW231102_071736,C00:NRSur7dq4    GW231104_133418,C00:Mixed
GW231108_125142,C00:Mixed        GW231110_040320,C00:Mixed
GW231113_122623,C00:Mixed        GW231113_200417,C00:Mixed
GW231114_043211,C00:Mixed        GW231118_005626,C00:Mixed
GW231118_071402,C00:NRSur7dq4    GW231118_090602,C00:Mixed
GW231119_075248,C00:NRSur7dq4    GW231123_135430,C00:NRSur7dq4
GW231127_165300,C00:NRSur7dq4    GW231129_081745,C00:NRSur7dq4
GW231206_233134,C00:NRSur7dq4    GW231206_233901,C00:NRSur7dq4
GW231213_111417,C00:NRSur7dq4    GW231221_135041,C00:Mixed
GW231223_032836,C00:NRSur7dq4    GW231223_075055,C00:Mixed
GW231223_202619,C00:Mixed        GW231224_024321,C00:Mixed
GW231226_101520,C00:NRSur7dq4    GW231230_170116,C00:Mixed
GW231231_154016,C00:Mixed        GW240104_164932,C00:NRSur7dq4
GW240107_013215,C00:Mixed        GW240109_050431,C00:Mixed
```

### Epoch 2: O4b / GWTC-5.0 BBH FAR < 1/yr (N = 104)

```
GW240413_022019  GW240414_054515  GW240420_175625  GW240426_031451
GW240428_225440  GW240501_033534  GW240505_133552  GW240507_041632
GW240511_031507  GW240512_024139  GW240513_183302  GW240514_121713
GW240515_005301  GW240519_012815  GW240520_213616  GW240525_031210
GW240526_093944  GW240527_183429  GW240527_230910  GW240530_012417
GW240531_040326  GW240531_075248  GW240601_061200  GW240601_231004
GW240612_081540  GW240615_113620  GW240615_160735  GW240618_071627
GW240621_195059  GW240621_200935  GW240621_214041  GW240622_004008
GW240627_131622  GW240629_145256  GW240630_101703  GW240703_191355
GW240705_053215  GW240716_034900  GW240824_205609  GW240825_055146
GW240830_211120  GW240902_143306  GW240907_153833  GW240908_082628
GW240908_125134  GW240910_103535  GW240915_001357  GW240915_105151
GW240916_184352  GW240919_061559  GW240920_073424  GW240920_124024
GW240921_201835  GW240922_142106  GW240923_204006  GW240924_000316
GW240925_005809  GW240930_035959  GW240930_234614  GW241002_030559
GW241006_015333  GW241007_082943  GW241009_022835  GW241009_084816
GW241009_220455  GW241011_233834  GW241101_220523  GW241102_124058
GW241102_144729  GW241109_033317  GW241109_115924  GW241110_124123
GW241111_111552  GW241113_163507  GW241114_024711  GW241114_235258
GW241116_151753  GW241124_024914  GW241125_010116  GW241127_061008
GW241129_021832  GW241130_034908  GW241130_110422  GW241201_055758
GW241210_060606  GW241210_120900  GW241225_042553  GW241225_082815
GW241229_155844  GW241230_084504  GW241230_233618  GW241231_054133
GW250101_011205  GW250104_015122  GW250108_152221  GW250109_010541
GW250109_074552  GW250114_082203  GW250116_015318  GW250118_023225
GW250118_055802  GW250118_170523  GW250119_025138  GW250119_190238
```

*(Authority order: if any transcription error is ever found in this appendix, the
byte-pinned source files (sha256s in §4) govern; such a correction is a labeled
non-substantive amendment.)*
