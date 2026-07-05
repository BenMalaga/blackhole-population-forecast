# Paper outline: arXiv astro-ph.HE preprint

**Working title:** *Did We See It Coming? A Pre-Registered Out-of-Sample Test of
Black-Hole Population Models Against GWTC-4.0 and GWTC-5.0*

**Target:** arXiv astro-ph.HE (cross-list gr-qc). The Method section below is the
locked design of `PRE_REGISTRATION.md` (committed 2026-06-10, before any predictive
score existed); the paper executes it and does not reinterpret it. All Results entries
are placeholders: nothing has been unblinded.

---

## Abstract (skeleton)

Population models for merging binary black holes are routinely *fit* to each new
catalog and shown to *describe* it. We test instead whether the field's fiducial model
*predicts*: we freeze the LVK GWTC-3 Power-Law+Peak population fit (O3 data through
2020-03-27, released 2021) exactly as released, forward-model the detectable BBH
population of two fully held-out epochs (O4a/GWTC-4.0, N = 84; O4b/GWTC-5.0, N = 104)
through the corresponding public sensitivity injections, and score the actual catalogs
with statistics pre-registered before any score was computed. [RESULT-H1-EPOCH1: PIT
uniformity verdict for m1, z.] [RESULT-H2-EPOCH1: count coverage verdict.]
[RESULT-EPOCH2: replication verdict.] [HEADLINE SENTENCE: validation or falsification,
calibrated against finite-catalog variance.] The scoring framework is selection-matched
to the LVK analysis cuts, calibrated against catalog variance, and re-runnable on every
future catalog release.

## 1. Introduction

- The fit-describe-refit cycle in GW population astronomy; each catalog release
  triggers a wave of cumulative reanalyses claiming new structure (mass-spectrum
  features, spin subpopulations, redshift evolution, see `related_work.md`).
- The unasked prior question: would the pre-O4 model have predicted the O4 data? A
  genuine out-of-sample test is the standard of evidence for "the model extrapolates"
  vs "the model interpolates."
- Why this is non-trivial: selection effects must be forward-modeled per epoch;
  predictive checks are weak when applied in-sample (Miller et al. 2026); finite-catalog
  ("catalog") variance can mimic model failure (Corelli et al. 2026).
- Contributions: (i) a pre-registered two-epoch verdict on the GWTC-3 fiducial model,
  with either outcome informative; (ii) a reusable cross-catalog predictive-scoring
  framework; (iii) full reproducibility from public data products only.
- Pre-registration statement: design, statistics, thresholds, and event lists locked
  and git-committed before any predictive score existed; deviations, if any, labeled
  amendments.

## 2. Method (locked design: PRE_REGISTRATION.md)

### 2.1 Frozen model (training epoch)

- GWTC-3 population data release (Zenodo 5655785, 2021-11-05); Power-Law+Peak mass ×
  iid spin magnitude/tilt × power-law-in-(1+z) redshift evolution hyperposterior, used
  exactly as released (all 11,469 samples; no refit, no subsampling beyond seeded
  chunking). Released rate posterior used as-is.
- Provenance: tarball md5 + per-member sha256 pinned; forward-modeling via
  `gwpopulation 1.3.1` evaluated at the released hyperparameter samples; convention
  mappings resolved against the release's own documentation.

### 2.2 Held-out epochs and event inclusion

- Epoch 1 (PRIMARY): O4a / GWTC-4.0, FAR < 1/yr BBH, N = 84, the LVK population
  analysis's own byte-pinned event list (Zenodo 16911563). Epoch 2 (SECONDARY,
  replication): O4b / GWTC-5.0, N = 104 (Zenodo 20292639). GWTC-4.0 vs 4.1 list
  difference handled as a labeled robustness check.
- Per-event PE: the LVK per-event population-analysis sample labels exactly as recorded
  in the pinned lists (primary); single-label variant as robustness (§6.2).
- Locked exclusion rules; any technically dropped event named.

### 2.3 Selection function

- Per-epoch public sensitivity-injection sets (O4a: Zenodo 16740117; O4ab joint:
  Zenodo 19500064); detection = min-over-searches FAR < 1/yr, matching the event cut.
- VT estimator verbatim from the release documentation (importance reweighting of the
  injection mixture); detected-population predictive densities from the same weights.
- O4b restriction recipe for the joint injection set + pre-specified consistency check
  against the O4a-only set, with a locked fallback.
- Effective-sample-size diagnostic (Neff > 4 N_events) with a locked reporting rule.

### 2.4 Test statistics

- **H1 (distribution shape):** per-event posterior-averaged PIT values against the
  frozen detected-population predictive CDF, per observable (primary: m1, z; supporting:
  q; secondary: spin parameters); KS and AD statistics against uniformity; falsification
  iff the observed statistic exceeds the 97.5th percentile of the mock-null
  distribution.
- **H2 (detected count):** Poisson mixture over the frozen hyperposterior with
  injection-estimated VT per epoch; falsification iff the observed count lies outside
  the central 95% of the predictive distribution. Numeric predictive distributions
  committed before any H1 unblinding; the disclosure that the observed counts are
  public catalog metadata at lock is stated plainly (§6.5).
- **Per-event surprise index** (descriptive): highest-predictive-density tail
  probability of each event's (m1, z), posterior-averaged, with a
  catalog-variance-calibrated reference band; no accept/reject role.
- Data-level and split/partial predictive checks (after Miller et al. 2026) as labeled
  supporting diagnostics, never naive uncalibrated PPC p-values.

### 2.5 Null calibration and the mock gate

- All nulls from mock catalogs drawn from the frozen model itself, pushed through the
  locked injection-based selection and a measurement-noise emulation, including
  finite-catalog fluctuations per the bootstrap logic of Corelli et al. (2026). An
  apparent failure inside the catalog-variance null is not a falsification.
- Blocking mock-calibration gate: the pipeline must recover uniform PITs and count
  coverage on its own mocks at the locked rates before any real-data score; calibration
  report committed first.
- Positive control: reproduce Miller et al.'s known in-sample GWTC-4 spin-model failure
  signature before believing any out-of-sample null.

### 2.6 PE prior treatment

- Standard LVK PE prior conversion to the source-frame (m1, q, z) measure; per-event
  documented deviations recomputed from the stored priors group and named.

### 2.7 Reproducibility

- Pinned toolchain (Python 3.14.2; gwpopulation 1.3.1, bilby 2.8.0; full pins in
  `requirements.txt`); seed 20260610 for all stochastic steps; laptop-scale CPU.

## 3. Data

- Table: the five public data products (frozen fit, two event lists, two injection
  sets) with Zenodo records, file checksums, and access dates; per-event PE releases.
- Acquisition discipline: fixed-column blind extraction of posterior files, raw files
  deleted, structure-only logging; acquisition manifests committed.

## 4. Results: ALL PLACEHOLDERS (nothing unblinded)

- 4.1 Mock-calibration gate report. [PLACEHOLDER: gate pass/fail rates, committed
  before unblinding.]
- 4.2 Positive control (in-sample spin-model failure recovery). [PLACEHOLDER.]
- 4.3 Epoch 1 (O4a), H1: PIT histograms + KS/AD vs mock nulls for m1, z (headline);
  q; spin observables. [PLACEHOLDER: verdict per observable.]
- 4.4 Epoch 1 (O4a), H2: frozen predictive count distribution vs N = 84.
  [PLACEHOLDER: coverage verdict.]
- 4.5 Epoch 2 (O4b), H1 + H2, identical machinery. [PLACEHOLDER: replication verdict.]
- 4.6 Per-event surprise table and (m1, z) map with calibrated reference band.
  [PLACEHOLDER: most/least surprising events.]
- 4.7 Robustness: GWTC-4.1 event list; single-label PE variant; seed variations;
  O4b-restriction consistency check. [PLACEHOLDER.]

## 5. Discussion

- What a validation means: the O3-era recipe extrapolates one/two epochs forward;
  reported new features (e.g. the 35 solar-mass overdensity's evolution, spin
  structure) are refinements within the frozen model's predictive envelope, consistent
  with the catalog-variance reading of Corelli et al. [CONDITIONAL ON RESULTS.]
- What a falsification means: which observable failed, in which direction; comparison
  with the failure modes the in-sample reanalysis literature proposes. [CONDITIONAL.]
- Relation to sequential updating (Wolfe et al. 2026): change-of-posterior vs
  surprise-under-forecast as complementary diagnostics.
- The framework as a standing benchmark: re-running at each future catalog release
  turns population inference into an iteratively scored forecasting problem.

## 6. Limitations (pre-empting the obvious objections)

1. **Catalog selection effects.** The selection function is estimated from the LVK's
   public injection campaigns, with the detection rule (min-search FAR < 1/yr) matched
   to the event-inclusion cut. Residual mismatch between the injection campaigns and the
   true search sensitivity (waveform coverage of the injections, per-pipeline FAR
   calibration, the clipped/joint structure of the O4ab set) propagates directly into
   the predicted detected distributions and counts. The mock gate tests internal
   consistency of our pipeline, not the fidelity of the injection sets themselves; H2
   is the more exposed of the two hypotheses.
2. **Waveform and PE systematics.** O3-era and O4-era posteriors come from different
   waveform families and PE pipelines (e.g. NRSur7dq4-labeled vs mixed-sample releases).
   We never make per-event waveform choices: primary scoring uses exactly the per-event
   samples the LVK population analysis used, with a single-label variant as robustness.
   A systematic shift common to O4 PE would still register as model failure; the
   honesty pass (§8.7 of the pre-registration) requires attempting this attribution
   before any misspecification claim.
3. **Comparing population fits made with different pipelines.** The frozen GWTC-3
   hyperposterior was produced by the 2021 LVK analysis chain; our forward model
   re-implements the released model family via `gwpopulation` at the released
   hyperparameter samples. Convention mismatches (parameterization, spin-magnitude
   variance vs standard-deviation conventions, rate normalization) are resolved against
   the release's own documentation and the released PPD grids are used as cross-checks;
   any residual implementation difference is a stated systematic, distinct from the
   model's predictive performance. Likewise, descriptive comparisons to the LVK's own
   GWTC-4/5 hyperposteriors juxtapose fits produced by different pipelines on different
   cumulative datasets and are presented as context, never as the test.
4. **PE prior treatment.** The analytic standard-prior conversion is assumed unless an
   event's release documents otherwise; an undocumented non-standard sampling prior
   would bias that event's PIT.
5. **Known counts at lock.** The held-out event counts (84, 104) are public catalog
   metadata and were known when the pre-registration was locked; the H2 claim is
   therefore that a locked formula applied to a 2021-frozen model covers (or not) a
   known number, an out-of-sample test of the model, not a blinded guess. Stated in
   the pre-registration header and here.
6. **Statistical power.** With N = 84/104, KS-type tests resolve gross
   mis-extrapolation (CDF distortions of order 0.15), not few-percent refinements; the
   count test inherits the width of the frozen rate posterior. A null is therefore "no
   gross failure," not "the model is correct."
7. **Catalog variance.** Conversely, the calibrated nulls deliberately absorb
   finite-catalog fluctuations; a real but small population drift inside that envelope
   will not be flagged. This is the designed trade-off, not an oversight.
8. **One frozen model.** We test the LVK fiducial Power-Law+Peak chain, not every
   O3-era model; conclusions attach to that recipe (the field's de facto default), and
   the framework accepts any other frozen fit as input.

## 7. Conclusion

- One-paragraph restatement of the verdict [PLACEHOLDER] + the standing-benchmark
  framework + invitation to score future catalogs (GWTC-6 and beyond) against any
  frozen fit.

## Data and code availability

- All inputs are public LVK data products (Zenodo records pinned in `data/README.md`).
- Code: this repository (scoring framework, manifests, locked pre-registration);
  derived scoring tables and figures in `results/` on completion; archival DOI to be
  minted at release.

## References

- Maintained working bibliography with full citations: `docs/related_work.md`.

## Figure and table plan (provisional)

- Fig. 1: design schematic (frozen 2021 fit → per-epoch selection → predictive
  distributions → scored catalogs).
- Fig. 2: detected-population predictive density (m1, z) per epoch with event overlay.
  [Drawn only after unblinding.]
- Fig. 3: PIT uniformity panels + mock-null envelopes per observable per epoch.
- Fig. 4: H2 predictive count distributions vs observed counts.
- Fig. 5: per-event surprise map with calibrated band.
- Table 1: data products + checksums. Table 2: locked statistics and thresholds.
  Table 3: per-epoch verdicts. Table 4: robustness variants.
