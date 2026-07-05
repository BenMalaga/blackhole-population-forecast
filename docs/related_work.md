# Related Work

*Literature review current as of 2026-06-14. All arXiv identifiers, titles, author
lists, and submission dates below were verified against arXiv on 2026-06-10/11, and
again on 2026-06-14. As of 2026-06-14 the literature contains no out-of-sample /
frozen-fit predictive-scoring work of the kind described here.*

This project asks a question that is logically prior to most current work on
gravitational-wave populations: does the binary-black-hole (BBH) population model the
LIGO–Virgo–KAGRA (LVK) collaboration fit to O3-era data *predict*, not merely describe,
the catalogs that came after it? We freeze the GWTC-3 Power-Law+Peak population fit
(released 2021), forward-model the detectable BBH population of two fully held-out
epochs (O4a/GWTC-4.0 and O4b/GWTC-5.0) through the public sensitivity injections, and
score the actual catalogs against that frozen prediction with statistics that were
pre-registered before any score existed. The relevant prior literature falls into four
strands: the collaboration's own catalogs and population analyses, the methodology of
predictive model checking, the statistics of finite catalogs, and the rapidly growing
body of in-sample reanalyses of the O4-era data.

## Collaboration catalogs and population analyses

The frozen model under test is the fiducial analysis of the GWTC-3 population paper
(LVK Collaboration, *The population of merging compact binaries inferred using
gravitational waves through GWTC-3*, arXiv:2111.03634; Phys. Rev. X 13, 011048): the
Power-Law+Peak mass function with iid spin magnitudes and tilts and power-law redshift
evolution, fit to events observed through 2020-03-27. The two held-out epochs are
defined by the GWTC-4.0 population analysis (LVK Collaboration, *GWTC-4.0: Population
Properties of Merging Compact Binaries*, arXiv:2508.18083, 2025-08-25), which analyzes
84 BBH events from O4a and reports a mass spectrum described by a broken power law with
overdensities near 10 and 35 solar masses, and by the GWTC-5.0 catalog and population
papers (LVK Collaboration, arXiv:2605.27225 and *GWTC-5.0: Population Properties of
Merging Compact Binaries*, arXiv:2605.27226, both 2026-05-26), which extend the
cumulative analysis to 267 mergers through O4b and report, among other findings,
evidence for a rapidly spinning subpopulation consistent with hierarchical mergers.
These are cumulative, in-sample fits by construction: each new catalog is folded into
the likelihood and the hyperposterior is re-inferred. Our analysis is deliberately
different in kind: the GWTC-3 hyperposterior is used exactly as released in 2021, and
the newer catalogs enter only as scoring data. We adopt the LVK's own event-inclusion
cuts (FAR < 1/yr BBH, byte-pinned event lists from the population-analysis data
releases) and sensitivity-injection products so that any disagreement we find is
attributable to the model, not to a mismatched comparison set.

## Predictive model checking for gravitational-wave populations

Posterior predictive checking in gravitational-wave astronomy was introduced
pedagogically by Romero-Shaw, Thrane & Lasky (*When models fail: an introduction to
posterior predictive checks and model misspecification in gravitational-wave
astronomy*, arXiv:2202.05479; PASA 39, e025, 2022). The methodological anchor for the
present work is Miller, Winney, Chatziioannou & Meyers (*Posterior Predictive Checks
for Gravitational-wave Populations: Limitations and Improvements*, arXiv:2604.06090,
2026-04-07), who systematically dissect how standard event-level posterior predictive
checks (PPCs) lose power for poorly measured single-event parameters: when single-event
inference is prior-dominated, traditional PPCs can fail to flag even badly misspecified
models. They evaluate alternatives (split and partial predictive checks, checks on
maximum-likelihood parameters) and, applying them to GWTC-4.0, find that the Gaussian
component-spin model under-predicts BBHs with large spin magnitudes and over-predicts
perfectly anti-aligned tilts. Two points of contact are central here. First, all of
Miller et al.'s checks are in-sample: the same catalog is used both to fit the
hyperposterior and to test it, which is exactly the regime in which predictive checks
are weakest and most conservative. Our design is the out-of-sample complement, the
scoring data never enter the fit, and so sidesteps the double-use-of-data problem
their paper characterizes. Second, their concrete GWTC-4.0 spin finding serves as a
positive control for our pipeline: a scoring framework that cannot reproduce a known
in-sample failure signature should not be trusted when it reports out-of-sample nulls.
Our pre-registration adopts both the lesson (no naive uncalibrated PPC p-values
anywhere) and the benchmark.

## Catalog variance: how much disagreement is just a finite catalog?

Corelli, Gerosa, Mould & Fabbri (*Variance of gravitational-wave populations*,
arXiv:2603.00239, 2026-02-27) quantify "catalog variance", the variability of inferred
population distributions that arises because the observed catalog is one finite
realization, by bootstrapping over both detected events and sensitivity injections.
Applied to GWTC-4, their analysis substantially broadens the effective uncertainties of
standard single-catalog inference; notably, the ~35 solar-mass peak in the primary-mass
distribution is largely absorbed by these statistical fluctuations. This result is a
direct warning for any out-of-sample comparison: a frozen model can "fail" against a
new catalog merely because that catalog is an unlucky draw. Our pre-registered design
internalizes the point: every test statistic is calibrated against null distributions
built from mock catalogs drawn from the frozen model itself, pushed through the same
injection-based selection, with finite-catalog fluctuations included following the
bootstrap logic of this paper, and an apparent discrepancy that sits inside the
catalog-variance null is explicitly not counted as a falsification.

## Inference with growing catalogs

The closest work in spirit is Wolfe, Mould, Veitch & Vitale (*Neural Bayesian updates
to populations with growing gravitational-wave catalogs*, arXiv:2602.20277,
2026-02-23), who use variational neural posterior estimation to update BBH population
posteriors sequentially as new data arrive, testing update cadences from
catalog-by-catalog to event-by-event across O4a. Their problem is computational
assimilation: how to incorporate new events efficiently and faithfully. Ours is
epistemic assessment: whether the pre-update model, held fixed, would have predicted
the new data at all. The two are complementary: sequential updating quantifies how
much each catalog *changes* the posterior, whereas a frozen-fit predictive score
quantifies whether the change was a *surprise* under the model's own predictive
distribution, but Wolfe et al. perform no held-out predictive scoring, and their
updates are still likelihood refits rather than locked forecasts.

## In-sample reanalyses of the O4-era catalogs

The months following each new catalog release have produced a substantial wave of
population reanalyses, all of which refit (or nonparametrically reconstruct) the
cumulative data. Following GWTC-4.0 in 2025: claims of a bimodal mass distribution
(arXiv:2508.20787), data-driven trend analyses (arXiv:2509.09876), evidence for three
mass-segregated subpopulations (arXiv:2509.15646), trimodality in chirp mass
(arXiv:2510.07573), flexible-mixture inference over 150 events (arXiv:2510.25579), and
binned-Gaussian-process characterizations of the 35 solar-mass feature
(arXiv:2511.22093). In the 2026 window most relevant to this project: Farah, Vijaykumar
& Fishbach (arXiv:2601.03456, 2026-01-06) attribute the observed redshift–effective-spin
correlation to a steeply evolving hierarchical-merger rate; data-driven subpopulation
searches with dimensionality reduction (arXiv:2603.06566, 2026-03); Gennari, Bertheas &
Tamanini (*Emergent structure in the binary black hole mass distribution and
implications for population-based cosmology*, arXiv:2604.14290, 2026-04-15) reconstruct
the primary-mass spectrum with B-splines, find a hierarchy of mass features, and trace
their impact on spectral-siren cosmology; Chatterjee (arXiv:2604.20941, 2026-04-22)
compresses GWTC-4 population fits into closed-form analytic laws via symbolic
regression; Agapito, De Renzis & Mancarella (arXiv:2605.20112, 2026-05-19) test
(putative) redshift evolution of the mass spectrum on GWTC-4.0 and find current
spectral-siren H0 constraints robust to it; Hussain, Isi & Zimmerman (*Evidence for
mass-dependent spin subpopulations in GWTC-4*, arXiv:2605.24281, 2026-05-22) find two
spin subpopulations whose mixing fraction varies with mass; formation-channel inference
on GWTC-4 (arXiv:2606.00234, 2026-05-29) and delay-time-distribution subpopulation
claims (arXiv:2606.02318, 2026-06-01) continue the same pattern. The first community
reanalysis of GWTC-5.0 has now appeared as well: Alvarez-Lopez, Heinzel & Vitale
(*Evidence for additional structure in the effective spin distribution hints at
multiple formation pathways in GWTC-5.0*, arXiv:2606.12205, 2026-06-10) jointly model
effective spin and primary mass and find structure beyond a non-skewed Gaussian bulk.

Five further GWTC-4-era in-sample studies (astro-ph.HE, abstract-mentions-GWTC, sorted
by submission date) round out this strand, all of which refit or reinterpret the
cumulative catalog rather than scoring a frozen forecast: Galaudage (*Compactness Peaks
and Subpopulations: Probing Stellar Physics and
Formation Channels of Merging Binary Black Holes*, arXiv:2605.25994, 2026-05-25); Godfrey
et al. (*A Strongly Parametrized Mass Ratio Model for the Stable Mass Transfer Channel: a
Case Study of the 10 solar-mass Peak*, arXiv:2605.23083, 2026-05-21);
Schiebelbein-Zwack et al. (*Forbidden Formation Histories: The Binary Black Hole Merger
Rate Disfavors Long Delay Times*, arXiv:2605.12858, 2026-05-13); Li et al.
(*Secondary-Mass Features improve Spectral-Siren H0 Constraints*, arXiv:2605.11474,
2026-05-12); and Islam (*Inference of recoil kicks from binary black hole mergers up to
GWTC-4 and their astrophysical implications*, arXiv:2604.04546, 2026-04-06). None performs
a held-out predictive test; they extend the same fit-describe-refit pattern.

Every entry in this strand answers the question "what does the cumulative catalog say
now?", and several of the headline features they report (the 35 solar-mass
overdensity, spin subpopulations, redshift-dependent structure) are precisely the kinds
of departures from the O3-era Power-Law+Peak description that an out-of-sample test is
sensitive to. Whether those features constitute a *predictive failure* of the O3-era
model, or sit comfortably within its forecast once selection effects and catalog
variance are accounted for, is not answered by refitting; it requires scoring the
frozen model. That is the question this project isolates.

## Positioning relative to prior work

To the best of our knowledge, no published or preprinted work freezes a pre-O4
population fit and scores the GWTC-4.0 or GWTC-5.0 catalogs against it prospectively,
with test statistics and falsification thresholds specified before any score was
computed. A search of arXiv (astro-ph.HE and gr-qc; arXiv API title/abstract queries
combining GWTC, out-of-sample, posterior predictive, forecast, held-out,
cross-validation, and frozen, covering all GWTC-4/GWTC-5-mentioning submissions in 2026,
last run 2026-06-14) returns in-sample predictive checks
(Miller et al.), sequential refits (Wolfe et al.), catalog-variance analyses (Corelli
et al.), and a large family of cumulative reanalyses, but no locked-forecast scoring
of the new catalogs. The contribution of this project is therefore twofold: the
pre-registered two-epoch verdict itself (validation or falsification, either being
informative), and a reusable, selection-matched, catalog-variance-calibrated
cross-catalog scoring framework that can be re-run on every future catalog release.

## Reference list

- LVK Collaboration, *The population of merging compact binaries inferred using
  gravitational waves through GWTC-3*, arXiv:2111.03634; Phys. Rev. X 13, 011048 (2023).
- I. M. Romero-Shaw, E. Thrane, P. D. Lasky, *When models fail: an introduction to
  posterior predictive checks and model misspecification in gravitational-wave
  astronomy*, arXiv:2202.05479; PASA 39, e025 (2022).
- LVK Collaboration, *GWTC-4.0: Population Properties of Merging Compact Binaries*,
  arXiv:2508.18083 (2025-08-25).
- A. M. Farah, A. Vijaykumar, M. Fishbach, *The steep redshift evolution of the
  hierarchical binary black hole merger rate may cause the z–χeff correlation*,
  arXiv:2601.03456 (2026-01-06).
- N. E. Wolfe, M. Mould, J. Veitch, S. Vitale, *Neural Bayesian updates to populations
  with growing gravitational-wave catalogs*, arXiv:2602.20277 (2026-02-23).
- A. Corelli, D. Gerosa, M. Mould, C. M. Fabbri, *Variance of gravitational-wave
  populations*, arXiv:2603.00239 (2026-02-27).
- *Data-Driven Trends and Subpopulations in the Gravitational Wave Binary Black Hole
  Merger Population with UMAP*, arXiv:2603.06566 (2026-03).
- S. J. Miller, S. Winney, K. Chatziioannou, P. M. Meyers, *Posterior Predictive Checks
  for Gravitational-wave Populations: Limitations and Improvements*, arXiv:2604.06090
  (2026-04-07).
- V. Gennari, T. Bertheas, N. Tamanini, *Emergent structure in the binary black hole
  mass distribution and implications for population-based cosmology*, arXiv:2604.14290
  (2026-04-15).
- C. Chatterjee, *Interpretable Analytic Formulae for GWTC-4 Binary Black Hole
  Population Properties via Symbolic Regression*, arXiv:2604.20941 (2026-04-22).
- A. Agapito, V. De Renzis, M. Mancarella, *Gravitational-wave constraints on H0 are
  robust to (putative) redshift evolution in the binary black hole mass spectrum at
  current sensitivity*, arXiv:2605.20112 (2026-05-19).
- A. Hussain, M. Isi, A. Zimmerman, *Evidence for mass-dependent spin subpopulations in
  GWTC-4*, arXiv:2605.24281 (2026-05-22).
- LVK Collaboration, *GWTC-5.0* catalog and population papers, arXiv:2605.27225,
  arXiv:2605.27226 (2026-05-26).
- *BBH-Genesis: Disentangling Binary Black Hole Formation Channels with GWTC-4*,
  arXiv:2606.00234 (2026-05-29).
- *The First Detection of Sub-Populations in the Delay-Time Distribution of Binary
  Black Holes*, arXiv:2606.02318 (2026-06-01).
- S. Alvarez-Lopez, J. Heinzel, S. Vitale, *Evidence for additional structure in the
  effective spin distribution hints at multiple formation pathways in GWTC-5.0*,
  arXiv:2606.12205 (2026-06-10).
- GWTC-4.0-era in-sample reanalyses (2025): arXiv:2508.20787, arXiv:2509.09876,
  arXiv:2509.15646, arXiv:2510.07573, arXiv:2510.25579, arXiv:2511.22093.
- Further GWTC-4.0-era in-sample studies: S. Galaudage,
  *Compactness Peaks and Subpopulations: Probing Stellar Physics and Formation Channels of
  Merging Binary Black Holes*, arXiv:2605.25994 (2026-05-25); A. Godfrey et al., *A Strongly
  Parametrized Mass Ratio Model for the Stable Mass Transfer Channel: a Case Study of the
  10 solar-mass Peak*, arXiv:2605.23083 (2026-05-21); Schiebelbein-Zwack et al., *Forbidden
  Formation Histories: The Binary Black Hole Merger Rate Disfavors Long Delay Times*,
  arXiv:2605.12858 (2026-05-13); Li et al., *Secondary-Mass Features improve Spectral-Siren
  H0 Constraints*, arXiv:2605.11474 (2026-05-12); T. Islam, *Inference of recoil kicks from
  binary black hole mergers up to GWTC-4 and their astrophysical implications*,
  arXiv:2604.04546 (2026-04-06).
