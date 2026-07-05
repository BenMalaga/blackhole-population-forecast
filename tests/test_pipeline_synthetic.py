"""End-to-end synthetic validation of the scoring pipeline (Stage 3).

NO real data is read anywhere in this suite, the world is fabricated in
``conftest.py``. The 'known answers' the pipeline must recover:

1. **VT / count machinery:** the locked importance-sampling VT estimator,
   evaluated at the truth Λ*, must agree with an independent direct Monte-Carlo
   estimate of the detected fraction (grid inverse-CDF sampling of the truth
   population + the same deterministic detection rule).
2. **H2 coverage:** a detected count realized from the truth must fall inside the
   central-95% predictive interval built from the jittered hyperposterior.
3. **H1 calibration:** observed events drawn FROM the truth (with emulated PE
   noise) must NOT falsify H1, KS/AD of the PITs stay below the 97.5th
   percentile of the mock-null.
4. **H1 power:** grossly shifted observed events (m1 x 1.6) MUST falsify H1 on m1.
5. Unit-level checks of the PIT operator, KS/AD statistics, the count pmf, and
   the pre-registration lock guard.

Locked seed 20260610; the whole suite is CPU-only and runs in minutes.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src import score as score_mod
from src.population_model import FrozenPopulationModel
from src.predict import PosteriorPredictor
from src.score import (H1_NULL_PERCENTILE, PreRegistrationError, ad_uniform,
                       ks_uniform, mock_null_distribution, pit_set,
                       posterior_averaged_pit, require_locked_prereg, score_h1,
                       score_h2, surprise_index)
from tests.conftest import (RATE_TRUE, SEED, T_ANALYSIS, TRUTH,
                            sample_truth_events_direct,
                            truth_detection_probability)

OBSERVABLES = ("mass1_source", "redshift")


# ===================================================================== unit level

def test_ks_ad_uniform_sanity():
    rng = np.random.default_rng(SEED)
    u = rng.uniform(size=5000)
    assert ks_uniform(u) < 0.03            # uniform sample -> tiny KS distance
    assert ad_uniform(u) < 4.0             # and small AD
    shifted = np.clip(u * 0.5, 0, 1)       # grossly non-uniform
    assert ks_uniform(shifted) > 0.3
    assert ad_uniform(shifted) > 50.0


def test_posterior_averaged_pit_matches_analytic_cdf():
    rng = np.random.default_rng(SEED)
    pred = np.sort(rng.normal(size=200_000))
    cum = np.arange(1, len(pred) + 1) / len(pred)
    from scipy.stats import norm
    for x in (-1.5, 0.0, 0.7):
        pit = posterior_averaged_pit(np.array([x]), pred, cum)
        assert abs(pit - norm.cdf(x)) < 5e-3


def test_count_pmf_is_poisson_mixture():
    mu = np.array([10.0, 20.0])
    pmf = PosteriorPredictor.count_pmf(mu, n_max=80)
    from scipy.stats import poisson
    expected = 0.5 * (poisson.pmf(15, 10.0) + poisson.pmf(15, 20.0))
    assert abs(pmf[15] - expected) < 1e-12
    assert abs(pmf.sum() - 1.0) < 1e-6


def test_h2_decision_rule():
    mu = np.full(64, 50.0)
    assert not score_h2(50, mu).falsified          # dead center
    assert score_h2(5, mu).falsified               # absurdly low
    assert score_h2(120, mu).falsified             # absurdly high
    res = score_h2(50, mu)
    lo, hi = res.interval
    assert lo < 50 < hi and abs(res.pmf.sum() - 1.0) < 1e-6


def test_spin_magnitude_convention_matches_scipy_beta(world):
    """gwpopulation's iid Beta magnitude must match scipy's Beta pdf (amax=1)."""
    from scipy.stats import beta as beta_dist
    from src.population_model import _import_gwpop

    _, _, spin_models = _import_gwpop()
    a = np.linspace(0.05, 0.95, 7)
    dataset = {"a_1": a, "a_2": a, "cos_tilt_1": np.zeros(7),
               "cos_tilt_2": np.zeros(7)}
    p = np.asarray(spin_models.iid_spin(
        dataset, xi_spin=0.0, sigma_spin=0.5, amax=1.0,
        alpha_chi=TRUTH["alpha_chi"], beta_chi=TRUTH["beta_chi"]), dtype=float)
    # xi=0 -> tilts isotropic (1/2 each); magnitudes iid Beta
    expected = (beta_dist.pdf(a, TRUTH["alpha_chi"], TRUTH["beta_chi"]) ** 2) / 4.0
    assert np.allclose(p, expected, rtol=5e-2)


def test_lock_guard_accepts_locked_prereg_and_blocks_real_scoring(monkeypatch):
    commit = require_locked_prereg()               # repo state: locked + committed
    assert len(commit) == 40
    # the real entrypoint must refuse without the explicit main-run go-signal
    monkeypatch.delenv("BHPF_ORCHESTRATOR_GO", raising=False)
    with pytest.raises(PreRegistrationError, match="blocked"):
        score_mod.score_epoch_real()


def test_lock_guard_rejects_missing_prereg(monkeypatch, tmp_path):
    monkeypatch.setattr(score_mod, "PREREG", tmp_path / "PRE_REGISTRATION.md")
    with pytest.raises(PreRegistrationError, match="does not exist"):
        require_locked_prereg()


# ===================================================================== end to end

@pytest.fixture(scope="module")
def prediction(world):
    """The forward model run once: per-sample VT/ESS + averaged predictive weights."""
    result = world["predictor"].run(n_events_for_ess=80)
    world["predictor"].expected_counts(result, world["model"].rate_samples)
    return result


@pytest.fixture(scope="module")
def truth_world(world):
    """Independent truth-side quantities (the 'known answer')."""
    rng = np.random.default_rng(SEED + 1)          # independent stream from fixtures
    p_det = truth_detection_probability(world["truth_model"], rng, n=60_000)
    mu_direct = RATE_TRUE * T_ANALYSIS * p_det
    return {"p_det": p_det, "mu_direct": mu_direct, "rng": rng}


def test_vt_estimator_recovers_direct_monte_carlo(world, prediction, truth_world):
    """Known answer #1: locked VT estimator == direct MC detected fraction at Λ*."""
    truth_predictor = PosteriorPredictor(model=world["truth_model"],
                                         selection=world["injections"])
    res = truth_predictor.run()
    vt_at_truth = res.vt[0]                        # T * P_det estimate
    p_det_injection = vt_at_truth / T_ANALYSIS
    p_det_direct = truth_world["p_det"]
    assert p_det_direct > 0.005, "synthetic detection rule selects almost nothing"
    assert abs(p_det_injection - p_det_direct) / p_det_direct < 0.15, (
        f"VT estimator P_det={p_det_injection:.4f} vs direct MC "
        f"P_det={p_det_direct:.4f}"
    )
    # diagnostics must be populated and sane
    assert res.neff[0] > 50, "importance ESS too small for a meaningful test"
    assert len(res.detected.df) > 200


def test_h2_count_coverage_under_truth(prediction, truth_world):
    """Known answer #2: a truth-realized count sits inside the central 95%."""
    n_obs = int(truth_world["rng"].poisson(truth_world["mu_direct"]))
    res = score_h2(n_obs, prediction.mu)
    assert not res.falsified, (
        f"H2 wrongly falsified: N={n_obs}, interval={res.interval}, "
        f"mu_direct={truth_world['mu_direct']:.1f}, "
        f"mu_pipeline median={np.median(prediction.mu):.1f}"
    )
    # and the pipeline's own mu must agree with the independent direct estimate
    ratio = np.median(prediction.mu) / truth_world["mu_direct"]
    assert 0.7 < ratio < 1.4, f"pipeline/direct count ratio {ratio:.2f}"


@pytest.fixture(scope="module")
def observed_catalogs(world, truth_world):
    """Fake observed events: truth draws + emulated PE posteriors (both epochs of
    the same world: one faithful, one with a gross m1 shift)."""
    rng = np.random.default_rng(SEED + 2)
    n_obs = max(int(truth_world["rng"].poisson(truth_world["mu_direct"])), 40)
    events = sample_truth_events_direct(world["truth_model"], rng, n_obs)
    emulator = world["emulator"]
    faithful = {obs: [] for obs in OBSERVABLES}
    shifted = {obs: [] for obs in OBSERVABLES}
    for _, row in events.iterrows():
        em = emulator.emulate(row, OBSERVABLES, rng)
        for obs in OBSERVABLES:
            faithful[obs].append(em[obs])
        row_shift = row.copy()
        row_shift["mass1_source"] = row["mass1_source"] * 1.6   # gross shift
        em_s = emulator.emulate(row_shift, OBSERVABLES, rng)
        for obs in OBSERVABLES:
            shifted[obs].append(em_s[obs])
    return {"faithful": faithful, "shifted": shifted, "n_obs": n_obs}


@pytest.fixture(scope="module")
def nulls(world, prediction, observed_catalogs):
    rng = np.random.default_rng(SEED + 3)
    return mock_null_distribution(
        world["predictor"], prediction, world["emulator"], OBSERVABLES,
        n_mocks=150, rng=rng, n_events=observed_catalogs["n_obs"],
    )


def test_h1_calibration_truth_events_pass(prediction, observed_catalogs, nulls):
    """Known answer #3: events from the truth do NOT falsify H1 (m1 and z)."""
    results = score_h1(observed_catalogs["faithful"], prediction, nulls)
    for r in results:
        assert not r.falsified, (
            f"H1 wrongly falsified for {r.observable}/{r.statistic}: "
            f"observed {r.observed:.3f} > {H1_NULL_PERCENTILE}th null pct "
            f"{r.null_threshold:.3f}"
        )


def test_h1_power_gross_m1_shift_detected(prediction, observed_catalogs, nulls):
    """Known answer #4: a 1.6x m1 shift MUST trip the locked m1 KS rule."""
    results = score_h1(observed_catalogs["shifted"], prediction, nulls)
    by_key = {(r.observable, r.statistic): r for r in results}
    r = by_key[("mass1_source", "ks")]
    assert r.falsified, (
        f"pipeline failed to detect a gross m1 shift: KS {r.observed:.3f} "
        f"<= threshold {r.null_threshold:.3f}"
    )
    assert by_key[("mass1_source", "ad")].falsified


def test_surprise_index_orders_outliers(world, prediction):
    """Descriptive surprise index: a wild outlier scores far more 'surprising'
    (smaller tail probability) than a bulk event."""
    det = prediction.detected.df
    pred_m1z = det[["mass1_source", "redshift"]].to_numpy()
    w = prediction.w_bar
    bulk = np.array([[np.median(det["mass1_source"]), np.median(det["redshift"])]])
    outlier = np.array([[140.0, 1.9]])
    s_bulk = surprise_index(np.repeat(bulk, 8, axis=0), pred_m1z, w)
    s_out = surprise_index(np.repeat(outlier, 8, axis=0), pred_m1z, w)
    assert s_out < s_bulk
    assert s_out < 0.05


def test_ess_diagnostic_reported(prediction):
    assert prediction.neff.shape == prediction.vt.shape
    assert np.all(prediction.neff >= 0)
    assert 0.0 <= prediction.flagged_fraction <= 1.0


def test_no_real_data_touched():
    """Stage-3 invariant: this suite must not read any real derived parquet."""
    import tests.conftest as c
    assert c.N_INJECTIONS <= 50_000
    # the synthetic world never points at data/derived
    import inspect
    src = inspect.getsource(c)
    assert "data/derived" not in src
