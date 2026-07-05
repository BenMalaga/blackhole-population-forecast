"""Synthetic world for end-to-end pipeline tests (Stage 3, NO real data touched).

Everything here is fabricated:

* a known truth Λ* for the locked model family (Power-Law+Peak + iid spins +
  power-law redshift),
* a tiny hyperposterior = Λ* + small jitter (and a degenerate all-truth variant),
* a fake injection set in the *exact* parquet schema of the real sensitivity
  releases (mass1_source/mass2_source/redshift/cartesian spins, a summed
  ``lnpdraw_*`` field, mixture ``weights``, per-search ``*_far`` columns,
  ``time_geocenter``), with an analytically known draw density,
* fake "observed" events drawn from the truth population *independently* of the
  injection set (grid inverse-CDF sampling), pushed through the same deterministic
  detection rule, with emulated PE posteriors.

The known answers the pipeline must recover are asserted in
``test_pipeline_synthetic.py``. Locked seed 20260610 (pre-reg §9).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.population_model import FrozenPopulationModel
from src.predict import LNPDRAW_MAIN, NoiseEmulator, PosteriorPredictor, SelectionSet

SEED = 20260610

#: Known truth hyperparameters (canonical gwpopulation names).
TRUTH = {
    "alpha": 3.0, "beta": 1.2, "mmin": 5.0, "mmax": 85.0,
    "lam": 0.04, "mpp": 34.0, "sigpp": 4.0, "delta_m": 4.5,
    "amax": 1.0, "alpha_chi": 2.0, "beta_chi": 5.0,
    "xi_spin": 0.7, "sigma_spin": 0.6, "lamb": 2.5,
}
Z_MAX = 2.0
RATE_TRUE = 60.0           # arbitrary synthetic rate units (consistent both sides)
N_INJECTIONS = 24_000      # all "generated" -> total_generated = N_INJECTIONS
T_ANALYSIS = 1.0           # synthetic analysis time (units cancel by construction)

# Draw-density bounds (broad reference distribution, covers the truth support)
M1_LO, M1_HI = 3.0, 150.0
Z_LO = 1e-3


def _draw_density_and_samples(rng: np.random.Generator, n: int):
    """Reference draws + their exact ln density in the lnpdraw measure.

    p_draw(m1) log-uniform on [M1_LO, M1_HI]; m2 | m1 uniform on [M1_LO, m1];
    z with p(z) ∝ (1+z)^2 on [Z_LO, Z_MAX]; spins isotropic with uniform
    magnitude on [0, 1) -> cartesian density 1/(4 pi a^2) per spin.
    """
    m1 = np.exp(rng.uniform(np.log(M1_LO), np.log(M1_HI), n))
    m2 = rng.uniform(M1_LO, m1)
    # z via inverse CDF of (1+z)^2
    u = rng.uniform(size=n)
    a3, b3 = (1 + Z_LO) ** 3, (1 + Z_MAX) ** 3
    z = ((a3 + u * (b3 - a3)) ** (1.0 / 3.0)) - 1.0

    def iso_spin(n):
        a = rng.uniform(0.0, 1.0, n)
        cos_t = rng.uniform(-1.0, 1.0, n)
        phi = rng.uniform(0.0, 2 * np.pi, n)
        sin_t = np.sqrt(1 - cos_t ** 2)
        return (a * sin_t * np.cos(phi), a * sin_t * np.sin(phi), a * cos_t, a)

    s1x, s1y, s1z, a1 = iso_spin(n)
    s2x, s2y, s2z, a2 = iso_spin(n)

    ln_p_m1 = -np.log(m1) - np.log(np.log(M1_HI / M1_LO))
    ln_p_m2 = -np.log(m1 - M1_LO)
    ln_p_z = 2 * np.log1p(z) - np.log((b3 - a3) / 3.0)
    ln_p_s = -np.log(4 * np.pi * a1 ** 2) - np.log(4 * np.pi * a2 ** 2)
    lnpdraw = ln_p_m1 + ln_p_m2 + ln_p_z + ln_p_s

    df = pd.DataFrame({
        "mass1_source": m1, "mass2_source": m2, "redshift": z,
        "spin1x": s1x, "spin1y": s1y, "spin1z": s1z,
        "spin2x": s2x, "spin2y": s2y, "spin2z": s2z,
        LNPDRAW_MAIN: lnpdraw,
        "weights": np.ones(n),
        "time_geocenter": rng.uniform(0.0, 1.0e6, n),
    })
    return df


def detection_far(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Deterministic two-search detection rule (mimics min-over-*_far semantics).

    'SNR proxy' grows with chirp mass and falls with distance; search A detects
    above 1.0, search B (slightly different weighting) above 1.05. far = 0.1
    (detected) or 1e6 (not), only the min < 1/yr cut matters.
    """
    mc = (df["mass1_source"] * df["mass2_source"]) ** 0.6 / \
         (df["mass1_source"] + df["mass2_source"]) ** 0.2
    snr = (mc / 15.0) ** (5.0 / 6.0) / (df["redshift"] / 0.5)
    far_a = np.where(snr > 1.0, 0.1, 1.0e6)
    far_b = np.where(0.95 * snr > 1.0, 0.1, 1.0e6)
    return np.asarray(far_a, dtype=float), np.asarray(far_b, dtype=float)


def make_injection_set(rng: np.random.Generator) -> SelectionSet:
    df = _draw_density_and_samples(rng, N_INJECTIONS)
    far_a, far_b = detection_far(df)
    df["mock_a_far"], df["mock_b_far"] = far_a, far_b
    return SelectionSet(df=df, total_analysis_time=T_ANALYSIS,
                        total_generated=float(N_INJECTIONS),
                        lnpdraw_columns=(LNPDRAW_MAIN,), name="synthetic")


def make_hyperposterior(rng: np.random.Generator, n_samples: int = 48,
                        jitter: float = 0.03, rate_scatter: float = 0.10
                        ) -> pd.DataFrame:
    """Λ* + small relative jitter, clipped to safe ranges; lognormal rate."""
    rows = {}
    for k, v in TRUTH.items():
        if k == "amax":
            rows[k] = np.full(n_samples, v)
        else:
            rows[k] = v * (1.0 + jitter * rng.standard_normal(n_samples))
    rows["lam"] = np.clip(rows["lam"], 0.005, 0.3)
    rows["xi_spin"] = np.clip(rows["xi_spin"], 0.05, 0.99)
    rows["mmin"] = np.clip(rows["mmin"], 3.5, 8.0)
    rows["rate"] = RATE_TRUE * np.exp(rate_scatter * rng.standard_normal(n_samples)
                                      - 0.5 * rate_scatter ** 2)
    return pd.DataFrame(rows)


def sample_truth_events_direct(model: FrozenPopulationModel, rng: np.random.Generator,
                               n_target: int) -> pd.DataFrame:
    """Sample events from the TRUTH population independently of the injection set.

    (m1, q) by 2D-grid inverse-CDF on the gwpopulation density at Λ* (sample 0 of
    a degenerate hyperposterior), z by 1D-grid inverse-CDF, spins analytically
    (Beta magnitudes; isotropic+truncated-Gaussian tilt mixture). Keeps only
    events passing the deterministic detection rule. Returns >= n_target rows
    (drawn in batches), truncated to n_target.
    """
    kept = []
    n_kept = 0
    while n_kept < n_target:
        batch = _sample_truth_batch(model, rng, 4000)
        far_a, far_b = detection_far(batch)
        det = batch.loc[np.minimum(far_a, far_b) < 1.0]
        kept.append(det)
        n_kept += len(det)
    return pd.concat(kept, ignore_index=True).iloc[:n_target].copy()


def truth_detection_probability(model: FrozenPopulationModel,
                                rng: np.random.Generator, n: int = 60_000) -> float:
    """Direct Monte-Carlo P_det under Λ*, the independent 'known answer' for VT."""
    batch = _sample_truth_batch(model, rng, n)
    far_a, far_b = detection_far(batch)
    return float(np.mean(np.minimum(far_a, far_b) < 1.0))


def _sample_truth_batch(model: FrozenPopulationModel, rng: np.random.Generator,
                        n: int) -> pd.DataFrame:
    from scipy.stats import beta as beta_dist

    # --- (m1, q) on a grid
    m1_grid = np.linspace(TRUTH["mmin"] * 0.8, TRUTH["mmax"] * 1.05, 480)
    q_grid = np.linspace(0.05, 1.0, 240)
    m1g, qg = np.meshgrid(m1_grid, q_grid, indexing="ij")
    mass_model = model._mass_model()
    p = np.asarray(mass_model(
        {"mass_1": m1g.ravel(), "mass_ratio": qg.ravel()},
        alpha=TRUTH["alpha"], beta=TRUTH["beta"], mmin=TRUTH["mmin"],
        mmax=TRUTH["mmax"], lam=TRUTH["lam"], mpp=TRUTH["mpp"],
        sigpp=TRUTH["sigpp"], delta_m=TRUTH["delta_m"]), dtype=float
    ).reshape(m1g.shape)
    p = np.where(np.isfinite(p) & (p > 0), p, 0.0)
    flat = p.ravel() / p.sum()
    idx = rng.choice(flat.size, size=n, p=flat)
    i1, i2 = np.unravel_index(idx, p.shape)
    dm1 = m1_grid[1] - m1_grid[0]
    dq = q_grid[1] - q_grid[0]
    m1 = m1_grid[i1] + rng.uniform(-0.5, 0.5, n) * dm1
    q = np.clip(q_grid[i2] + rng.uniform(-0.5, 0.5, n) * dq, 0.02, 1.0)
    m2 = q * m1

    # --- z on a grid
    z_grid = np.linspace(1e-3, Z_MAX, 600)
    pz = np.asarray(model._redshift_model()({"redshift": z_grid},
                                            lamb=TRUTH["lamb"]), dtype=float)
    pz = np.where(np.isfinite(pz) & (pz > 0), pz, 0.0)
    zi = rng.choice(z_grid.size, size=n, p=pz / pz.sum())
    z = np.clip(z_grid[zi] + rng.uniform(-0.5, 0.5, n) * (z_grid[1] - z_grid[0]),
                1e-4, Z_MAX)

    # --- spins: iid Beta magnitudes; EVENT-level tilt mixture (matches the LVK
    # Default / gwpopulation model: with prob xi BOTH tilts ~ TruncNorm(1, sigma),
    # else BOTH isotropic, not an independent per-spin mixture)
    a1 = beta_dist.rvs(TRUTH["alpha_chi"], TRUTH["beta_chi"], size=n,
                       random_state=rng)
    a2 = beta_dist.rvs(TRUTH["alpha_chi"], TRUTH["beta_chi"], size=n,
                       random_state=rng)
    cos_t1, cos_t2 = _tilt_mixture_pair(rng, n)
    phi1 = rng.uniform(0, 2 * np.pi, n)
    phi2 = rng.uniform(0, 2 * np.pi, n)
    s1 = np.sqrt(1 - cos_t1 ** 2)
    s2 = np.sqrt(1 - cos_t2 ** 2)
    return pd.DataFrame({
        "mass1_source": m1, "mass2_source": m2, "redshift": z,
        "spin1x": a1 * s1 * np.cos(phi1), "spin1y": a1 * s1 * np.sin(phi1),
        "spin1z": a1 * cos_t1,
        "spin2x": a2 * s2 * np.cos(phi2), "spin2y": a2 * s2 * np.sin(phi2),
        "spin2z": a2 * cos_t2,
    })


def _tilt_mixture_pair(rng: np.random.Generator, n: int
                       ) -> tuple[np.ndarray, np.ndarray]:
    from scipy.stats import truncnorm

    sig = TRUTH["sigma_spin"]

    def aligned(k):
        return truncnorm.rvs((-1 - 1) / sig, (1 - 1) / sig, loc=1.0, scale=sig,
                             size=k, random_state=rng)

    pick = rng.uniform(size=n) < TRUTH["xi_spin"]   # one pick per EVENT
    c1 = np.where(pick, aligned(n), rng.uniform(-1, 1, n))
    c2 = np.where(pick, aligned(n), rng.uniform(-1, 1, n))
    return c1, c2


# ===================================================================== fixtures

@pytest.fixture(scope="session")
def rng():
    return np.random.default_rng(SEED)


@pytest.fixture(scope="session")
def world(rng):
    """The full synthetic world, built once per test session."""
    injections = make_injection_set(rng)
    hyper = make_hyperposterior(rng)
    model = FrozenPopulationModel(hyperposterior=hyper, z_max=Z_MAX, mmax_grid=160.0)
    truth_model = FrozenPopulationModel(
        hyperposterior=pd.DataFrame({k: [v] for k, v in {**TRUTH, "rate": RATE_TRUE}.items()}),
        z_max=Z_MAX, mmax_grid=160.0)
    predictor = PosteriorPredictor(model=model, selection=injections)
    emulator = NoiseEmulator(
        sigmas={"mass1_source": 0.08, "redshift": 0.12},
        transforms={"mass1_source": "log", "redshift": "log"},
        n_samples=96,
    )
    return {
        "injections": injections,
        "hyper": hyper,
        "model": model,
        "truth_model": truth_model,
        "predictor": predictor,
        "emulator": emulator,
    }
