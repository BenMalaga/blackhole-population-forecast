"""Pre-registered scoring rules (PRE_REGISTRATION.md §1, §7) + the integrity guard.

H1 (distribution shape): posterior-averaged PITs per event/observable against the
frozen detected-population predictive CDF; KS + AD statistics of the PIT set against
uniformity; decision threshold = 97.5th percentile of the mock-null distribution
(mocks drawn from the frozen model through the locked selection, with measurement-
noise emulation and finite-catalog variance, pre-reg §1.3, §7).

H2 (detected count): observed N versus the central 95% of
``P(N) = (1/S) sum_s Poisson(N; R_s * VT_s)``.

Per-event surprise index (descriptive, §1.6): highest-density tail probability of
each event's (m1, z) under the frozen detected-population density, posterior-averaged.

INTEGRITY GUARD (pre-reg header + §10): any *real-data* score requires
``PRE_REGISTRATION.md`` to exist, carry no draft markers, be committed to git, and
be clean in the worktree, checked at runtime by :func:`require_locked_prereg`.
In the current stage (Stage 3, synthetic-only), the real entrypoint additionally
refuses to run without an explicit main-run go-signal: no real predictive score is
permitted yet.
"""

from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .predict import NoiseEmulator, PosteriorPredictor, PredictionResult

ROOT = Path(__file__).resolve().parent.parent
PREREG = ROOT / "PRE_REGISTRATION.md"

#: Locked decision thresholds (PRE_REGISTRATION.md §1).
H1_NULL_PERCENTILE = 97.5
H2_CENTRAL_COVERAGE = 0.95
LOCKED_SEED = 20260610


# ===================================================================== lock guard

class PreRegistrationError(RuntimeError):
    """Raised when scoring is attempted without a locked, committed pre-registration."""


def require_locked_prereg() -> str:
    """Verify PRE_REGISTRATION.md is locked + committed; return its commit hash.

    Mirrors src/fetch_data.py's firewall and the pre-reg header requirement:
    'Scoring code must re-verify at runtime that this file is committed before
    producing any score.'
    """
    if not PREREG.exists():
        raise PreRegistrationError("PRE_REGISTRATION.md does not exist, no real "
                                   "score may be produced (pre-reg firewall).")
    head = PREREG.read_text(encoding="utf-8", errors="replace").splitlines()[0].upper()
    if "DRAFT" in head or "NOT LOCKED" in head:
        raise PreRegistrationError("PRE_REGISTRATION.md is marked draft/unlocked.")
    try:
        commit = subprocess.run(
            ["git", "log", "-1", "--format=%H", "--", PREREG.name],
            cwd=ROOT, capture_output=True, text=True, check=True,
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "status", "--porcelain", "--", PREREG.name],
            cwd=ROOT, capture_output=True, text=True, check=True,
        ).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        raise PreRegistrationError(f"cannot verify git state of PRE_REGISTRATION.md: {exc}")
    if not commit:
        raise PreRegistrationError("PRE_REGISTRATION.md is not committed to git.")
    if dirty:
        raise PreRegistrationError("PRE_REGISTRATION.md has uncommitted local "
                                   f"changes ({dirty!r}), refusing to score.")
    return commit


# ===================================================================== statistics

def ks_uniform(pits: np.ndarray) -> float:
    """One-sample Kolmogorov–Smirnov distance of ``pits`` against U(0,1)."""
    u = np.sort(np.asarray(pits, dtype=float))
    n = len(u)
    if n == 0:
        raise ValueError("empty PIT set")
    i = np.arange(1, n + 1)
    return float(np.max(np.maximum(i / n - u, u - (i - 1) / n)))


def ad_uniform(pits: np.ndarray) -> float:
    """One-sample Anderson–Darling statistic of ``pits`` against U(0,1)."""
    eps = 1e-12
    u = np.clip(np.sort(np.asarray(pits, dtype=float)), eps, 1 - eps)
    n = len(u)
    if n == 0:
        raise ValueError("empty PIT set")
    i = np.arange(1, n + 1)
    s = np.sum((2 * i - 1) * (np.log(u) + np.log1p(-u[::-1])))
    return float(-n - s / n)


def posterior_averaged_pit(event_samples: np.ndarray, pred_sorted: np.ndarray,
                           pred_cum: np.ndarray) -> float:
    """``PIT_i = (1/K) sum_k F_pred(u_ik)`` (pre-reg §1.1).

    ``pred_sorted``/``pred_cum`` come from ``PredictionResult.predictive_cdf``.
    """
    x = np.asarray(event_samples, dtype=float)
    idx = np.searchsorted(pred_sorted, x, side="right")
    cdf = np.where(idx > 0, pred_cum[np.minimum(idx, len(pred_cum)) - 1], 0.0)
    return float(np.mean(cdf))


def pit_set(event_sample_sets: list[np.ndarray], pred_sorted: np.ndarray,
            pred_cum: np.ndarray) -> np.ndarray:
    return np.array([posterior_averaged_pit(s, pred_sorted, pred_cum)
                     for s in event_sample_sets])


def surprise_index(event_samples_m1z: np.ndarray, pred_m1z: np.ndarray,
                   pred_weights: np.ndarray) -> float:
    """Descriptive per-event surprise (pre-reg §1.6).

    ``s_i = (1/K) sum_k Pr_{x~pred}[ f_pred(x) <= f_pred(x_ik) ]`` with f_pred a
    weighted Gaussian KDE of the detected-population predictive in
    (log m1, z), the density-estimation detail is internal/descriptive (no
    accept/reject threshold is attached to this number).
    """
    from scipy.stats import gaussian_kde

    pred = np.column_stack([np.log(pred_m1z[:, 0]), pred_m1z[:, 1]])
    w = np.asarray(pred_weights, dtype=float)
    w = w / w.sum()
    kde = gaussian_kde(pred.T, weights=w)
    f_pred_at_pred = kde(pred.T)
    ev = np.column_stack([np.log(event_samples_m1z[:, 0]), event_samples_m1z[:, 1]])
    f_at_event = kde(ev.T)
    # Pr[f(x) <= f(x_k)] under the predictive, weighted
    order = np.argsort(f_pred_at_pred)
    cum_w = np.cumsum(w[order])
    idx = np.searchsorted(f_pred_at_pred[order], f_at_event, side="right")
    tail = np.where(idx > 0, cum_w[np.minimum(idx, len(cum_w)) - 1], 0.0)
    return float(np.mean(tail))


# ===================================================================== H1 machinery

@dataclass
class H1Observable:
    observable: str
    statistic: str               # "ks" | "ad"
    observed: float
    null_threshold: float        # 97.5th percentile of the mock-null
    null_values: np.ndarray
    falsified: bool


def mock_null_distribution(predictor: PosteriorPredictor, result: PredictionResult,
                           emulator: NoiseEmulator, observables: tuple[str, ...],
                           n_mocks: int, rng: np.random.Generator,
                           n_events: int | None = None,
                           ) -> dict[str, dict[str, np.ndarray]]:
    """Null distributions of KS/AD per observable from frozen-model mock catalogs.

    Each mock draws a hyper-sample (hyper variance), a catalog size (Poisson if
    ``n_events`` is None, catalog variance per §1.3), event truths through the
    locked selection, and emulated PE posteriors; the identical PIT operator is
    applied (pre-reg §1.1: 'the null distribution of every statistic is exact by
    construction').
    """
    cdfs = {obs: result.predictive_cdf(obs) for obs in observables}
    nulls = {obs: {"ks": np.empty(n_mocks), "ad": np.empty(n_mocks)}
             for obs in observables}
    for m in range(n_mocks):
        catalog, _ = predictor.draw_mock_catalog(result, rng, n_events=n_events)
        emulated = [emulator.emulate(row, observables, rng)
                    for _, row in catalog.iterrows()]
        for obs in observables:
            pred_sorted, pred_cum = cdfs[obs]
            pits = pit_set([e[obs] for e in emulated], pred_sorted, pred_cum)
            nulls[obs]["ks"][m] = ks_uniform(pits)
            nulls[obs]["ad"][m] = ad_uniform(pits)
    return nulls


def score_h1(observed_event_samples: dict[str, list[np.ndarray]],
             result: PredictionResult,
             nulls: dict[str, dict[str, np.ndarray]]) -> list[H1Observable]:
    """Apply the locked H1 decision rule to a set of events.

    ``observed_event_samples``: observable -> list of per-event sample arrays
    (real PE samples post-unblinding; emulated samples in synthetic tests).
    """
    out: list[H1Observable] = []
    for obs, sample_sets in observed_event_samples.items():
        pred_sorted, pred_cum = result.predictive_cdf(obs)
        pits = pit_set(sample_sets, pred_sorted, pred_cum)
        for stat_name, fn in (("ks", ks_uniform), ("ad", ad_uniform)):
            value = fn(pits)
            null = nulls[obs][stat_name]
            thresh = float(np.percentile(null, H1_NULL_PERCENTILE))
            out.append(H1Observable(observable=obs, statistic=stat_name,
                                    observed=value, null_threshold=thresh,
                                    null_values=null,
                                    falsified=bool(value > thresh)))
    return out


# ===================================================================== H2 machinery

@dataclass
class H2Result:
    observed_n: int
    interval: tuple[int, int]    # central-95% predictive interval [lo, hi]
    pmf: np.ndarray
    falsified: bool


def score_h2(observed_n: int, mu_samples: np.ndarray) -> H2Result:
    """Locked H2 rule: falsified iff observed N outside the central 95% of P(N)."""
    pmf = PosteriorPredictor.count_pmf(np.asarray(mu_samples, dtype=float))
    cdf = np.cumsum(pmf)
    alpha = (1.0 - H2_CENTRAL_COVERAGE) / 2.0
    lo = int(np.searchsorted(cdf, alpha))
    hi = int(np.searchsorted(cdf, 1.0 - alpha))
    falsified = bool(observed_n < lo or observed_n > hi)
    return H2Result(observed_n=int(observed_n), interval=(lo, hi), pmf=pmf,
                    falsified=falsified)


# ===================================================================== real entrypoint

def score_epoch_real(*args, **kwargs):  # pragma: no cover - guarded, not used in Stage 3
    """Real-data scoring entrypoint: BLOCKED until the explicit main-run go.

    Stage-3 rule (and analysis-sequence §8): the mock-calibration gate (§7) must
    pass and the H2 predictive distributions must be committed before unblinding.
    The guard chain:

    1. ``require_locked_prereg()``: locked + committed + clean.
    2. ``BHPF_ORCHESTRATOR_GO=1`` in the environment, set only by the deliberate
       main-run directive, never by tests.
    3. The mock-calibration report must exist (``results/mock_calibration/``).
    """
    commit = require_locked_prereg()
    if os.environ.get("BHPF_ORCHESTRATOR_GO") != "1":
        raise PreRegistrationError(
            "real-data scoring is blocked in this stage: the mock-calibration "
            "gate (§7) has not been run. "
            "Set BHPF_ORCHESTRATOR_GO=1 only at the deliberate main-run go."
        )
    gate = ROOT / "results" / "mock_calibration"
    if not gate.exists():
        raise PreRegistrationError(
            "mock-calibration gate report (results/mock_calibration/) missing, "
            "pre-reg §7: if the gate has not passed, the main run does not start."
        )
    raise NotImplementedError(
        f"main-run wiring lands at analysis step §8.4 (lock commit {commit[:8]}); "
        "this stage ships the machinery + synthetic validation only."
    )
