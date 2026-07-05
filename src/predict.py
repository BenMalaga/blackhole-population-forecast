"""Posterior-predictive machinery (PRE_REGISTRATION.md §5–§6, §1).

Frozen hyperposterior draws -> (a) the predicted detected-event-count distribution
per epoch (H2) and (b) detected-population predictive distributions in observable
space (H1), both through the public sensitivity-injection selection function.

Locked recipes implemented here, verbatim from the pre-registration:

* **Detection definition (§5):** an injection is detected iff
  ``min_over_searches(<search>_far) < 1/yr`` over all per-search ``*_far`` fields
  the release provides (absent detections -> ``inf``).
* **VT estimator (§5):**
  ``VT(Λ) = total_analysis_time * sum_sel[ weights * exp(lnp(Λ) - lnpdraw) ] / total_generated``
  with ``lnpdraw`` the release's summed draw-density fields and ``weights`` the
  release's mixture weights; root attrs supply ``total_analysis_time`` /
  ``total_generated``.
* **Count predictive (§1 H2):** ``P(N) = (1/S) * sum_s Poisson(N; R_s * VT_s)``.
* **Detected-population predictive density (§5):** the importance-reweighted
  injection distribution (normalization cancels); the H1 predictive CDF is the
  hyperposterior-averaged weighted ECDF.
* **ESS diagnostic (§5):** flag hyperposterior samples with
  ``Neff <= 4 * N_events``; report if > 1% fail.

Real-data conventions still to be pinned at forward-model build (correctness
questions, not tunable choices, pre-reg §2/§3/§5): the rate-normalization unit
chain (release rate posterior x VT units) and whether the inclination draw density
(`lnpdraw_inclination`) pairs with an isotropic population factor. Both are isolated
in single, documented spots below.

NO real predictive score is produced by importing or unit-testing this module; the
real-data entrypoints live in ``src/score.py`` behind the lock guard.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .population_model import FrozenPopulationModel

#: Seconds in a Julian year, FAR thresholds are quoted per year (LVK convention).
SECONDS_PER_YEAR = 365.25 * 24 * 3600.0

#: The locked combined draw-density field of the O3-era mixture / O4a / O4ab releases.
LNPDRAW_MAIN = ("lnpdraw_mass1_source_mass2_source_redshift_"
                "spin1x_spin1y_spin1z_spin2x_spin2y_spin2z")


# ===================================================================== selection set

@dataclass
class SelectionSet:
    """An injection (sensitivity) release: rows + the root attrs the VT formula needs.

    ``lnpdraw_columns`` are summed (locked: "the release's summed draw-density
    fields"). Which fields participate for the real epochs is pinned at
    forward-model build against the release notes; synthetic tests pass their own.
    """

    df: pd.DataFrame
    total_analysis_time: float          # seconds (release root attr, as released)
    total_generated: float
    lnpdraw_columns: tuple[str, ...] = (LNPDRAW_MAIN,)
    weights_column: str | None = "weights"
    name: str = "selection"

    @classmethod
    def from_parquet(cls, parquet_path: Path | str, manifest_attrs: dict,
                     **kwargs) -> "SelectionSet":
        df = pd.read_parquet(parquet_path)
        return cls(df=df,
                   total_analysis_time=float(manifest_attrs["total_analysis_time"]),
                   total_generated=float(manifest_attrs["total_generated"]),
                   name=str(parquet_path), **kwargs)

    @property
    def far_columns(self) -> list[str]:
        return [c for c in self.df.columns if c.endswith("_far")]

    def far_min(self) -> np.ndarray:
        """Min over per-search FAR fields; NaN/absent -> inf (locked §5)."""
        cols = self.far_columns
        if not cols:
            raise ValueError(f"{self.name}: no *_far columns found")
        vals = self.df[cols].to_numpy(dtype=float)
        vals = np.where(np.isfinite(vals) & (vals >= 0), vals, np.inf)
        return vals.min(axis=1)

    def detected(self, far_threshold_per_yr: float = 1.0) -> "SelectionSet":
        """Subset to detected injections (FAR < threshold). total_* attrs unchanged."""
        mask = self.far_min() < far_threshold_per_yr
        return SelectionSet(df=self.df.loc[mask].reset_index(drop=True),
                            total_analysis_time=self.total_analysis_time,
                            total_generated=self.total_generated,
                            lnpdraw_columns=self.lnpdraw_columns,
                            weights_column=self.weights_column,
                            name=f"{self.name}[detected]")

    def lnpdraw(self) -> np.ndarray:
        out = np.zeros(len(self.df))
        for col in self.lnpdraw_columns:
            out = out + self.df[col].to_numpy(dtype=float)
        return out

    def weights(self) -> np.ndarray:
        if self.weights_column is None:
            return np.ones(len(self.df))
        return self.df[self.weights_column].to_numpy(dtype=float)

    def restrict_gps_window(self, gps_start: float, gps_end: float,
                            column: str = "time_geocenter") -> "SelectionSet":
        """Epoch-2 primary recipe (§5): restrict the joint O4ab set to a GPS window.

        NOTE: ``total_analysis_time``/``total_generated`` are NOT rescaled here,
        the locked recipe relies on the month-wise mixture ``weights`` partitioning
        the joint VT sum, so the restricted weighted sum *is* the window's share.
        """
        t = self.df[column].to_numpy(dtype=float)
        mask = (t >= gps_start) & (t < gps_end)
        return SelectionSet(df=self.df.loc[mask].reset_index(drop=True),
                            total_analysis_time=self.total_analysis_time,
                            total_generated=self.total_generated,
                            lnpdraw_columns=self.lnpdraw_columns,
                            weights_column=self.weights_column,
                            name=f"{self.name}[gps {gps_start}-{gps_end}]")


# ===================================================================== predictor

@dataclass
class PredictionResult:
    """Everything the scoring stage consumes, per epoch."""

    vt: np.ndarray                       # per-hyper-sample VT (estimator units)
    ln_sum_w: np.ndarray                 # per-sample ln sum of importance weights
    neff: np.ndarray                     # per-sample importance ESS
    w_bar: np.ndarray                    # hyper-averaged normalized weights over detected rows
    detected: SelectionSet               # the detected injection subset used
    mu: np.ndarray | None = None         # per-sample expected detected count (needs rate)
    flagged_fraction: float = 0.0        # ESS-rule failures (pre-reg §5)

    def predictive_cdf(self, column: str) -> tuple[np.ndarray, np.ndarray]:
        """Sorted values + cumulative hyper-averaged weights for observable ``column``.

        This is the frozen detected-population predictive CDF ``F_pred`` of
        pre-reg §1.1 (hyperposterior-averaged).
        """
        u = self.detected.df[column].to_numpy(dtype=float)
        order = np.argsort(u)
        cum = np.cumsum(self.w_bar[order])
        total = cum[-1]
        if not np.isfinite(total) or total <= 0:
            raise ValueError(f"degenerate predictive weights for {column}")
        return u[order], cum / total


class PosteriorPredictor:
    """Forward-model the frozen fit through an injection selection function."""

    def __init__(self, model: FrozenPopulationModel, selection: SelectionSet,
                 far_threshold_per_yr: float = 1.0, chunk_size: int = 64):
        self.model = model
        self.selection = selection
        self.far_threshold = far_threshold_per_yr
        self.chunk_size = int(chunk_size)

    def run(self, n_events_for_ess: int | None = None) -> PredictionResult:
        """Compute per-hyper-sample VT, ESS, and the averaged predictive weights.

        Memory-bounded: iterates hyperposterior samples in chunks; holds one
        float64 vector over detected injections per chunk row (8 GB RAM rule,
        pre-reg §9).
        """
        det = self.selection.detected(self.far_threshold)
        n_det = len(det.df)
        if n_det == 0:
            raise ValueError("no detected injections under the FAR cut")
        lnpdraw = det.lnpdraw()
        ln_mix_w = np.log(det.weights())
        S = self.model.n_samples

        vt = np.empty(S)
        ln_sum_w = np.empty(S)
        neff = np.empty(S)
        w_bar = np.zeros(n_det)

        for s in range(S):
            lnp = self.model.log_prob_injection_vars(det.df, s)
            lnw = ln_mix_w + lnp - lnpdraw
            m = np.max(lnw)
            if not np.isfinite(m):
                vt[s], ln_sum_w[s], neff[s] = 0.0, -np.inf, 0.0
                continue
            w = np.exp(lnw - m)
            sw = w.sum()
            ln_sum_w[s] = m + np.log(sw)
            # Locked VT estimator (§5), verbatim shape:
            vt[s] = (self.selection.total_analysis_time
                     * np.exp(ln_sum_w[s]) / self.selection.total_generated)
            neff[s] = sw ** 2 / np.square(w).sum()
            w_bar += w / sw
        w_bar /= S

        flagged = 0.0
        if n_events_for_ess is not None:
            flagged = float(np.mean(neff <= 4 * n_events_for_ess))
        return PredictionResult(vt=vt, ln_sum_w=ln_sum_w, neff=neff, w_bar=w_bar,
                                detected=det, flagged_fraction=flagged)

    # ------------------------------------------------------------------ H2 count
    def expected_counts(self, result: PredictionResult,
                        rate_samples: np.ndarray,
                        rate_to_vt_unit_factor: float = 1.0) -> np.ndarray:
        """``mu_s = R_s * VT_s`` (pre-reg §1 H2).

        ``rate_to_vt_unit_factor`` is the single pinned spot where the released
        rate posterior's units meet the injection release's VT units for the real
        epochs (resolved from the release READMEs at forward-model build and
        documented; pre-reg §2 'no third option'). Synthetic tests construct both
        sides consistently and use 1.0.
        """
        rate = np.asarray(rate_samples, dtype=float)
        if rate.shape != result.vt.shape:
            raise ValueError("rate_samples must align 1:1 with hyperposterior samples")
        mu = rate * result.vt * float(rate_to_vt_unit_factor)
        result.mu = mu
        return mu

    @staticmethod
    def count_pmf(mu: np.ndarray, n_max: int | None = None) -> np.ndarray:
        """``P(N) = (1/S) sum_s Poisson(N; mu_s)`` for N = 0..n_max (locked §1 H2)."""
        from scipy.stats import poisson

        mu = np.asarray(mu, dtype=float)
        if n_max is None:
            n_max = int(np.ceil(mu.max() + 10 * np.sqrt(mu.max() + 1) + 20))
        ns = np.arange(n_max + 1)
        pmf = poisson.pmf(ns[:, None], mu[None, :]).mean(axis=1)
        return pmf

    # ------------------------------------------------------------------ mocks (§7)
    def draw_mock_catalog(self, result: PredictionResult, rng: np.random.Generator,
                          n_events: int | None = None) -> tuple[pd.DataFrame, int]:
        """One mock catalog from the frozen model through the locked selection.

        Hyper-sample index is drawn per mock (hyper variance); the catalog size is
        Poisson(mu_s) when ``n_events`` is None (finite-catalog 'catalog variance'
        per arXiv:2603.00239 / pre-reg §1.3), requires ``expected_counts`` first.
        Event parameters are importance-resampled detected injections under Λ_s.
        """
        s = int(rng.integers(self.model.n_samples))
        det = result.detected
        lnp = self.model.log_prob_injection_vars(det.df, s)
        lnw = np.log(det.weights()) + lnp - det.lnpdraw()
        lnw -= lnw.max()
        w = np.exp(lnw)
        w /= w.sum()
        if n_events is None:
            if result.mu is None:
                raise ValueError("call expected_counts() before Poisson-size mocks")
            n = int(rng.poisson(result.mu[s]))
        else:
            n = int(n_events)
        n = max(n, 1)
        idx = rng.choice(len(det.df), size=n, replace=True, p=w)
        return det.df.iloc[idx].reset_index(drop=True), s


# ===================================================================== noise emulation

@dataclass
class NoiseEmulator:
    """Measurement-noise emulation for mock catalogs (pre-reg §7).

    Internal to null construction and validated solely by the mock-calibration
    gate. Emulates a per-event 'PE posterior' by scattering the true value and
    drawing samples around the scattered center, in a transformed space per
    observable: ``log`` (positive quantities like m1) or ``identity``.

    The PIT operator applied to these emulated posteriors is the identical
    operator applied to real events, so the null distribution of every statistic
    is exact by construction (pre-reg §1.1).
    """

    sigmas: dict[str, float] = field(default_factory=dict)   # observable -> rel/abs sigma
    transforms: dict[str, str] = field(default_factory=dict)  # observable -> log|identity
    n_samples: int = 128

    def emulate(self, truth_row: pd.Series, observables: tuple[str, ...],
                rng: np.random.Generator) -> dict[str, np.ndarray]:
        out: dict[str, np.ndarray] = {}
        for obs in observables:
            x = float(truth_row[obs])
            sig = self.sigmas.get(obs, 0.1)
            tf = self.transforms.get(obs, "log")
            if tf == "log":
                center = np.log(max(x, 1e-12)) + rng.normal(0, sig)
                out[obs] = np.exp(center + rng.normal(0, sig, size=self.n_samples))
            else:
                center = x + rng.normal(0, sig)
                out[obs] = center + rng.normal(0, sig, size=self.n_samples)
        return out
