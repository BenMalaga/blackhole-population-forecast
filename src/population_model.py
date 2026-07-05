"""The frozen O3 population model (PRE_REGISTRATION.md §2-§3).

Loads the GWTC-3 Power-Law+Peak + iid-spin + power-law-redshift hyperposterior
*as released* (Zenodo 5655785, member
``analyses/PowerLawPeak/o1o2o3_mass_c_iid_mag_iid_tilt_powerlaw_redshift_result.json``)
and evaluates the population density at those hyperparameter samples using
``gwpopulation``'s implementations of the same models (pre-reg §3: we do not refit;
any naming/parameter-convention mismatch is a correctness question resolved against
the release README/tutorial, never a tunable choice).

Density convention (matches the sensitivity releases' ``lnpdraw`` measure)
--------------------------------------------------------------------------
``log_prob_injection_vars`` returns ``ln p(m1_src, m2_src, z, s1_vec, s2_vec | Λ)``:

* mass: gwpopulation's ``SinglePeakSmoothedMassDistribution`` gives ``p(m1, q)``;
  the Jacobian ``q -> m2`` contributes ``1/m1``.
* spins: the iid magnitude (Beta) x iid tilt (isotropic + truncated-Gaussian
  aligned mixture) density over ``(a, cos_tilt)`` converts to cartesian components
  via ``p_cart = p(a) * p(cos_tilt) / (2*pi*a^2)`` (azimuth uniform).
* redshift: gwpopulation's ``PowerLawRedshift`` normalized pdf on ``[0, z_max]``.

NO out-of-sample predictive score is computed in this module. It evaluates input
(training-epoch) data only; scoring lives in ``src/score.py`` behind the lock guard.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

#: The pinned primary frozen artifact (pre-reg §2 member path, from the release
#: tutorial). Stream observation (2026-06-10): the tarball ALSO carries a
#: byte-identical copy (sha256-verified,
#: ``03a759f07c39dbd967f1adb888c1e3e9e9763f836162b9d6058dad31cfe4836d``) under
#: ``analyses/PowerLawPeak/OtherModels/`` alongside the mass-model c/d/e/f
#: variants. The pinned top-level path governs; :func:`locate_frozen_artifact`
#: prefers it and tolerates only byte-identical duplicates elsewhere.
GWTC3_PLPEAK_RESULT = (
    ROOT / "data" / "gwtc3_pop" / "analyses" / "PowerLawPeak"
    / "o1o2o3_mass_c_iid_mag_iid_tilt_powerlaw_redshift_result.json"
)

GWTC3_PLPEAK_FILENAME = "o1o2o3_mass_c_iid_mag_iid_tilt_powerlaw_redshift_result.json"


def locate_frozen_artifact(filename: str = GWTC3_PLPEAK_FILENAME) -> Path:
    """Resolve the frozen result file among the streamed PowerLawPeak members.

    Preference order: the exact §2-pinned path; otherwise a *unique* exact-filename
    match elsewhere under the PowerLawPeak tree (excluding the GWTC2PaperPrior
    variant directory). Raises on zero matches or on ambiguity.
    """
    if filename == GWTC3_PLPEAK_FILENAME and GWTC3_PLPEAK_RESULT.exists():
        return GWTC3_PLPEAK_RESULT
    root = ROOT / "data" / "gwtc3_pop" / "analyses" / "PowerLawPeak"
    if not root.exists():
        raise FileNotFoundError(f"{root} missing, run fetch_data --tier gwtc3-plp")
    hits = sorted(p for p in root.rglob(filename)
                  if "GWTC2PaperPrior" not in p.parts)
    if len(hits) != 1:
        raise FileNotFoundError(
            f"expected exactly one {filename} under {root}, found {len(hits)}: "
            f"{[str(h) for h in hits]}"
        )
    return hits[0]

#: Canonical hyperparameter names (gwpopulation conventions) -> aliases that LVK
#: release files have used. Resolution against the release README/tutorial is a
#: correctness question (pre-reg §3); unknown columns raise rather than guess.
PARAM_ALIASES: dict[str, tuple[str, ...]] = {
    "alpha": ("alpha", "alpha_m"),
    "beta": ("beta", "beta_q"),
    "mmax": ("mmax", "m_max"),
    "mmin": ("mmin", "m_min"),
    "lam": ("lam", "lambda_peak", "lam_m"),
    "mpp": ("mpp", "mu_m", "mean_m"),
    "sigpp": ("sigpp", "sigma_m"),
    "delta_m": ("delta_m", "dm"),
    "amax": ("amax", "a_max"),
    "alpha_chi": ("alpha_chi",),
    "beta_chi": ("beta_chi",),
    "mu_chi": ("mu_chi",),
    # gwpopulation convention (verified against gwpopulation 1.3.1
    # ``convert_to_beta_parameters`` -> ``mu_var_max_to_alpha_beta_max`` and the
    # released posterior's value ranges): the parameter named ``sigma_chi`` is the
    # VARIANCE of the Beta magnitude distribution, not a standard deviation.
    "sigma_chi": ("sigma_chi", "sigma_chi_sq", "var_chi"),
    "xi_spin": ("xi_spin", "xi"),
    "sigma_spin": ("sigma_spin", "sigma_t", "sigma_1"),
    "lamb": ("lamb", "lambda_z", "kappa"),
    "rate": ("rate", "R"),
}

#: Parameters the density evaluation requires (after any mu/sigma -> alpha/beta
#: conversion). ``rate`` is carried separately for H2.
REQUIRED_PARAMS = (
    "alpha", "beta", "mmax", "mmin", "lam", "mpp", "sigpp", "delta_m",
    "amax", "alpha_chi", "beta_chi", "xi_spin", "sigma_spin", "lamb",
)


def _import_gwpop():
    """Import gwpopulation lazily with the numpy backend (no GPU, pre-reg §9)."""
    import gwpopulation

    try:  # gwpopulation >= 1.0 backend switch; numpy is also the default
        gwpopulation.set_backend("numpy")
    except Exception:  # pragma: no cover - older/newer API drift, default is numpy
        pass
    from gwpopulation.models.mass import SinglePeakSmoothedMassDistribution
    from gwpopulation.models.redshift import PowerLawRedshift
    from gwpopulation.models import spin as spin_models

    return SinglePeakSmoothedMassDistribution, PowerLawRedshift, spin_models


def mu_var_to_alpha_beta(mu: np.ndarray, var: np.ndarray, amax: np.ndarray):
    """Beta (mu, var, amax) -> (alpha, beta), via the locked library when available.

    Prefers ``gwpopulation.conversions.mu_var_max_to_alpha_beta_max`` (the same
    function the released analysis chain used, pre-reg §3); the closed form below
    is the documented fallback and matches it analytically for amax = 1.
    """
    try:
        from gwpopulation.conversions import mu_var_max_to_alpha_beta_max

        out = mu_var_max_to_alpha_beta_max(np.asarray(mu, dtype=float),
                                           np.asarray(var, dtype=float),
                                           np.asarray(amax, dtype=float))
        return np.asarray(out[0], dtype=float), np.asarray(out[1], dtype=float)
    except Exception:  # pragma: no cover - fallback if the import surface drifts
        mu = np.asarray(mu, dtype=float) / np.asarray(amax, dtype=float)
        var = np.asarray(var, dtype=float) / np.asarray(amax, dtype=float) ** 2
        nu = mu * (1 - mu) / var - 1
        return mu * nu, (1 - mu) * nu


def _resolve_columns(columns: list[str]) -> dict[str, str]:
    """Map canonical parameter names to the columns actually present."""
    found: dict[str, str] = {}
    for canonical, aliases in PARAM_ALIASES.items():
        for alias in aliases:
            if alias in columns:
                found[canonical] = alias
                break
    return found


@dataclass
class FrozenPopulationModel:
    """Power-Law+Peak + iid spins + power-law redshift at frozen hyperposterior samples.

    ``hyperposterior`` columns are canonical names (REQUIRED_PARAMS, plus optionally
    ``rate``). Use :meth:`from_result_json` for the released file or build the frame
    directly for synthetic tests.
    """

    hyperposterior: pd.DataFrame
    z_max: float = 2.3            # GWTC-3 analysis redshift bound (release tutorial)
    mmax_grid: float = 300.0      # mass normalization grid bound (gwpopulation default)
    source_path: str | None = None
    _models: dict = field(default_factory=dict, repr=False)

    def __post_init__(self) -> None:
        missing = [p for p in REQUIRED_PARAMS if p not in self.hyperposterior.columns]
        if missing:
            raise ValueError(
                f"hyperposterior is missing required parameters {missing}; "
                f"have {sorted(self.hyperposterior.columns)}"
            )
        if not len(self.hyperposterior):
            raise ValueError("hyperposterior has zero samples")

    # ------------------------------------------------------------------ loading
    @classmethod
    def from_result_json(cls, path: Path | str = GWTC3_PLPEAK_RESULT,
                         **kwargs) -> "FrozenPopulationModel":
        """Load a bilby result JSON (the released frozen artifact) as-is.

        Handles both bilby posterior layouts (``posterior.content`` dict-of-lists
        and a plain dict). All released samples are used, no subsampling here
        (pre-reg §2).
        """
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(
                f"{path} not on disk. The gwtc3-plp stream (src/fetch_data.py) has "
                f"not delivered the pinned member yet, re-run "
                f"`.venv/bin/python -m src.fetch_data --tier gwtc3-plp` when the "
                f"disk window opens (data/README.md)."
            )
        raw = json.loads(path.read_text())
        posterior = raw.get("posterior", raw)
        if isinstance(posterior, dict) and "content" in posterior:
            posterior = posterior["content"]
        frame = pd.DataFrame({k: np.asarray(v, dtype=float)
                              for k, v in posterior.items()
                              if np.asarray(v).ndim == 1})
        resolved = _resolve_columns(list(frame.columns))
        out = pd.DataFrame(index=frame.index)
        for canonical, col in resolved.items():
            out[canonical] = frame[col].to_numpy(dtype=float)
        # amax is fixed (not sampled) in the released iid-magnitude analysis
        if "amax" not in out.columns:
            out["amax"] = 1.0
        # mu/sigma -> alpha/beta if the release uses the (mu, var) parametrization.
        # ``sigma_chi`` is the VARIANCE in the gwpopulation/LVK convention (see
        # PARAM_ALIASES note; verified against gwpopulation 1.3.1's own
        # convert_to_beta_parameters and the released posterior ranges).
        if "alpha_chi" not in out.columns and "mu_chi" in out.columns:
            if "sigma_chi" not in out.columns:
                raise ValueError(
                    "release file has mu_chi but no sigma_chi variance column; "
                    f"columns: {sorted(frame.columns)}"
                )
            a, b = mu_var_to_alpha_beta(out["mu_chi"].to_numpy(),
                                        out["sigma_chi"].to_numpy(),
                                        out["amax"].to_numpy())
            out["alpha_chi"], out["beta_chi"] = a, b
        missing = [p for p in REQUIRED_PARAMS if p not in out.columns]
        if missing:
            raise ValueError(
                "could not resolve release hyperparameter names -> "
                f"missing {missing}. Release columns: {sorted(frame.columns)}. "
                "Extend PARAM_ALIASES against the release README/tutorial "
                "(correctness fix per PRE_REGISTRATION.md §3, document in repo)."
            )
        return cls(hyperposterior=out.reset_index(drop=True),
                   source_path=str(path), **kwargs)

    # ------------------------------------------------------------------ pieces
    @property
    def n_samples(self) -> int:
        return len(self.hyperposterior)

    @property
    def rate_samples(self) -> np.ndarray:
        """Merger-rate posterior, as released (pre-reg §2 locked rule).

        Raises if the released file carries no rate column, reconstruction then
        follows the release recipe via a labeled amendment, never silently.
        """
        if "rate" not in self.hyperposterior.columns:
            raise KeyError(
                "frozen hyperposterior carries no 'rate' column. Pre-reg §2: "
                "reconstruct via the release README recipe + labeled amendment "
                "BEFORE any H2 evaluation."
            )
        return self.hyperposterior["rate"].to_numpy(dtype=float)

    def _mass_model(self):
        if "mass" not in self._models:
            SinglePeak, _, _ = _import_gwpop()
            self._models["mass"] = SinglePeak(mmin=1.0, mmax=self.mmax_grid)
        return self._models["mass"]

    def _redshift_model(self):
        if "redshift" not in self._models:
            _, PowerLawRedshift, _ = _import_gwpop()
            self._models["redshift"] = PowerLawRedshift(z_max=self.z_max)
        return self._models["redshift"]

    # ------------------------------------------------------------------ density
    @staticmethod
    def _spins_to_polar(df: pd.DataFrame) -> tuple[np.ndarray, ...]:
        """(a_1, a_2, cos_tilt_1, cos_tilt_2) from cartesian or polar columns.

        Cartesian components are the locked convention for the injection mixture
        (pre-reg §10); PE parquets carry polar columns directly.
        """
        if {"spin1x", "spin1y", "spin1z", "spin2x", "spin2y", "spin2z"} <= set(df.columns):
            a1 = np.sqrt(df["spin1x"] ** 2 + df["spin1y"] ** 2 + df["spin1z"] ** 2)
            a2 = np.sqrt(df["spin2x"] ** 2 + df["spin2y"] ** 2 + df["spin2z"] ** 2)
            with np.errstate(divide="ignore", invalid="ignore"):
                c1 = np.where(a1 > 0, df["spin1z"] / a1, 0.0)
                c2 = np.where(a2 > 0, df["spin2z"] / a2, 0.0)
            return (np.asarray(a1), np.asarray(a2),
                    np.asarray(c1), np.asarray(c2))
        if {"a_1", "a_2", "cos_tilt_1", "cos_tilt_2"} <= set(df.columns):
            return (df["a_1"].to_numpy(dtype=float),
                    df["a_2"].to_numpy(dtype=float),
                    df["cos_tilt_1"].to_numpy(dtype=float),
                    df["cos_tilt_2"].to_numpy(dtype=float))
        raise KeyError("need cartesian (spin1x..spin2z) or polar (a_1..cos_tilt_2) spins")

    def _dataset(self, df: pd.DataFrame) -> dict[str, np.ndarray]:
        if {"mass1_source", "mass2_source"} <= set(df.columns):
            m1 = df["mass1_source"].to_numpy(dtype=float)
            m2 = df["mass2_source"].to_numpy(dtype=float)
        elif {"mass_1_source", "mass_ratio"} <= set(df.columns):
            m1 = df["mass_1_source"].to_numpy(dtype=float)
            m2 = m1 * df["mass_ratio"].to_numpy(dtype=float)
        else:
            raise KeyError("need (mass1_source, mass2_source) or "
                           "(mass_1_source, mass_ratio) columns")
        a1, a2, c1, c2 = self._spins_to_polar(df)
        return {
            "mass_1": m1,
            "mass_ratio": m2 / m1,
            "a_1": a1, "a_2": a2,
            "cos_tilt_1": c1, "cos_tilt_2": c2,
            "redshift": df["redshift"].to_numpy(dtype=float),
        }

    def log_prob_injection_vars(self, df: pd.DataFrame, sample_index: int,
                                spin_measure: str = "auto") -> np.ndarray:
        """``ln p(m1_src, m2_src, z, spins | Λ_s)`` in the injection (lnpdraw) measure.

        ``spin_measure``:
          * ``"cartesian"``: density over the six cartesian components
            (matches the injection sets' lnpdraw measure): includes the
            ``1/(2*pi*a^2)`` factors per spin.
          * ``"polar"``: density over ``(a_1, a_2, cos_tilt_1, cos_tilt_2)``
            (PE-sample measure; azimuths marginalized).
          * ``"auto"``: cartesian iff cartesian columns are present.
        """
        params = self.hyperposterior.iloc[int(sample_index)].to_dict()
        data = self._dataset(df)
        if spin_measure == "auto":
            spin_measure = ("cartesian"
                            if "spin1x" in df.columns else "polar")

        _, _, spin_models = _import_gwpop()
        mass_model = self._mass_model()
        z_model = self._redshift_model()

        p_mass = mass_model(
            {"mass_1": data["mass_1"], "mass_ratio": data["mass_ratio"]},
            alpha=params["alpha"], beta=params["beta"], mmin=params["mmin"],
            mmax=params["mmax"], lam=params["lam"], mpp=params["mpp"],
            sigpp=params["sigpp"], delta_m=params["delta_m"],
        )
        p_spin = spin_models.iid_spin(
            {"a_1": data["a_1"], "a_2": data["a_2"],
             "cos_tilt_1": data["cos_tilt_1"], "cos_tilt_2": data["cos_tilt_2"]},
            xi_spin=params["xi_spin"], sigma_spin=params["sigma_spin"],
            amax=params["amax"], alpha_chi=params["alpha_chi"],
            beta_chi=params["beta_chi"],
        )
        p_z = z_model({"redshift": data["redshift"]}, lamb=params["lamb"])

        with np.errstate(divide="ignore", invalid="ignore"):
            lnp = (np.log(np.asarray(p_mass, dtype=float))
                   + np.log(np.asarray(p_spin, dtype=float))
                   + np.log(np.asarray(p_z, dtype=float))
                   - np.log(data["mass_1"]))           # Jacobian q -> m2
            if spin_measure == "cartesian":
                lnp -= np.log(2 * np.pi * data["a_1"] ** 2)
                lnp -= np.log(2 * np.pi * data["a_2"] ** 2)
        return np.where(np.isfinite(lnp), lnp, -np.inf)

    def log_prob_marginal(self, df: pd.DataFrame, sample_index: int,
                          observables: tuple[str, ...] = ("m1", "z")) -> np.ndarray:
        """ln of the (m1, z) population density (mass1 marginal x redshift), for the
        surprise index (pre-reg §1.6). Marginalizes q on a grid; spins integrate to 1.
        """
        params = self.hyperposterior.iloc[int(sample_index)].to_dict()
        if tuple(observables) != ("m1", "z"):
            raise NotImplementedError("surprise density is defined on (m1, z)")
        m1 = (df["mass1_source"] if "mass1_source" in df.columns
              else df["mass_1_source"]).to_numpy(dtype=float)
        z = df["redshift"].to_numpy(dtype=float)
        mass_model = self._mass_model()
        q_grid = np.linspace(1e-3, 1.0, 256)
        m1_g, q_g = np.meshgrid(m1, q_grid, indexing="ij")
        p_joint = mass_model(
            {"mass_1": m1_g.ravel(), "mass_ratio": q_g.ravel()},
            alpha=params["alpha"], beta=params["beta"], mmin=params["mmin"],
            mmax=params["mmax"], lam=params["lam"], mpp=params["mpp"],
            sigpp=params["sigpp"], delta_m=params["delta_m"],
        )
        p_m1 = np.trapezoid(np.asarray(p_joint, dtype=float).reshape(m1_g.shape),
                            q_grid, axis=1)
        p_z = np.asarray(self._redshift_model()({"redshift": z},
                                                lamb=params["lamb"]), dtype=float)
        with np.errstate(divide="ignore"):
            lnp = np.log(p_m1) + np.log(p_z)
        return np.where(np.isfinite(lnp), lnp, -np.inf)
