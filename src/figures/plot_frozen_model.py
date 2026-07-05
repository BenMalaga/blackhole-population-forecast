"""Render the locked O3-era primary-mass model family (the design's INPUT model).

This figure shows the Power-Law+Peak primary-mass density that the GWTC-3
population fit defines, the model that this project freezes and then tests
out-of-sample against later catalogs. It plots ONLY the input model: a handful
of hyperparameter draws around the published GWTC-3 Power-Law+Peak values, plus
the median curve. No new-catalog data is read and no out-of-sample score is
computed here.

Hyperparameter values are the published GWTC-3 Power-Law+Peak maximum-a-posteriori
/ median estimates (Abbott et al. 2023, "Population of Merging Compact Binaries
Inferred Using Gravitational Waves through GWTC-3", Phys. Rev. X 13, 011048;
LVK GWTC-3 population data release, Zenodo 5655785). They are inputs to this
project's design, used here to draw the model family, not results.

Run:
    .venv/bin/python -m src.figures.plot_frozen_model
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless, deterministic
import matplotlib.pyplot as plt
import numpy as np

from gwpopulation.models.mass import SinglePeakSmoothedMassDistribution

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / "docs" / "figures" / "frozen_o3_mass_model.png"

SEED = 20260610  # the project's locked seed (PRE_REGISTRATION.md)

# Published GWTC-3 Power-Law+Peak primary-mass hyperparameters (medians) and
# approximate 1-sigma widths, used only to draw a few model-family curves.
# Source: Abbott et al. 2023 (GWTC-3 population), Table-level summaries; the
# released hyperposterior (Zenodo 5655785) is the canonical object this project
# loads at run time. These numbers are design INPUTS, not project results.
GWTC3_PLPEAK_MEDIAN = dict(
    alpha=3.5,      # power-law slope of the primary mass
    mmin=5.0,       # minimum mass (Msun)
    mmax=88.0,      # maximum mass (Msun)
    lam=0.038,      # fraction of systems in the Gaussian peak
    mpp=34.0,       # peak location (Msun)
    sigpp=3.6,      # peak width (Msun)
    delta_m=4.9,    # low-mass smoothing scale (Msun)
    beta=1.1,       # mass-ratio slope (unused for the m1 marginal)
)
# 1-sigma-ish spreads for the illustrative parameter draws (qualitative widths
# consistent with the published credible intervals).
GWTC3_PLPEAK_SPREAD = dict(
    alpha=0.6, mmin=1.3, mmax=10.0, lam=0.02,
    mpp=2.5, sigpp=1.5, delta_m=2.0,
)

COLORBLIND_DRAWS = "#9ecae1"   # light blue for the ensemble draws
MEDIAN_COLOR = "#08519c"       # deep blue for the median curve
PEAK_COLOR = "#d95f02"         # orange accent for the ~34 Msun peak marker


def primary_mass_pdf(m1: np.ndarray, params: dict) -> np.ndarray:
    """p(m1) for the Power-Law+Peak model, marginalized over mass ratio.

    Uses the project's own mass-model library (gwpopulation) so the curve is the
    exact functional form the analysis evaluates. We integrate the joint
    p(m1, q) over q to get the primary-mass marginal.
    """
    model = SinglePeakSmoothedMassDistribution(mmin=1.0, mmax=300.0)
    q_grid = np.linspace(1e-3, 1.0, 400)
    m1_g, q_g = np.meshgrid(m1, q_grid, indexing="ij")
    p_joint = np.asarray(
        model(
            {"mass_1": m1_g.ravel(), "mass_ratio": q_g.ravel()},
            alpha=params["alpha"], beta=params["beta"],
            mmin=params["mmin"], mmax=params["mmax"],
            lam=params["lam"], mpp=params["mpp"],
            sigpp=params["sigpp"], delta_m=params["delta_m"],
        ),
        dtype=float,
    ).reshape(m1_g.shape)
    return np.trapezoid(p_joint, q_grid, axis=1)


def draw_params(rng: np.random.Generator) -> dict:
    """One illustrative hyperparameter draw around the published medians."""
    p = dict(GWTC3_PLPEAK_MEDIAN)
    for key, sd in GWTC3_PLPEAK_SPREAD.items():
        p[key] = float(rng.normal(GWTC3_PLPEAK_MEDIAN[key], sd))
    # keep draws physical
    p["alpha"] = max(p["alpha"], 0.5)
    p["mmin"] = float(np.clip(p["mmin"], 3.0, 8.0))
    p["mmax"] = float(np.clip(p["mmax"], 60.0, 120.0))
    p["lam"] = float(np.clip(p["lam"], 0.0, 0.2))
    p["sigpp"] = max(p["sigpp"], 1.0)
    p["delta_m"] = float(np.clip(p["delta_m"], 1.0, 10.0))
    return p


def main() -> None:
    rng = np.random.default_rng(SEED)
    m1 = np.linspace(3.0, 100.0, 600)  # coarse-ish grid, renders in well under a second

    n_draws = 24
    draws = [primary_mass_pdf(m1, draw_params(rng)) for _ in range(n_draws)]
    median_curve = primary_mass_pdf(m1, GWTC3_PLPEAK_MEDIAN)

    plt.rcParams.update({"font.size": 11, "axes.linewidth": 0.8})
    fig, ax = plt.subplots(figsize=(8.0, 5.0), dpi=150)

    for i, d in enumerate(draws):
        ax.plot(
            m1, d, color=COLORBLIND_DRAWS, lw=1.0, alpha=0.55,
            zorder=1, label="hyperposterior draws" if i == 0 else None,
        )
    ax.plot(
        m1, median_curve, color=MEDIAN_COLOR, lw=2.6, zorder=3,
        label="published GWTC-3 median fit",
    )

    # Mark the ~34 Msun Gaussian peak the model encodes.
    ax.axvline(GWTC3_PLPEAK_MEDIAN["mpp"], color=PEAK_COLOR, lw=1.2,
               ls="--", alpha=0.8, zorder=2)
    ax.annotate(
        f"~{GWTC3_PLPEAK_MEDIAN['mpp']:.0f} M$_\\odot$ peak",
        xy=(GWTC3_PLPEAK_MEDIAN["mpp"], np.interp(GWTC3_PLPEAK_MEDIAN["mpp"], m1, median_curve)),
        xytext=(46, np.interp(GWTC3_PLPEAK_MEDIAN["mpp"], m1, median_curve) * 2.4),
        color=PEAK_COLOR, fontsize=10,
        arrowprops=dict(arrowstyle="->", color=PEAK_COLOR, lw=1.0),
    )

    ax.set_yscale("log")
    ax.set_xlim(3, 100)
    ax.set_ylim(1e-5, 5e-1)
    ax.set_xlabel("Primary black-hole mass  $m_1$  (M$_\\odot$)")
    ax.set_ylabel("Population density  $p(m_1)$")
    ax.set_title(
        "The locked O3-era model we test out-of-sample\n"
        "Power-Law+Peak primary-mass distribution (frozen GWTC-3 fit)",
        fontsize=12.5,
    )
    ax.grid(True, which="major", ls=":", lw=0.5, alpha=0.4)
    ax.legend(frameon=False, loc="upper right", fontsize=10)

    fig.text(
        0.012, 0.012,
        "Input model only: the design freezes this O3-era family, then scores later catalogs against it. "
        "No new-catalog data shown.",
        fontsize=7.5, color="#555555",
    )
    fig.tight_layout(rect=(0, 0.03, 1, 1))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=150, bbox_inches="tight")
    plt.close(fig)
    size_kb = OUT.stat().st_size / 1024
    print(f"wrote {OUT} ({size_kb:.0f} KB)")


if __name__ == "__main__":
    main()
