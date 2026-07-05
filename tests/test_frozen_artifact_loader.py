"""Loader smoke test against the REAL frozen GWTC-3 artifact (input data only).

PRE_REGISTRATION.md §2 classifies the GWTC-3 hyperposterior as *training-epoch
input data*, acquiring and reading it carries no firewall implications. This test
loads it and checks structure/convention resolution. It computes NO predictive
quantity: no selection function, no VT, no PIT, no count, nothing involving
GWTC-4/5 data. Skipped cleanly while the gwtc3-plp stream has not yet delivered
the member.
"""

from __future__ import annotations

import numpy as np
import pytest

from src.population_model import (REQUIRED_PARAMS, FrozenPopulationModel,
                                  locate_frozen_artifact)


def _artifact_or_skip():
    try:
        return locate_frozen_artifact()
    except FileNotFoundError as exc:
        pytest.skip(f"frozen artifact not streamed yet: {exc}")


def test_locate_and_load_frozen_hyperposterior():
    path = _artifact_or_skip()
    model = FrozenPopulationModel.from_result_json(path)
    # the released analysis: all samples used, no subsampling (pre-reg §2)
    assert model.n_samples > 1000
    for p in REQUIRED_PARAMS:
        assert p in model.hyperposterior.columns, p
        vals = model.hyperposterior[p].to_numpy()
        assert np.all(np.isfinite(vals)), f"non-finite {p}"
    # Beta conversion must yield a proper (positive-parameter) Beta distribution
    assert (model.hyperposterior["alpha_chi"] > 0).all()
    assert (model.hyperposterior["beta_chi"] > 0).all()
    # the released rate posterior is present and positive (locked H2 rule: as-is)
    rate = model.rate_samples
    assert (rate > 0).all()


def test_frozen_density_evaluates_on_synthetic_points():
    """Evaluate the *frozen* density at a handful of synthetic parameter points.

    Synthetic inputs, frozen-model output, still zero contact with GWTC-4/5 data.
    """
    import pandas as pd

    path = _artifact_or_skip()
    model = FrozenPopulationModel.from_result_json(path)
    df = pd.DataFrame({
        "mass1_source": [10.0, 35.0, 60.0],
        "mass2_source": [8.0, 30.0, 30.0],
        "redshift": [0.2, 0.4, 0.8],
        "a_1": [0.2, 0.4, 0.7], "a_2": [0.1, 0.3, 0.5],
        "cos_tilt_1": [0.9, 0.5, -0.2], "cos_tilt_2": [0.8, 0.1, 0.3],
    })
    lnp = model.log_prob_injection_vars(df, sample_index=0, spin_measure="polar")
    assert lnp.shape == (3,)
    # in-support points must be finite; a 35 Msun primary near the peak should not
    # be wildly improbable under the frozen fit
    assert np.all(np.isfinite(lnp))
