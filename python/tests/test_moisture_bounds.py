"""Tests for post-diffusion soil moisture bounds enforcement."""

import logging
from types import SimpleNamespace

import numpy as np
import pytest

from utahlsm.core import UtahLSM
from utahlsm.exceptions import SolverError


def _make_model() -> UtahLSM:
    model = UtahLSM.__new__(UtahLSM)
    model.logger = logging.getLogger("test")
    model.input = SimpleNamespace(
        numerics=SimpleNamespace(
            tolerances=SimpleNamespace(smb_flux=1e-6)
        )
    )
    model.sfc_state = SimpleNamespace(moisture=np.array([0.20], dtype=float))
    model.soil_state = SimpleNamespace(
        moisture=np.array([[0.20], [0.25], [0.30]], dtype=float)
    )
    model.soil = SimpleNamespace(
        properties=SimpleNamespace(
            residual=np.array([0.10, 0.10, 0.10], dtype=float),
            porosity=np.array([0.40, 0.40, 0.40], dtype=float),
        ),
        _expand_profile_property=(
            lambda prop, soil_q: prop[:, None]
            if np.asarray(soil_q).ndim == 2 else prop
        ),
        diffusivity_moisture=lambda soil_q: soil_q,
        conductivity_gradient=lambda soil_q: soil_q,
    )
    return model


def test_solve_diffusion_mois_clips_tiny_overshoots():
    model = _make_model()

    def fake_solve_diffusion(**_kwargs):
        model.soil_state.moisture[:] = np.array(
            [[0.10 - 5e-7], [0.25], [0.40 + 5e-7]],
            dtype=float,
        )

    model._solve_diffusion = fake_solve_diffusion  # type: ignore[assignment]

    model._solve_diffusion_mois()

    assert np.allclose(
        model.soil_state.moisture,
        np.array([[0.10], [0.25], [0.40]], dtype=float),
    )


def test_solve_diffusion_mois_raises_on_material_overshoot():
    model = _make_model()

    def fake_solve_diffusion(**_kwargs):
        model.soil_state.moisture[:] = np.array(
            [[0.10 - 2e-4], [0.25], [0.40 + 2e-4]],
            dtype=float,
        )

    model._solve_diffusion = fake_solve_diffusion  # type: ignore[assignment]

    with pytest.raises(SolverError, match="Soil moisture left physical bounds"):
        model._solve_diffusion_mois()
