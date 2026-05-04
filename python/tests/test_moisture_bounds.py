"""Tests for post-diffusion soil moisture bounds enforcement."""

import logging
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest

from utahlsm.core import UtahLSM
from utahlsm.exceptions import SolverError
from utahlsm.physics.soil.soil import Soil


def _make_model() -> UtahLSM:
    model = UtahLSM.__new__(UtahLSM)  # type: ignore[attr-defined]
    model.logger = logging.getLogger("test")
    model.input = SimpleNamespace(  # type: ignore[assignment]
        numerics=SimpleNamespace(
            tolerances=SimpleNamespace(moisture_bounds=1e-6)
        )
    )
    model.sfc_state = SimpleNamespace(moisture=np.array([0.20], dtype=float))  # type: ignore[assignment]
    model.soil_state = SimpleNamespace(  # type: ignore[assignment]
        moisture=np.array([[0.20], [0.25], [0.30]], dtype=float)
    )
    fake_soil = SimpleNamespace(
        properties=SimpleNamespace(
            residual=np.array([0.10, 0.10, 0.10], dtype=float),
            porosity=np.array([0.40, 0.40, 0.40], dtype=float),
        ),
        expand_profile_property=(
            lambda prop, soil_q: prop[:, None]  # type: ignore[arg-type]
            if np.asarray(soil_q).ndim == 2 else prop  # type: ignore[arg-type]
        ),
        diffusivity_moisture=lambda soil_q: soil_q,  # type: ignore[arg-type]
        conductivity_gradient=lambda soil_q: soil_q,  # type: ignore[arg-type]
    )
    fake_soil.logger = logging.getLogger("test")
    fake_soil.enforce_moisture_bounds = (
        lambda moisture, tol=1e-8: Soil.enforce_moisture_bounds(  # type: ignore[arg-type]
            fake_soil, moisture, tol  # type: ignore[arg-type]
        )
    )
    model.soil = fake_soil  # type: ignore[assignment]
    return model


def test_solve_diffusion_mois_clips_tiny_overshoots():
    model = _make_model()

    def fake_solve_mixed_moisture(**_kwargs: Any) -> None:
        model.soil_state.moisture[:] = np.array(
            [[0.10 - 5e-7], [0.25], [0.40 + 5e-7]],
            dtype=float,
        )

    model._solve_mixed_moisture = fake_solve_mixed_moisture  # type: ignore[assignment]

    model._solve_diffusion_mois()  # type: ignore[attr-defined]

    assert np.allclose(
        model.soil_state.moisture,
        np.array([[0.10], [0.25], [0.40]], dtype=float),
    )


def test_solve_diffusion_mois_raises_on_material_overshoot():
    model = _make_model()

    def fake_solve_mixed_moisture(**_kwargs: Any) -> None:
        model.soil_state.moisture[:] = np.array(
            [[0.10 - 2e-4], [0.25], [0.40 + 2e-4]],
            dtype=float,
        )

    model._solve_mixed_moisture = fake_solve_mixed_moisture  # type: ignore[assignment]

    with pytest.raises(SolverError, match="Soil moisture left physical bounds"):
        model._solve_diffusion_mois()  # type: ignore[attr-defined]
