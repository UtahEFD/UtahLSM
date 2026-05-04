"""Regression tests for diffusion boundary conditions."""

import logging
from types import SimpleNamespace
from typing import Any

import numpy as np

from utahlsm.core import UtahLSM
from utahlsm.data_models import (
    IterationsConfig,
    NumericsConfig,
    SoilState,
    SolverState,
    TolerancesConfig,
)


def test_bottom_neumann_boundary_couples_deepest_heat_layer() -> None:
    """The bottom zero-gradient row must not decouple the deepest layer."""
    model: Any = UtahLSM.__new__(UtahLSM)
    model.logger = logging.getLogger("test")
    model.tstep = 1.0
    model.ncol = 1
    model.input = SimpleNamespace(
        grid=SimpleNamespace(
            nz=3, nx=1, ny=1, z=np.array([0.0, -1.0, -2.0])
        ),
        numerics=NumericsConfig(
            heat_diffusion_back_weight=1.0,
            iterations=IterationsConfig(
                sfc_flux=1,
                seb_bracket=1,
                seb_root=1,
                smb_flux=1,
                moisture_picard=1,
                coupling=1,
            ),
            tolerances=TolerancesConfig(
                sfc_flux=1e-6,
                seb_root=1e-12,
                smb_flux=1e-12,
                moisture_picard=1e-8,
                moisture_bounds=1e-12,
                coupling_temp=1e-6,
                coupling_mois=1e-12,
            ),
        ),
    )
    model.solver_state = SolverState()
    model.soil_state = SoilState(
        temperature=np.array([[0.0], [0.0], [10.0]]),
        moisture=np.full((3, 1), 0.25),
        type=np.array(["clay", "clay", "clay"], dtype=object),
    )

    model._solve_diffusion(
        state_field=model.soil_state.temperature,
        get_diffusivity=lambda moisture: np.ones_like(moisture, dtype=float),
        sfc_boundary=np.array([0.0]),
        field_name="temperature",
    )

    expected = np.array([[0.0], [10.0 / 7.0], [30.0 / 7.0]])
    np.testing.assert_allclose(model.soil_state.temperature, expected)
