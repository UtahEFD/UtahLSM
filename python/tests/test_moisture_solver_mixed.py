"""Smoke tests for the mixed-form soil moisture solver."""

import logging
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest

from utahlsm.core import UtahLSM
from utahlsm.data_models import (
    IterationsConfig,
    NumericsConfig,
    SolverState,
    TolerancesConfig,
)
from utahlsm.physics.soil.soil_brookscorey import BrooksCorey
from utahlsm.util.io.soil_properties_loader import SoilPropertiesLoader


def _make_model(theta: float = 0.25) -> Any:
    model: Any = UtahLSM.__new__(UtahLSM)
    model.logger = logging.getLogger("test")
    model.ncol = 1
    model.tstep = 3600.0
    model.input = SimpleNamespace(
        grid=SimpleNamespace(nx=1, ny=1, nz=3, z=np.array([0.05, 0.20, 0.35])),
        numerics=NumericsConfig(
            heat_diffusion_back_weight=1.0,
            iterations=IterationsConfig(
                sfc_flux=10,
                seb_bracket=10,
                seb_root=10,
                smb_flux=50,
                moisture_picard=50,
                coupling=10,
            ),
            tolerances=TolerancesConfig(
                sfc_flux=1e-3,
                seb_root=1e-6,
                smb_flux=1e-6,
                moisture_picard=1e-8,
                moisture_bounds=1e-6,
                coupling_temp=1e-3,
                coupling_mois=1e-8,
            ),
        ),
    )
    model.solver_state = SolverState()
    model.sfc_state = SimpleNamespace(moisture=np.array([theta], dtype=float))
    model.soil_state = SimpleNamespace(
        moisture=np.array([theta, theta, theta], dtype=float)
    )
    props = SoilPropertiesLoader.load("cosby")
    model.soil = BrooksCorey(props, ["loam", "loam", "loam"], "cosby")
    return model


def test_mixed_moisture_solver_preserves_uniform_equilibrium() -> None:
    model = _make_model(theta=0.25)
    initial = np.array(model.soil_state.moisture, copy=True)

    # Testing internal solver behavior directly
    model._solve_diffusion_mois()

    assert np.allclose(model.soil_state.moisture, initial, atol=1e-8)


def test_mixed_moisture_solver_applies_sink_term() -> None:
    model = _make_model(theta=0.25)
    source = np.array([0.0, -1.0e-7, -1.0e-7], dtype=float)
    initial = np.array(model.soil_state.moisture, copy=True)

    # Testing internal solver behavior directly
    model._solve_mixed_moisture(source_term=source)

    assert np.all(model.soil_state.moisture[1:] < initial[1:])


def test_mixed_moisture_solver_conserves_prescribed_surface_flux() -> None:
    """Richards uses the SMB flux without recomputing the upper face."""
    model = _make_model(theta=0.25)
    initial = np.array(model.soil_state.moisture, copy=True)
    top_mass_flux = 3.0e-4
    model._matrix_top_flux = np.array([top_mass_flux])

    model._solve_diffusion_mois()

    final = np.asarray(model.soil_state.moisture, dtype=float)
    dz = abs(model.input.grid.z[1] - model.input.grid.z[0])
    storage_rate = dz * np.sum(final[1:] - initial[1:]) / model.tstep
    bottom_mass_flux = float(model._bottom_drainage[0])
    expected_rate = (top_mass_flux - bottom_mass_flux) / 1000.0
    assert storage_rate == pytest.approx(expected_rate, abs=1e-10)
    assert bottom_mass_flux > 0.0


def test_face_conductivity_uses_darcy_direction_at_texture_break() -> None:
    """Texture jumps must not infer flow direction from θ ordering."""
    K_upper = np.array([2.0e-6])
    K_lower = np.array([1.0e-8])

    downward = UtahLSM._moisture_face_conductivity(
        K_upper,
        K_lower,
        psi_upper=np.array([-0.50]),
        psi_lower=np.array([-0.55]),
        dz=0.10,
    )
    upward = UtahLSM._moisture_face_conductivity(
        K_upper,
        K_lower,
        psi_upper=np.array([-0.70]),
        psi_lower=np.array([-0.50]),
        dz=0.10,
    )

    assert downward[0] == K_upper[0]
    assert upward[0] == np.sqrt(K_upper[0] * K_lower[0])
