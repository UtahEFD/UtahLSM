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
    model.sfc_state = SimpleNamespace(soil_top_temperature=np.array([0.0]))
    model.soil = SimpleNamespace(
        conductivity_thermal=lambda moisture: np.ones_like(moisture),
        heat_capacity=lambda moisture: np.ones_like(moisture),
    )

    model._solve_soil_heat()

    expected = np.array([[0.0], [2.0], [6.0]])
    np.testing.assert_allclose(model.soil_state.temperature, expected)


def test_soil_heat_conserves_energy_across_property_jumps() -> None:
    """One face flux must close energy across discontinuous lambda and C."""
    model: Any = UtahLSM.__new__(UtahLSM)
    model.logger = logging.getLogger("test")
    model.tstep = 600.0
    model.ncol = 1
    theta = 0.5
    dz = 0.1
    model.input = SimpleNamespace(
        grid=SimpleNamespace(
            nz=5, nx=1, ny=1, z=-dz * np.arange(5, dtype=float)
        ),
        numerics=NumericsConfig(
            heat_diffusion_back_weight=theta,
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
    temperature_old = np.array([[300.0], [295.0], [289.0], [292.0], [285.0]])
    capacity = np.array([[1.0e6], [1.5e6], [3.0e6], [1.2e6], [4.0e6]])
    conductivity = np.array([[0.4], [1.6], [0.2], [2.4], [0.6]])
    model.soil_state = SoilState(
        temperature=temperature_old.copy(),
        moisture=np.full((5, 1), 0.25),
        type=np.array(["clay"] * 5, dtype=object),
    )
    surface_new = np.array([301.0])
    model.sfc_state = SimpleNamespace(soil_top_temperature=surface_new)
    model.soil = SimpleNamespace(
        conductivity_thermal=lambda moisture: conductivity,
        heat_capacity=lambda moisture: capacity,
    )
    source = np.array([[0.0], [1.0e-4], [-2.0e-5], [5.0e-5], [-1.0e-5]])

    model._solve_soil_heat(source_term=source)

    temperature_new = np.asarray(model.soil_state.temperature)
    energy_change = float(
        np.sum(capacity[1:] * dz * (temperature_new[1:] - temperature_old[1:]))
    )
    lambda_top = float(
        2.0 * conductivity[0, 0] * conductivity[1, 0]
        / (conductivity[0, 0] + conductivity[1, 0])
    )
    flux_old = lambda_top * (
        temperature_old[0, 0] - temperature_old[1, 0]
    ) / dz
    flux_new = lambda_top * (
        surface_new[0] - temperature_new[1, 0]
    ) / dz
    boundary_energy = model.tstep * (
        (1.0 - theta) * flux_old + theta * flux_new
    )
    source_energy = float(
        model.tstep * np.sum(capacity[1:] * dz * source[1:])
    )
    np.testing.assert_allclose(
        energy_change, boundary_energy + source_energy, atol=1e-7
    )
