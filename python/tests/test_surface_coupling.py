"""Tests for the coupled surface SEB/SMB iteration.

These tests validate that the coupled surface solver leaves the diagnostic
fluxes consistent with the final converged surface state.
"""

import logging
from types import SimpleNamespace
from typing import Any

import numpy as np

from utahlsm.core import UtahLSM
from utahlsm.data_models import (
    AtmosphericState,
    IterationsConfig,
    NumericsConfig,
    SoilState,
    SurfaceState,
    TolerancesConfig,
)
from utahlsm.exceptions import SolverError


def _make_minimal_model(
    *,
    coupling_iterations: int = 3,
    tol_temp: float = 1e-6,
    tol_mois: float = 1e-12,
) -> Any:
    """Create a minimal UtahLSM instance for unit-testing private methods."""
    model: Any = UtahLSM.__new__(UtahLSM)
    model.logger = logging.getLogger("test")
    model.input = SimpleNamespace(
        numerics=NumericsConfig(
            heat_diffusion_back_weight=0.5,
            iterations=IterationsConfig(
                sfc_flux=1,
                seb_bracket=1,
                seb_root=1,
                smb_flux=1,
                moisture_picard=1,
                coupling=coupling_iterations,
            ),
            tolerances=TolerancesConfig(
                sfc_flux=1e-6,
                seb_root=1e-12,
                smb_flux=1e-12,
                moisture_picard=1e-8,
                moisture_bounds=1e-12,
                coupling_temp=tol_temp,
                coupling_mois=tol_mois,
            ),
        )
    )
    model.sfc_state = SurfaceState(temperature=300.0, moisture=0.2)
    return model


def test_compute_seb_vec_is_deterministic_and_pure() -> None:
    """SEB residual evaluation should be deterministic and not mutate state."""
    model = _make_minimal_model()
    model.ncol = 1

    # Provide minimal atmosphere and state containers as arrays.
    model.atm_state = AtmosphericState(
        wind_speed=np.array([3.0]),
        temperature=np.array([290.0]),
        specific_humidity=np.array([0.01]),
        pressure=np.array([101000.0]),
        radiation_net=np.array([100.0]),
        seb_storage=np.array([0.0]),
    )
    model.soil_state = SoilState(
        temperature=np.array([[280.0], [285.0]]),
        moisture=np.array([[0.25], [0.25]]),
        type=np.array(["clay", "clay"], dtype=object),
    )
    model.sfc_state = SurfaceState(
        temperature=np.array([280.0]),
        moisture=np.array([0.25]),
    )
    model.sfc_state.fluxes.ground_heat[0] = 1.0
    model.sfc_state.fluxes.sensible_heat[0] = 2.0
    model.sfc_state.fluxes.latent_heat[0] = 3.0
    model.sfc_state.fluxes.kinematic_heat[0] = 4.0
    model.sfc_state.fluxes.kinematic_moisture[0] = 5.0
    model.sfc_state.turbulence.obukhov_length[0] = 7.0
    model.sfc_state.turbulence.friction_velocity[0] = 8.0

    model.solver_state = SimpleNamespace(conductivity_thermal_mid=np.array([1.0]))
    model.soil = SimpleNamespace(surface_specific_humidity=lambda _T, _q, _p: np.zeros_like(_T))
    model.sfc = SimpleNamespace(
        fm=lambda _z1, _z0, _L: np.full_like(_L, 0.1),
        fh=lambda _z1, _z0h, _L: np.full_like(_L, 0.1),
    )
    model.input.surface = SimpleNamespace(
        z_m=10.0, z_o=0.1, z_s=2.0, z_t=0.01, zeta_max=5.0,
        gustiness=0.0, gustiness_stable_only=True,
    )
    model.input.grid = SimpleNamespace(z=np.array([0.0, 0.05]), nz=2, nx=1, ny=1)

    saved = (
        float(model.sfc_state.turbulence.obukhov_length[0]),
        float(model.sfc_state.turbulence.friction_velocity[0]),
        float(model.sfc_state.fluxes.ground_heat[0]),
        float(model.sfc_state.fluxes.sensible_heat[0]),
        float(model.sfc_state.fluxes.latent_heat[0]),
        float(model.sfc_state.fluxes.kinematic_heat[0]),
        float(model.sfc_state.fluxes.kinematic_moisture[0]),
    )

    # _compute_seb_vec is a pure function - same inputs give same outputs
    sfc_T = np.array([280.0])
    initial_L = np.array([2.0])
    # Testing internal solver behavior directly
    r1 = model._compute_seb_vec(sfc_T, initial_L)
    r2 = model._compute_seb_vec(sfc_T, initial_L)
    assert np.allclose(r1, r2)

    model.atm_state.seb_storage = np.array([10.0])
    r_stored = model._compute_seb_vec(sfc_T, initial_L)
    assert np.allclose(r_stored, r1 - 10.0)
    model.atm_state.seb_storage = np.array([0.0])

    # Pure function should not mutate state
    restored = (
        float(model.sfc_state.turbulence.obukhov_length[0]),
        float(model.sfc_state.turbulence.friction_velocity[0]),
        float(model.sfc_state.fluxes.ground_heat[0]),
        float(model.sfc_state.fluxes.sensible_heat[0]),
        float(model.sfc_state.fluxes.latent_heat[0]),
        float(model.sfc_state.fluxes.kinematic_heat[0]),
        float(model.sfc_state.fluxes.kinematic_moisture[0]),
    )
    assert restored == saved


def test_solve_seb_left_bracket_respects_temperature_floor() -> None:
    """SEB bracket expansion should not evaluate sub-physical temperatures."""
    model = _make_minimal_model()
    model.ncol = 1
    model.sfc_state = SurfaceState(
        temperature=np.array([101.0]),
        moisture=np.array([0.2]),
    )
    model.sfc_state.turbulence.obukhov_length[0] = 1.0
    model.soil_state = SoilState(
        temperature=np.array([[101.0], [101.0]]),
        moisture=np.array([[0.2], [0.2]]),
        type=np.array(["clay", "clay"], dtype=object),
    )
    model.solver_state = SimpleNamespace(conductivity_thermal_mid=np.array([1.0]))
    model.soil = SimpleNamespace(
        conductivity_thermal=lambda _moisture: np.ones((2, 1))
    )

    evaluated: list[np.ndarray] = []

    def fake_compute_seb_vec(
        sfc_T: np.ndarray,
        _initial_L: np.ndarray,
        cols: np.ndarray | None = None,
    ) -> np.ndarray:
        del cols
        temps = np.asarray(sfc_T, dtype=float)
        evaluated.append(temps.copy())
        return np.where(temps <= 100.0, 1.0, 2.0)

    model._compute_seb_vec = fake_compute_seb_vec

    with np.errstate(all="raise"):
        try:
            model._solve_seb()
        except SolverError:
            pass

    assert evaluated
    assert min(float(np.min(temps)) for temps in evaluated) == 100.0


def test_surface_coupling_recomputes_fluxes_on_convergence() -> None:
    """Recomputes fluxes after the final SMB update before returning."""
    model = _make_minimal_model(coupling_iterations=5)

    calls: list[tuple[float, float]] = []

    def fake_compute_fluxes(sfc_T: float, sfc_q: float) -> None:
        calls.append((sfc_T, sfc_q))

    # First iteration changes T and q, second iteration converges.
    state = {"iter": 0}

    def fake_solve_seb() -> None:
        model.sfc_state.temperature = 301.0

    def fake_solve_smb() -> None:
        if state["iter"] == 0:
            model.sfc_state.moisture = 0.1
        state["iter"] += 1

    model._compute_fluxes = fake_compute_fluxes
    model._solve_seb = fake_solve_seb
    model._solve_smb = fake_solve_smb

    model._solve_surface_coupling()

    assert calls, "Expected coupled solver to recompute fluxes on exit."
    assert calls[-1][0] == 301.0
    assert abs(calls[-1][1] - 0.1) < 1e-12


def test_surface_coupling_refreshes_canopy_diagnostics_on_convergence() -> None:
    """Final canopy diagnostics should reflect the converged state."""
    model = _make_minimal_model(coupling_iterations=5, tol_mois=1e-3)

    refreshes: list[tuple[float, float]] = []

    def fake_refresh() -> None:
        refreshes.append((
            float(np.asarray(model.sfc_state.temperature)),
            float(np.asarray(model.sfc_state.moisture)),
        ))

    state = {"iter": 0}

    def fake_solve_seb() -> None:
        model.sfc_state.temperature = 301.0

    def fake_solve_smb() -> None:
        model.sfc_state.moisture = 0.1 if state["iter"] == 0 else 0.1005
        state["iter"] += 1

    model._refresh_canopy_diagnostics = fake_refresh
    model._compute_fluxes = lambda *_args, **_kwargs: None
    model._solve_seb = fake_solve_seb
    model._solve_smb = fake_solve_smb

    model._solve_surface_coupling()

    assert len(refreshes) == 3
    assert refreshes[-1] == (301.0, 0.1005)


def test_surface_coupling_recomputes_fluxes_on_nonconvergence() -> None:
    """Recomputes fluxes even when the coupling loop hits its iteration cap."""
    model = _make_minimal_model(coupling_iterations=2, tol_temp=0.0, tol_mois=0.0)

    calls: list[tuple[float, float]] = []

    def fake_compute_fluxes(sfc_T: float, sfc_q: float) -> None:
        calls.append((sfc_T, sfc_q))

    # Never converges within 2 iterations because T keeps changing.
    def fake_solve_seb() -> None:
        model.sfc_state.temperature += 1.0

    def fake_solve_smb() -> None:
        model.sfc_state.moisture -= 0.01

    model._compute_fluxes = fake_compute_fluxes
    model._solve_seb = fake_solve_seb
    model._solve_smb = fake_solve_smb

    model._solve_surface_coupling()

    assert calls, "Expected coupled solver to recompute fluxes on exit."


def test_surface_coupling_refreshes_canopy_diagnostics_on_nonconvergence() -> None:
    """Final canopy diagnostics should refresh even on iteration-cap exit."""
    model = _make_minimal_model(coupling_iterations=2, tol_temp=0.0, tol_mois=0.0)

    refreshes: list[tuple[float, float]] = []

    def fake_refresh() -> None:
        refreshes.append((
            float(np.asarray(model.sfc_state.temperature)),
            float(np.asarray(model.sfc_state.moisture)),
        ))

    def fake_solve_seb() -> None:
        model.sfc_state.temperature += 1.0

    def fake_solve_smb() -> None:
        model.sfc_state.moisture -= 0.01

    model._refresh_canopy_diagnostics = fake_refresh
    model._compute_fluxes = lambda *_args, **_kwargs: None
    model._solve_seb = fake_solve_seb
    model._solve_smb = fake_solve_smb

    model._solve_surface_coupling()

    assert len(refreshes) == 3
    assert refreshes[-1] == (302.0, 0.18)
