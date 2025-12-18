"""Tests for the coupled surface SEB/SMB iteration.

These tests validate that the coupled surface solver leaves the diagnostic
fluxes consistent with the final converged surface state.
"""

import logging
from types import SimpleNamespace

from utahlsm.core import UtahLSM
from utahlsm.data_models import (
    IterationsConfig,
    NumericsConfig,
    SurfaceState,
    TolerancesConfig,
)


def _make_minimal_model(
    *,
    coupling_iterations: int = 3,
    tol_temp: float = 1e-6,
    tol_mois: float = 1e-12,
) -> UtahLSM:
    """Create a minimal UtahLSM instance for unit-testing private methods."""
    model = UtahLSM.__new__(UtahLSM)
    model.logger = logging.getLogger("test")
    model.input = SimpleNamespace(
        numerics=NumericsConfig(
            diffusion_back_weight=0.5,
            iterations=IterationsConfig(
                sfc_flux=1,
                seb_bracket=1,
                seb_root=1,
                smb_flux=1,
                coupling=coupling_iterations,
            ),
            tolerances=TolerancesConfig(
                sfc_flux=1e-6,
                seb_root=1e-12,
                smb_flux=1e-12,
                coupling_temp=tol_temp,
                coupling_mois=tol_mois,
            ),
        )
    )
    model.sfc_state = SurfaceState(temperature=300.0, moisture=0.2)
    return model


def test_compute_seb_is_deterministic_and_restores_state():
    """SEB residual evaluation should not leak mutated internal state."""
    model = _make_minimal_model()

    # Provide a minimal atmosphere and flux containers.
    model.atm_state = SimpleNamespace(radiation_net=100.0)
    model.sfc_state.fluxes.ground_heat[0] = 1.0
    model.sfc_state.fluxes.sensible_heat[0] = 2.0
    model.sfc_state.fluxes.latent_heat[0] = 3.0
    model.sfc_state.fluxes.kinematic_heat[0] = 4.0
    model.sfc_state.fluxes.kinematic_moisture[0] = 5.0
    model.sfc_state.turbulence.obukhov_length[0] = 7.0
    model.sfc_state.turbulence.friction_velocity[0] = 8.0

    # `_compute_fluxes` mutates state as the real one does; `_compute_seb`
    # should still be deterministic and restore these values after returning.
    def fake_compute_fluxes(_sfc_T: float, _sfc_q: float) -> None:
        model.sfc_state.fluxes.ground_heat[0] = 10.0
        model.sfc_state.fluxes.sensible_heat[0] = 20.0
        model.sfc_state.fluxes.latent_heat[0] = 30.0
        model.sfc_state.fluxes.kinematic_heat[0] = 40.0
        model.sfc_state.fluxes.kinematic_moisture[0] = 50.0
        model.sfc_state.turbulence.obukhov_length[0] = 70.0
        model.sfc_state.turbulence.friction_velocity[0] = 80.0

    model._compute_fluxes = fake_compute_fluxes  # type: ignore[attr-defined]

    saved = (
        float(model.sfc_state.turbulence.obukhov_length[0]),
        float(model.sfc_state.turbulence.friction_velocity[0]),
        float(model.sfc_state.fluxes.ground_heat[0]),
        float(model.sfc_state.fluxes.sensible_heat[0]),
        float(model.sfc_state.fluxes.latent_heat[0]),
        float(model.sfc_state.fluxes.kinematic_heat[0]),
        float(model.sfc_state.fluxes.kinematic_moisture[0]),
    )

    r1 = model._compute_seb(280.0, initial_obukhov_length=2.0)
    r2 = model._compute_seb(280.0, initial_obukhov_length=2.0)
    assert r1 == r2

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


def test_surface_coupling_recomputes_fluxes_on_convergence():
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

    model._compute_fluxes = fake_compute_fluxes  # type: ignore[attr-defined]
    model._solve_seb = fake_solve_seb  # type: ignore[attr-defined]
    model._solve_smb = fake_solve_smb  # type: ignore[attr-defined]

    model._solve_surface_coupling()

    assert calls, "Expected coupled solver to recompute fluxes on exit."
    assert calls[-1] == (301.0, 0.1)


def test_surface_coupling_recomputes_fluxes_on_nonconvergence():
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

    model._compute_fluxes = fake_compute_fluxes  # type: ignore[attr-defined]
    model._solve_seb = fake_solve_seb  # type: ignore[attr-defined]
    model._solve_smb = fake_solve_smb  # type: ignore[attr-defined]

    model._solve_surface_coupling()

    assert calls, "Expected coupled solver to recompute fluxes on exit."
