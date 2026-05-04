"""Tests for the optional gustiness wind-speed adjustment."""

import logging
from types import SimpleNamespace
from typing import Any

import numpy as np

from utahlsm.core import UtahLSM
from utahlsm.data_models import AtmosphericState, SoilState, SurfaceState


def _make_minimal_model(*, gustiness: float, stable_only: bool, L0: float) -> Any:
    model: Any = UtahLSM.__new__(UtahLSM)
    model.logger = logging.getLogger("test")
    model.ncol = 1
    model.atm_state = AtmosphericState(
        wind_speed=np.array([3.0]),
        temperature=np.array([290.0]),
        specific_humidity=np.array([0.0]),
        pressure=np.array([101000.0]),
        radiation_net=np.array([-50.0]),
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
    model.sfc_state.turbulence.obukhov_length[0] = L0

    model.soil = SimpleNamespace(surface_specific_humidity=lambda _T, _q, _p: np.zeros_like(_T))
    model.solver_state = SimpleNamespace(conductivity_thermal_mid=np.array([1.0]))
    model.sfc = SimpleNamespace(
        fm=lambda _z1, _z0, _L: np.full_like(_L, 0.1),
        fh=lambda _z1, _z0h, _L: np.full_like(_L, 0.1),
    )
    model.input = SimpleNamespace(
        surface=SimpleNamespace(
            z_m=10.0,
            z_o=0.1,
            z_s=2.0,
            z_t=0.01,
            zeta_max=5.0,
            gustiness=gustiness,
            gustiness_stable_only=stable_only,
        ),
        grid=SimpleNamespace(z=np.array([0.0, 0.05]), nz=2, nx=1, ny=1),
        numerics=SimpleNamespace(
            tolerances=SimpleNamespace(sfc_flux=0.0),
            iterations=SimpleNamespace(sfc_flux=1),
        ),
    )
    return model


def test_gustiness_increases_ustar_when_stable() -> None:
    """Applies gustiness when L>=0 and stable_only is True."""
    model = _make_minimal_model(gustiness=1.0, stable_only=True, L0=10.0)
    model._compute_fluxes(model.sfc_state.temperature, model.sfc_state.moisture)
    expected = np.hypot(3.0, 1.0) * 0.1
    assert np.isclose(model.sfc_state.turbulence.friction_velocity[0], expected)


def test_gustiness_not_applied_when_unstable_if_stable_only() -> None:
    """Does not apply gustiness when L<0 and stable_only is True."""
    model = _make_minimal_model(gustiness=1.0, stable_only=True, L0=-10.0)
    model._compute_fluxes(model.sfc_state.temperature, model.sfc_state.moisture)
    expected = 3.0 * 0.1
    assert np.isclose(model.sfc_state.turbulence.friction_velocity[0], expected)


def test_gustiness_always_applied_when_not_stable_only() -> None:
    """Applies gustiness regardless of stability when stable_only is False."""
    model = _make_minimal_model(gustiness=1.0, stable_only=False, L0=-10.0)
    model._compute_fluxes(model.sfc_state.temperature, model.sfc_state.moisture)
    expected = np.hypot(3.0, 1.0) * 0.1
    assert np.isclose(model.sfc_state.turbulence.friction_velocity[0], expected)

