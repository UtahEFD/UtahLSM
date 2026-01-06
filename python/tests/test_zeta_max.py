"""Tests for configurable MOST zeta (z/L) limits."""

import logging
from types import SimpleNamespace

import numpy as np

from utahlsm.core import UtahLSM
from utahlsm.data_models import AtmosphericState, SoilState, SurfaceState


def test_zeta_max_changes_obukhov_length_clamp():
    """Changing zeta_max changes the applied L clamp under strong stability."""
    model = UtahLSM.__new__(UtahLSM)
    model.logger = logging.getLogger("test")
    model.ncol = 1
    model.tstep = 600.0
    model.atm_state = AtmosphericState(
        wind_speed=np.array([3.0]),
        temperature=np.array([290.0]),
        specific_humidity=np.array([0.005]),
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

    # Provide minimal soil/surface functions used by _compute_fluxes.
    model.soil = SimpleNamespace(
        surface_mixing_ratio=lambda _T, _q, _p: np.zeros_like(_T),
    )
    model.solver_state = SimpleNamespace(conductivity_thermal_mid=np.array([1.0]))
    model.sfc = SimpleNamespace(
        fm=lambda _z1, _z0, _L: np.full_like(_L, 0.05),
        fh=lambda _z1, _z0h, _L: np.full_like(_L, 0.1),
    )

    # High stability: make flux_wTv negative and non-zero.
    # With fm=0.1 and U=3, u*=0.3. With fh=0.1 and (Ts-Ta)=-10K, flux_wT=-0.3.
    # That yields a small positive L, which should then be clamped by zeta_max.
    model.input = SimpleNamespace(
        surface=SimpleNamespace(
            z_m=10.0,
            z_o=0.1,
            z_s=2.0,
            z_t=0.01,
            zeta_max=2.0,
            gustiness=0.0,
            gustiness_stable_only=True,
        ),
        grid=SimpleNamespace(z=np.array([0.0, 0.05]), nz=2),
        numerics=SimpleNamespace(
            tolerances=SimpleNamespace(sfc_flux=0.0),
            iterations=SimpleNamespace(sfc_flux=1),
        ),
    )

    model._compute_fluxes(model.sfc_state.temperature, model.sfc_state.moisture)
    assert np.isclose(model.sfc_state.turbulence.obukhov_length[0], 5.0)  # 10/2
