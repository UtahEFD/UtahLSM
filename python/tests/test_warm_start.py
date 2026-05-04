"""Tests for optional warm-start behavior on the first predictive step."""

import logging
from types import SimpleNamespace

import numpy as np

from utahlsm.core import UtahLSM
from utahlsm.data_models import (
    IterationsConfig,
    NumericsConfig,
    SoilState,
    SurfaceState,
    TolerancesConfig,
)


def _make_model_for_run(*, warm_start_turbulence: bool) -> UtahLSM:
    """Create a minimal UtahLSM instance for unit-testing `run()` behavior."""
    model = UtahLSM.__new__(UtahLSM)  # type: ignore[attr-defined]
    model.logger = logging.getLogger("test")

    model.input = SimpleNamespace(  # type: ignore[assignment]
        numerics=NumericsConfig(
            heat_diffusion_back_weight=0.5,
            warm_start_turbulence=warm_start_turbulence,
            initialize_surface_temperature_from_seb=False,
            iterations=IterationsConfig(
                sfc_flux=10,
                seb_bracket=10,
                seb_root=10,
                smb_flux=10,
                moisture_picard=10,
                coupling=2,
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
        )
    )

    model.soil_state = SoilState(
        temperature=np.array([290.0, 289.0]),
        moisture=np.array([0.25, 0.25]),
        type=np.array(["clay", "clay"], dtype=object),
    )
    model.sfc_state = SurfaceState()
    model._did_warm_start_turbulence = False  # type: ignore[attr-defined]

    # No-op the rest of the time step to isolate warm-start call conditions.
    model._solve_surface_coupling = lambda: None  # type: ignore[assignment]
    model._solve_diffusion_heat = lambda: None  # type: ignore[assignment]
    model._solve_diffusion_mois = lambda: None  # type: ignore[assignment]

    return model


def test_run_does_not_warm_start_by_default():
    """Does not warm-start unless explicitly enabled in numerics config."""
    model = _make_model_for_run(warm_start_turbulence=False)

    called = {"count": 0}

    def fake_warm_start() -> None:
        called["count"] += 1

    model._warm_start_turbulence = fake_warm_start  # type: ignore[attr-defined]

    model.run()
    assert called["count"] == 0


def test_run_warm_starts_only_once_when_enabled():
    """Warm-start executes only on the first call to `run()`."""
    model = _make_model_for_run(warm_start_turbulence=True)

    called = {"count": 0}

    def fake_warm_start() -> None:
        called["count"] += 1

    model._warm_start_turbulence = fake_warm_start  # type: ignore[attr-defined]

    model.run()
    model.run()
    assert called["count"] == 1


def test_run_calls_warm_start_before_coupling():
    """Warm-start happens before the coupled surface solve."""
    model = _make_model_for_run(warm_start_turbulence=True)

    sequence: list[str] = []

    def fake_warm_start() -> None:
        sequence.append("warm_start")

    def fake_coupling() -> None:
        sequence.append("coupling")

    model._warm_start_turbulence = fake_warm_start  # type: ignore[attr-defined]
    model._solve_surface_coupling = fake_coupling  # type: ignore[assignment]

    model.run()
    assert sequence[:2] == ["warm_start", "coupling"]
