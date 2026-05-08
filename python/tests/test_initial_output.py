"""Tests for the initial (time=0) output snapshot."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np

from utahlsm.core import UtahLSM
from utahlsm.data_models import (
    AtmosphericState,
    IterationsConfig,
    NumericsConfig,
    OutputConfig,
    SoilState,
    SurfaceState,
    TolerancesConfig,
)
from utahlsm.util.io.output import Output


@dataclass
class _DummyNamelist:
    """Minimal input container for `_setup_output()` testing."""

    grid: SimpleNamespace
    numerics: NumericsConfig
    output: OutputConfig
    soil_type_names: list[str]
    forcing: SimpleNamespace | None = None


class _DummyNetcdf:
    """Dummy NetCDF handle used by the Output stub."""

    def __init__(self) -> None:
        self.attrs: dict[str, str] = {}

    def setncattr(self, name: str, value: str) -> None:
        self.attrs[name] = value


class _DummyOutput:
    """Output stub that records the initial save snapshot."""

    def __init__(self) -> None:
        self.saved_initial: dict[str, np.ndarray | float] | None = None
        self.configured_fields: list[str] = []
        self.outfile = _DummyNetcdf()

    def set_dims(self, _dims: dict[str, int]) -> None:
        return

    def set_fields(self, fields: dict[str, object]) -> None:
        self.configured_fields = list(fields)

    def save(
        self,
        fields: dict[str, object],
        _tidx: int,
        _time: float,
        initial: bool = False,
    ) -> None:
        if initial:
            scalar_fields = ("ust", "obl", "shf", "lhf", "ghf", "seb_res")
            self.saved_initial = {
                name: float(np.asarray(fields[name])[0])
                for name in scalar_fields
                if name in fields
            }


def _make_model(
    *,
    has_forcing: bool,
    output_save: bool = True,
    output_fields: list[str] | None = None,
) -> Any:
    """Create a minimal UtahLSM instance for `_setup_output()` tests."""
    model: Any = UtahLSM.__new__(UtahLSM)
    model.ncol = 1
    model.output = _DummyOutput()
    model.sfc_state = SurfaceState(
        temperature=np.array([290.0]),
        moisture=np.array([0.25]),
    )
    model.soil_state = SoilState(
        temperature=np.array([[290.0], [289.0]]),
        moisture=np.array([[0.25], [0.25]]),
        type=np.array(["clay", "clay"], dtype=object),
    )
    model.atm_state = AtmosphericState(
        wind_speed=np.array([0.0]),
        temperature=np.array([0.0]),
        specific_humidity=np.array([0.0]),
        pressure=np.array([0.0]),
        radiation_net=np.array([0.0]),
    )
    model.solver_state = SimpleNamespace(conductivity_thermal_mid=np.array([1.0]))

    model.input = _DummyNamelist(
        grid=SimpleNamespace(nz=2, nx=1, ny=1, z=np.array([0.1, 0.2])),
        soil_type_names=["clay", "clay"],
        numerics=NumericsConfig(
            heat_diffusion_back_weight=0.5,
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
        ),
        output=OutputConfig(
            save=output_save,
            fields=output_fields or ["all"],
        ),
        forcing=(
            SimpleNamespace(
                atmos=[AtmosphericState(wind_speed=5.0)], tstep=600.0)
            if has_forcing
            else None
        ),
    )
    return model


def test_initial_output_runs_full_coupling_when_forced() -> None:
    """Drives SEB+SMB coupling at forcing[0] and writes converged diagnostics."""
    model = _make_model(has_forcing=True)

    def fake_coupling() -> None:
        model.sfc_state.temperature[:] = 280.0
        model.sfc_state.soil_top_temperature[:] = 280.0
        model.sfc_state.turbulence.friction_velocity[0] = 0.2
        model.sfc_state.turbulence.obukhov_length[0] = 50.0
        model.sfc_state.fluxes.sensible_heat[0] = -10.0
        model.sfc_state.fluxes.latent_heat[0] = 5.0
        model.sfc_state.fluxes.ground_heat[0] = -40.0

    model._solve_surface_coupling = fake_coupling

    model._setup_output()

    assert model.soil_state.temperature[0, 0] == 280.0
    assert model.output.saved_initial == {
        "ust": 0.2,
        "obl": 50.0,
        "shf": -10.0,
        "lhf": 5.0,
        "ghf": -40.0,
        "seb_res": 0.0,
    }
    assert (
        model.output.outfile.attrs.get("initial_state")
        == "SEB+SMB coupling using forcing[0]"
    )


def test_initial_output_skips_coupling_without_forcing() -> None:
    """Leaves the initial output at zeros when no forcing is available."""
    model = _make_model(has_forcing=False)

    called = {"count": 0}

    def fake_coupling() -> None:
        called["count"] += 1

    model._solve_surface_coupling = fake_coupling
    model._setup_output()

    assert called["count"] == 0
    assert model.output.saved_initial == {
        "ust": 0.0,
        "obl": 0.0,
        "shf": 0.0,
        "lhf": 0.0,
        "ghf": 0.0,
        "seb_res": 0.0,
    }
    assert "initial_state" not in model.output.outfile.attrs


def test_setup_output_filters_requested_fields() -> None:
    """Writes only the explicitly requested output fields."""
    model = _make_model(
        has_forcing=False,
        output_fields=["soil_z", "ust", "soil_T"],
    )

    model._setup_output()

    assert list(model.output_fields) == ["soil_z", "ust", "soil_T"]
    assert model.output.configured_fields == ["soil_z", "ust", "soil_T"]
    assert model.output.saved_initial == {"ust": 0.0}


def test_setup_output_skips_normal_fields_when_save_disabled() -> None:
    """Disables normal output when the namelist turns saving off."""
    model = _make_model(
        has_forcing=False,
        output_save=False,
    )

    model._setup_output()

    assert model.output_fields == {}
    assert model.output.configured_fields == []
    assert model.output.saved_initial == {}


def test_setup_output_rejects_unknown_requested_field() -> None:
    """Raises when `output.fields` contains an unknown variable name."""
    model = _make_model(
        has_forcing=False,
        output_fields=["bogus_field"],
    )

    with np.testing.assert_raises_regex(ValueError, "Unknown output field"):
        model._setup_output()


def test_disabled_output_does_not_create_netcdf_file(tmp_path: Path) -> None:
    """Leaves no NetCDF file behind when output is disabled."""
    outfile = tmp_path / "disabled.nc"
    output = Output(str(outfile), enabled=False)

    output.set_dims({"t": 0, "z": 1})
    output.set_fields({"soil_z": np.array([0.1])})
    output.save({"soil_z": np.array([0.1])}, 0, 0.0, initial=True)
    output.close()

    assert not outfile.exists()
