"""Tests for the initial (time=0) output snapshot."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

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
            scalar_fields = ("ust", "obl", "shf", "lhf", "ghf")
            self.saved_initial = {
                name: float(np.asarray(fields[name])[0])
                for name in scalar_fields
                if name in fields
            }


def _make_model(
    *,
    warm_start: bool,
    has_forcing: bool,
    output_save: bool = True,
    output_fields: list[str] | None = None,
) -> UtahLSM:
    """Create a minimal UtahLSM instance for `_setup_output()` tests."""
    model = UtahLSM.__new__(UtahLSM)
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
    model._did_warm_start_turbulence = False

    model.input = _DummyNamelist(
        grid=SimpleNamespace(nz=2, nx=1, ny=1, z=np.array([0.1, 0.2])),
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
            warm_start_turbulence=warm_start,
            initialize_surface_temperature_from_seb=False,
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

    def fake_warm_start() -> None:
        model.sfc_state.turbulence.friction_velocity[0] = 0.123
        model.sfc_state.turbulence.obukhov_length[0] = 456.0
        model.sfc_state.fluxes.sensible_heat[0] = 7.0
        model.sfc_state.fluxes.latent_heat[0] = 8.0
        model.sfc_state.fluxes.ground_heat[0] = 9.0

    model._warm_start_turbulence = fake_warm_start  # type: ignore[attr-defined]
    return model


def test_initial_output_defaults_to_zeros():
    """Keeps the initial output snapshot at zeros by default."""
    model = _make_model(warm_start=False, has_forcing=True)
    model._setup_output()
    assert model.output.configured_fields == [
        "ust", "obl", "shf", "lhf", "ghf", "soil_z", "soil_T", "soil_q"
    ]
    assert model.output.saved_initial == {
        "ust": 0.0,
        "obl": 0.0,
        "shf": 0.0,
        "lhf": 0.0,
        "ghf": 0.0,
    }


def test_initial_output_warm_starts_when_enabled_and_forced():
    """Writes warm-started diagnostics into the initial output snapshot."""
    model = _make_model(warm_start=True, has_forcing=True)
    model._setup_output()
    assert model.output.saved_initial == {
        "ust": 0.123,
        "obl": 456.0,
        "shf": 7.0,
        "lhf": 8.0,
        "ghf": 9.0,
    }
    assert model._did_warm_start_turbulence is True
    assert (
        model.output.outfile.attrs.get("initial_diagnostics")
        == "warm_start_turbulence using forcing[0]"
    )


def test_initial_output_does_not_warm_start_without_forcing():
    """Does not warm-start initial output when forcing is unavailable."""
    model = _make_model(warm_start=True, has_forcing=False)
    model._setup_output()
    assert model.output.saved_initial == {
        "ust": 0.0,
        "obl": 0.0,
        "shf": 0.0,
        "lhf": 0.0,
        "ghf": 0.0,
    }


def test_initial_output_can_initialize_surface_temperature_from_seb():
    """When enabled, SEB initialization can populate initial diagnostics."""
    model = _make_model(warm_start=False, has_forcing=True)
    model.input.numerics = NumericsConfig(
        heat_diffusion_back_weight=(
            model.input.numerics.heat_diffusion_back_weight
        ),
        iterations=model.input.numerics.iterations,
        tolerances=model.input.numerics.tolerances,
        warm_start_turbulence=False,
        initialize_surface_temperature_from_seb=True,
    )

    def fake_solve_seb() -> None:
        model.sfc_state.temperature[:] = 280.0
        model.sfc_state.soil_top_temperature[:] = 280.0
        model.sfc_state.turbulence.friction_velocity[0] = 0.2
        model.sfc_state.turbulence.obukhov_length[0] = 50.0
        model.sfc_state.fluxes.sensible_heat[0] = -10.0
        model.sfc_state.fluxes.latent_heat[0] = 5.0
        model.sfc_state.fluxes.ground_heat[0] = -40.0

    model._solve_seb = fake_solve_seb  # type: ignore[assignment]

    model._setup_output()

    assert model.soil_state.temperature[0, 0] == 280.0
    assert model.output.saved_initial == {
        "ust": 0.2,
        "obl": 50.0,
        "shf": -10.0,
        "lhf": 5.0,
        "ghf": -40.0,
    }
    assert (
        model.output.outfile.attrs.get("initial_surface_temperature")
        == "initialized from SEB using forcing[0]"
    )


def test_setup_output_filters_requested_fields() -> None:
    """Writes only the explicitly requested output fields."""
    model = _make_model(
        warm_start=False,
        has_forcing=True,
        output_fields=["soil_z", "ust", "soil_T"],
    )

    model._setup_output()

    assert list(model.output_fields) == ["soil_z", "ust", "soil_T"]
    assert model.output.configured_fields == ["soil_z", "ust", "soil_T"]
    assert model.output.saved_initial == {"ust": 0.0}


def test_setup_output_skips_normal_fields_when_save_disabled() -> None:
    """Disables normal output when the namelist turns saving off."""
    model = _make_model(
        warm_start=False,
        has_forcing=True,
        output_save=False,
    )

    model._setup_output()

    assert model.output_fields == {}
    assert model.output.configured_fields == []
    assert model.output.saved_initial == {}


def test_setup_output_rejects_unknown_requested_field() -> None:
    """Raises when `output.fields` contains an unknown variable name."""
    model = _make_model(
        warm_start=False,
        has_forcing=True,
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
