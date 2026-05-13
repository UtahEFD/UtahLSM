"""Tests for liquid precipitation forcing and saturation-excess runoff.

Exercises:
- The SMB residual with a P_infil term (sign convention, runoff branch).
- The NetCDF forcing loader: optional precip variable, validation,
  cold-temperature warning.
"""

from __future__ import annotations

import logging
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import netCDF4 as nc
import numpy as np
import pytest

from utahlsm.core import UtahLSM
from utahlsm.data_models import (
    AtmosphericState,
    IterationsConfig,
    NumericsConfig,
    SoilState,
    SurfaceFluxes,
    SurfaceState,
    TolerancesConfig,
    TurbulenceScales,
)
from utahlsm.physics.soil.soil_brookscorey import BrooksCorey
from utahlsm.util.io.input import Input
from utahlsm.util.io.soil_properties_loader import SoilPropertiesLoader


def _make_smb_model(
    *,
    theta_profile: float = 0.30,
    precipitation: float = 0.0,
    soil_type: str = "loam",
) -> Any:
    """Build a minimal UtahLSM with the dependencies needed by _solve_smb."""
    model: Any = UtahLSM.__new__(UtahLSM)
    model.logger = logging.getLogger("test.precip")
    model.ncol = 1

    model.input = SimpleNamespace(
        grid=SimpleNamespace(nx=1, ny=1, nz=3, z=np.array([0.05, 0.20, 0.35])),
        surface=SimpleNamespace(z_s=2.0, z_t=0.01),
        numerics=NumericsConfig(
            heat_diffusion_back_weight=1.0,
            iterations=IterationsConfig(
                sfc_flux=10,
                seb_bracket=10,
                seb_root=10,
                smb_flux=100,
                moisture_picard=50,
                coupling=10,
            ),
            tolerances=TolerancesConfig(
                sfc_flux=1e-3,
                seb_root=1e-6,
                smb_flux=1e-10,
                moisture_picard=1e-8,
                moisture_bounds=1e-6,
                coupling_temp=1e-3,
                coupling_mois=1e-8,
            ),
        ),
    )

    props = SoilPropertiesLoader.load("cosby")
    model.soil = BrooksCorey(props, [soil_type] * 3, "cosby")

    model.soil_state = SoilState(
        temperature=np.array([[293.15], [293.15], [293.15]]),
        moisture=np.array([[theta_profile], [theta_profile], [theta_profile]]),
        type=np.array([soil_type] * 3, dtype=object),
    )

    fluxes = SurfaceFluxes(
        kinematic_heat=np.zeros(1),
        kinematic_moisture=np.zeros(1),
        sensible_heat=np.zeros(1),
        latent_heat=np.zeros(1),
        ground_heat=np.zeros(1),
        runoff=np.zeros(1),
    )
    turb = TurbulenceScales(
        friction_velocity=np.array([0.3]),
        obukhov_length=np.array([1e6]),
    )
    model.sfc_state = SurfaceState(
        temperature=np.array([293.15]),
        moisture=np.array([theta_profile]),
        air_density=np.array([1.205]),
        fluxes=fluxes,
        turbulence=turb,
    )

    model.atm_state = AtmosphericState(
        wind_speed=np.array([3.0]),
        temperature=np.array([290.0]),
        specific_humidity=np.array([0.006]),
        pressure=np.array([101325.0]),
        precipitation=np.array([precipitation]),
    )

    # Neutral stability stub: fh constant. f_veg=0 (no canopy).
    model.sfc = SimpleNamespace(
        fh=lambda _z_s, _z_t, L: np.full_like(L, 0.1),
    )
    model.canopy = None

    return model


# ---------------------------------------------------------------------------
# SMB physics
# ---------------------------------------------------------------------------


@pytest.mark.precip
@pytest.mark.smb
def test_zero_precip_baseline_runoff_is_zero() -> None:
    """With P=0, SMB behaves as before: runoff stays zero, theta_sfc in bounds."""
    model = _make_smb_model(theta_profile=0.30, precipitation=0.0)
    porosity = float(model.soil.properties.porosity[0])
    residual_q = float(model.soil.properties.residual[0])

    model._solve_smb()

    assert np.all(model.sfc_state.fluxes.runoff == 0.0)
    assert np.all(model.sfc_state.moisture >= residual_q)
    assert np.all(model.sfc_state.moisture <= porosity)


@pytest.mark.precip
@pytest.mark.smb
def test_light_rain_no_runoff() -> None:
    """A light rain rate produces no runoff (soil absorbs it all)."""
    model = _make_smb_model(theta_profile=0.30, precipitation=1e-5)
    porosity = float(model.soil.properties.porosity[0])

    model._solve_smb()

    assert np.all(model.sfc_state.fluxes.runoff == 0.0)
    assert float(model.sfc_state.moisture[0]) <= porosity


@pytest.mark.precip
@pytest.mark.smb
def test_heavy_rain_produces_runoff() -> None:
    """When P exceeds the saturated-drainage capacity, the excess is runoff."""
    probe = _make_smb_model()
    K_sat_top = float(probe.soil.properties.K_sat[0])
    RHO_W = 1000.0
    infil_cap = RHO_W * K_sat_top  # kg/m^2/s

    P_heavy = infil_cap * 10.0  # 10x the soil's drainage capacity
    model = _make_smb_model(theta_profile=0.30, precipitation=P_heavy)
    model._solve_smb()

    runoff = float(model.sfc_state.fluxes.runoff[0])
    assert runoff == pytest.approx(P_heavy - infil_cap, rel=1e-12)
    # Sanity: runoff cannot exceed total precipitation.
    assert 0.0 < runoff <= P_heavy


@pytest.mark.precip
@pytest.mark.smb
def test_rain_below_capacity_no_runoff() -> None:
    """Rain at or below the saturated-drainage capacity produces zero runoff."""
    probe = _make_smb_model()
    K_sat_top = float(probe.soil.properties.K_sat[0])
    RHO_W = 1000.0
    P_modest = 0.5 * RHO_W * K_sat_top  # half of capacity

    model = _make_smb_model(theta_profile=0.30, precipitation=P_modest)
    model._solve_smb()

    assert float(model.sfc_state.fluxes.runoff[0]) == 0.0


@pytest.mark.precip
@pytest.mark.smb
def test_precip_field_propagates_via_load_atm_state() -> None:
    """update() -> _load_atm_state copies precipitation onto self.atm_state."""
    model = _make_smb_model(precipitation=0.0)
    new_atm = AtmosphericState(
        wind_speed=np.array([3.0]),
        temperature=np.array([290.0]),
        specific_humidity=np.array([0.006]),
        pressure=np.array([101325.0]),
        sw_in=np.array([0.0]),
        lw_in=np.array([300.0]),
        seb_storage=np.array([0.0]),
        precipitation=np.array([2.5e-4]),
    )

    model._load_atm_state(new_atm)

    assert np.allclose(model.atm_state.precipitation, 2.5e-4)


# ---------------------------------------------------------------------------
# NetCDF forcing loader
# ---------------------------------------------------------------------------


def _make_input_loader() -> Input:
    """Build a bare Input instance with just the attrs the loader needs."""
    obj = Input.__new__(Input)
    obj.logger = logging.getLogger("tests.precip.input")
    obj.grid = SimpleNamespace(nx=1, ny=1)
    return obj


def _write_forcing_nc(
    path: Path,
    *,
    atm_T: np.ndarray,
    precip: np.ndarray | None,
    ntime: int = 3,
) -> None:
    """Write a minimal offline forcing NetCDF for loader tests."""
    with nc.Dataset(path, "w") as ds:
        ds.createDimension("t", ntime)
        ds.createDimension("scalar", 1)
        ds.createVariable("tstep", "f8", ("scalar",))[:] = [60.0]
        for name, vals in (
            ("atm_U", np.full(ntime, 3.0)),
            ("atm_T", atm_T),
            ("atm_q", np.full(ntime, 0.006)),
            ("atm_p", np.full(ntime, 101325.0)),
            ("sw_in", np.full(ntime, 300.0)),
            ("lw_in", np.full(ntime, 350.0)),
        ):
            ds.createVariable(name, "f8", ("t",))[:] = vals
        if precip is not None:
            ds.createVariable("precip", "f8", ("t",))[:] = precip


@pytest.mark.precip
@pytest.mark.unit
def test_precip_optional_in_netcdf(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """A forcing file without `precip` loads cleanly with zero defaults."""
    path = tmp_path / "lsm_offline.nc"
    _write_forcing_nc(path, atm_T=np.full(3, 290.0), precip=None)

    loader = _make_input_loader()
    with caplog.at_level(logging.INFO, logger=loader.logger.name):
        loader._load_offline_data(str(path))

    assert loader.forcing is not None
    assert loader.forcing.ntime == 3
    for atm in loader.forcing.atmos:
        assert np.allclose(np.asarray(atm.precipitation), 0.0)
    assert any("precip not in forcing" in rec.getMessage() for rec in caplog.records)


@pytest.mark.precip
@pytest.mark.unit
def test_precip_validation_rejects_negative(tmp_path: Path) -> None:
    """Negative precipitation in the NetCDF must raise."""
    path = tmp_path / "lsm_offline.nc"
    precip = np.array([0.0, -1.0, 0.0])
    _write_forcing_nc(path, atm_T=np.full(3, 290.0), precip=precip)

    loader = _make_input_loader()
    with pytest.raises(ValueError, match="precip"):
        loader._load_offline_data(str(path))


@pytest.mark.precip
@pytest.mark.unit
def test_cold_temp_precip_warning(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """T < 273.15 K with P > 0 fires a single load-time WARNING."""
    path = tmp_path / "lsm_offline.nc"
    atm_T = np.array([270.0, 268.0, 285.0])
    precip = np.array([1e-5, 2e-5, 1e-5])
    _write_forcing_nc(path, atm_T=atm_T, precip=precip)

    loader = _make_input_loader()
    with caplog.at_level(logging.WARNING, logger=loader.logger.name):
        loader._load_offline_data(str(path))

    warnings = [r for r in caplog.records
                if r.levelno == logging.WARNING
                and "Cold-temperature precipitation" in r.getMessage()]
    assert len(warnings) == 1
    assert "2 sample" in warnings[0].getMessage()
