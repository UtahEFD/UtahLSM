"""Input validation tests for namelist schema and forcing bounds."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import jsonschema
import netCDF4 as nc
import numpy as np
import pytest
from numpy.testing import assert_allclose

from utahlsm.util.io.input import Input


def _base_namelist() -> dict[str, Any]:
    """Returns a minimal namelist that satisfies the input schema."""
    return {
        "general": {"log_level": "info"},
        "numerics": {
            "heat_diffusion_back_weight": 0.5,
            "iterations": {
                "sfc_flux": 100,
                "seb_bracket": 100,
                "seb_root": 100,
                "smb_flux": 100,
                "moisture_picard": 25,
                "coupling": 50,
            },
            "tolerances": {
                "sfc_flux": 1e-3,
                "seb_root": 1e-6,
                "smb_flux": 1e-6,
                "moisture_picard": 1e-8,
                "moisture_bounds": 1e-6,
                "coupling_temp": 1e-2,
                "coupling_mois": 1e-5,
            },
        },
        "time": {
            "utc_start": 0,
            "utc_year": 2006,
            "julian_day": 183,
        },
        "grid": {
            "nx": 1,
            "ny": 1,
            "nz": 3,
        },
        "surface": {
            "z_o": 0.03,
            "z_t": 0.0046,
            "z_m": 10.0,
            "z_s": 2.0,
            "albedo": 0.33,
            "emissivity": 0.99,
            "model": 1,
        },
        "soil": {
            "properties": "rawls-brakensiek",
            "model": 2,
        },
        "radiation": {
            "model": 0,
            "latitude": 40.0,
            "longitude": -111.0,
        },
        "output": {
            "save": False,
            "fields": ["all"],
        },
    }


def _make_input() -> Input:
    """Builds a bare Input instance for direct validator tests."""
    input_obj = Input.__new__(Input)
    input_obj.logger = logging.getLogger("tests.input_validation")
    return input_obj


def test_load_and_validate_namelist_rejects_unknown_keys(
        tmp_path: Path) -> None:
    """Rejects namelist options that are not defined by the schema."""
    namelist = _base_namelist()
    namelist["numerics"]["warm_start_coupling"] = 10
    namelist_path = tmp_path / "lsm_namelist.json"
    namelist_path.write_text(json.dumps(namelist), encoding="utf-8")

    input_obj = _make_input()

    with pytest.raises(jsonschema.ValidationError, match="warm_start_coupling"):
        input_obj._load_and_validate_namelist(str(namelist_path))


def test_validate_forcing_data_clips_small_boundary_excursions() -> None:
    """Clips only small out-of-range forcing values back to valid bounds."""
    input_obj = _make_input()
    atm_U = np.array([[0.0, 50.2]])
    atm_T = np.array([[199.8, 350.2]])
    atm_q = np.array([[-5e-5, 0.0502]])
    atm_p = np.array([[49950.0, 110050.0]])
    sw_in = np.array([[-2.0, 1410.0]])
    sw_out = np.array([[-2.0, 1410.0]])
    lw_in = np.array([[95.0, 605.0]])
    lw_out = np.array([[95.0, 705.0]])
    r_net = np.array([[-210.0, 1210.0]])

    input_obj._validate_forcing_data(
        atm_U, atm_T, atm_q, atm_p,
        sw_in, sw_out, lw_in, lw_out,
        r_net, _ntime=1)

    assert_allclose(atm_U, np.array([[1e-4, 50.0]]))
    assert_allclose(atm_T, np.array([[200.0, 350.0]]))
    assert_allclose(atm_q, np.array([[0.0, 0.05]]))
    assert_allclose(atm_p, np.array([[50000.0, 110000.0]]))
    assert_allclose(sw_in, np.array([[0.0, 1400.0]]))
    assert_allclose(sw_out, np.array([[0.0, 1400.0]]))
    assert_allclose(lw_in, np.array([[100.0, 600.0]]))
    assert_allclose(lw_out, np.array([[100.0, 700.0]]))
    assert_allclose(r_net, np.array([[-200.0, 1200.0]]))


def test_validate_forcing_data_raises_on_large_violation() -> None:
    """Raises when forcing data are well outside the supported bounds."""
    input_obj = _make_input()

    with pytest.raises(ValueError, match="temperature"):
        input_obj._validate_forcing_data(
            atm_U=np.array([[2.0]]),
            atm_T=np.array([[190.0]]),
            atm_q=np.array([[0.01]]),
            atm_p=np.array([[101325.0]]),
            sw_in=np.array([[300.0]]),
            sw_out=np.array([[60.0]]),
            lw_in=np.array([[350.0]]),
            lw_out=np.array([[400.0]]),
            r_net=np.array([[190.0]]),
            _ntime=1,
        )


def test_load_and_validate_namelist_requires_general_and_numerics(
        tmp_path: Path) -> None:
    """Rejects namelists that omit sections required by the parser."""
    namelist = _base_namelist()
    namelist.pop("general")
    namelist.pop("numerics")
    namelist_path = tmp_path / "lsm_namelist.json"
    namelist_path.write_text(json.dumps(namelist), encoding="utf-8")

    input_obj = _make_input()

    with pytest.raises(jsonschema.ValidationError):
        input_obj._load_and_validate_namelist(str(namelist_path))


def test_load_and_validate_namelist_rejects_invalid_canopy_bounds(
        tmp_path: Path) -> None:
    """Rejects canopy parameters that violate physical bounds."""
    namelist = _base_namelist()
    namelist["canopy"] = {
        "model": "jarvis",
        "veg_fraction": 1.1,
    }
    namelist_path = tmp_path / "lsm_namelist.json"
    namelist_path.write_text(json.dumps(namelist), encoding="utf-8")

    input_obj = _make_input()

    with pytest.raises(jsonschema.ValidationError, match="maximum of 1"):
        input_obj._load_and_validate_namelist(str(namelist_path))


def test_load_initial_conditions_rejects_nonuniform_soil_z(
        tmp_path: Path) -> None:
    """Rejects init files whose soil_z spacing is not uniform."""
    init_path = tmp_path / "lsm_init.nc"
    with nc.Dataset(init_path, "w") as ds:
        ds.createDimension("z", 3)
        ds.createVariable("soil_z", "f8", ("z",))[:] = [0.0, 0.05, 0.15]
        ds.createVariable("soil_T", "f8", ("z",))[:] = [290.0, 289.0, 288.0]
        ds.createVariable("soil_q", "f8", ("z",))[:] = [0.25, 0.25, 0.25]
        ds.createVariable("soil_type", str, ("z",))[:] = np.asarray(
            ["clay", "clay", "clay"], dtype=object
        )

    input_obj = _make_input()

    with pytest.raises(ValueError, match="uniform spacing"):
        input_obj._load_initial_conditions(str(init_path), nx=1, ny=1, nz=3)
