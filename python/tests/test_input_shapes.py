"""Input shape normalization tests for single- and multi-column cases."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import utahlsm
from tests.gabls3_case_factory import write_2x2_case


@pytest.mark.integration
def test_input_shapes_single_column() -> None:
    """Single-column inputs should load with explicit column dimensions."""
    repo_root = Path(__file__).resolve().parents[1]
    cases_root = (repo_root / ".." / "cases").resolve()
    case = cases_root / "gabls3"

    if not case.exists():
        pytest.skip("GABLS3 case files not available.")

    input_lsm = utahlsm.Input(
        str(case / "lsm_namelist.json"),
        str(case / "lsm_init.nc"),
        str(case / "lsm_offline.nc"),
    )

    nz = input_lsm.grid.nz
    assert np.asarray(input_lsm.initial.temperature).shape == (nz, 1)
    assert np.asarray(input_lsm.initial.moisture).shape == (nz, 1)

    assert input_lsm.forcing is not None
    first = input_lsm.forcing.atmos[0]
    assert np.asarray(first.wind_speed).shape == (1,)
    assert np.asarray(first.temperature).shape == (1,)
    assert np.asarray(first.specific_humidity).shape == (1,)
    assert np.asarray(first.pressure).shape == (1,)
    assert np.asarray(first.sw_in).shape == (1,)
    assert np.asarray(first.lw_in).shape == (1,)
    assert np.asarray(first.seb_storage).shape == (1,)


@pytest.mark.integration
def test_input_shapes_multi_column(tmp_path: Path) -> None:
    """Generated multi-column inputs should load with flattened columns."""
    repo_root = Path(__file__).resolve().parents[1]
    cases_root = (repo_root / ".." / "cases").resolve()
    case_single = cases_root / "gabls3"

    if not case_single.exists():
        pytest.skip("GABLS3 case files not available.")

    case = tmp_path / "gabls3_2x2_generated"
    write_2x2_case(case_single, case)

    input_lsm = utahlsm.Input(
        str(case / "lsm_namelist.json"),
        str(case / "lsm_init.nc"),
        str(case / "lsm_offline.nc"),
    )

    nz = input_lsm.grid.nz
    ncol = input_lsm.grid.nx * input_lsm.grid.ny
    assert np.asarray(input_lsm.initial.temperature).shape == (nz, ncol)
    assert np.asarray(input_lsm.initial.moisture).shape == (nz, ncol)

    assert input_lsm.forcing is not None
    first = input_lsm.forcing.atmos[0]
    assert np.asarray(first.wind_speed).shape == (ncol,)
    assert np.asarray(first.temperature).shape == (ncol,)
    assert np.asarray(first.specific_humidity).shape == (ncol,)
    assert np.asarray(first.pressure).shape == (ncol,)
    assert np.asarray(first.sw_in).shape == (ncol,)
    assert np.asarray(first.lw_in).shape == (ncol,)
    assert np.asarray(first.seb_storage).shape == (ncol,)
