"""Helpers for constructing temporary multi-column GABLS3 test cases."""

from __future__ import annotations

import json
from pathlib import Path

import netCDF4 as nc
import numpy as np

import utahlsm


def write_2x2_case(single_case: Path, multi_case: Path) -> None:
    """Write a 2x2 clone of the single-column GABLS3 case."""
    multi_case.mkdir()

    with open(single_case / "lsm_namelist.json", encoding="utf-8") as infile:
        namelist = json.load(infile)
    namelist["grid"]["nx"] = 2
    namelist["grid"]["ny"] = 2
    with open(multi_case / "lsm_namelist.json", "w", encoding="utf-8") as outfile:
        json.dump(namelist, outfile, indent=4)
        outfile.write("\n")

    with (
        nc.Dataset(single_case / "lsm_init.nc") as src,
        nc.Dataset(multi_case / "lsm_init.nc", "w") as dst,
    ):
        nz = len(src.dimensions["z"])
        dst.createDimension("z", nz)
        dst.createDimension("y", 2)
        dst.createDimension("x", 2)

        soil_z = np.asarray(src.variables["soil_z"][:], dtype=float)
        dst.createVariable("soil_z", "f8", ("z",))[:] = soil_z

        for name in ("soil_T", "soil_q"):
            data = np.asarray(src.variables[name][:], dtype=float)
            tiled = np.broadcast_to(data[:, None, None], (nz, 2, 2))
            dst.createVariable(name, "f8", ("z", "y", "x"))[:] = tiled

        soil_type = np.asarray(src.variables["soil_type"][:], dtype=object)
        tiled_type = np.broadcast_to(soil_type[:, None, None], (nz, 2, 2))
        dst.createVariable("soil_type", str, ("z", "y", "x"))[:] = tiled_type

    with (
        nc.Dataset(single_case / "lsm_offline.nc") as src,
        nc.Dataset(multi_case / "lsm_offline.nc", "w") as dst,
    ):
        ntime = len(src.dimensions["t"])
        dst.createDimension("t", ntime)
        dst.createDimension("y", 2)
        dst.createDimension("x", 2)
        dst.createDimension("scalar", 1)

        tstep = float(np.asarray(src.variables["tstep"][:], dtype=float))
        dst.createVariable("tstep", "f8", ("scalar",))[:] = [tstep]

        for name in ("atm_U", "atm_T", "atm_q", "atm_p",
                     "sw_in", "sw_out", "lw_in", "lw_out"):
            data = np.asarray(src.variables[name][:], dtype=float)
            tiled = np.broadcast_to(data[:, None, None], (ntime, 2, 2))
            dst.createVariable(name, "f8", ("t", "y", "x"))[:] = tiled


def run_case(case_dir: Path, outfile: Path) -> None:
    """Run an offline case through the public UtahLSM Python API."""
    input_lsm = utahlsm.Input(
        str(case_dir / "lsm_namelist.json"),
        str(case_dir / "lsm_init.nc"),
        str(case_dir / "lsm_offline.nc"),
    )
    output_lsm = utahlsm.Output(str(outfile))
    try:
        lsm = utahlsm.UtahLSM(input_lsm, output_lsm)
        runtime = 0.0
        tstep = input_lsm.forcing.tstep
        for step_count, atm_state in enumerate(input_lsm.forcing.atmos):
            runtime += tstep
            lsm.update(tstep, runtime, atm_state)
            lsm.run()
            lsm.save(step_count + 1, runtime)
    finally:
        output_lsm.close()
