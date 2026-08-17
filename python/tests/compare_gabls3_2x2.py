"""Helpers for comparing GABLS3 single vs 2x2 outputs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import netCDF4 as nc
import numpy as np


@dataclass(frozen=True)
class ComparisonResult:
    """Summary of a single variable comparison."""

    name: str
    matches: bool
    max_abs_diff: float
    max_abs_col_diff: float | None


def compare_outputs(
    *,
    single_path: Path,
    multi_path: Path,
    rtol: float = 0.0,
    atol: float = 0.0,
) -> list[ComparisonResult]:
    """Compare single-column output with every column in the 2x2 output."""
    fields = (
        "soil_T",
        "soil_q",
        "ust",
        "obl",
        "shf",
        "lhf",
        "ghf",
        "T_skin",
        "sw_in",
        "sw_out",
        "lw_in",
        "lw_out",
        "rnet",
        "bottom_drainage",
    )
    results: list[ComparisonResult] = []

    with nc.Dataset(single_path) as single, nc.Dataset(multi_path) as multi:
        for name in fields:
            single_data = np.asarray(single.variables[name][:])
            multi_data = np.asarray(multi.variables[name][:])

            if multi_data.ndim == single_data.ndim:
                max_abs_diff = float(np.max(np.abs(single_data - multi_data)))
                matches = np.allclose(
                    single_data, multi_data, rtol=rtol, atol=atol
                )
                max_abs_col_diff = None
            elif multi_data.ndim == single_data.ndim + 2:
                ny, nx = multi_data.shape[-2:]
                max_abs_diff = 0.0
                max_abs_col_diff = 0.0
                matches = True
                ref_col = multi_data[..., 0, 0]
                for y in range(ny):
                    for x in range(nx):
                        col_data = multi_data[..., y, x]
                        diff = np.max(np.abs(single_data - col_data))
                        max_abs_diff = max(max_abs_diff, float(diff))
                        if not np.allclose(
                            single_data, col_data, rtol=rtol, atol=atol
                        ):
                            matches = False
                        col_diff = np.max(np.abs(ref_col - col_data))
                        max_abs_col_diff = max(
                            max_abs_col_diff, float(col_diff)
                        )
            else:
                raise ValueError(
                    f"{name} dims {multi_data.shape} do not align with "
                    f"{single_data.shape}."
                )

            results.append(
                ComparisonResult(
                    name=name,
                    matches=bool(matches),
                    max_abs_diff=max_abs_diff,
                    max_abs_col_diff=max_abs_col_diff,
                )
            )

    return results
