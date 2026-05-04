#!/usr/bin/env python
#
# UtahLSM
#
# Copyright (c) 2017–2026 Jeremy A. Gibbs
# Copyright (c) 2017–2026 Rob Stoll
# Copyright (c) 2017–2026 Eric Pardyjak
# Copyright (c) 2017–2026 Pete Willemsen
#
# This file is part of UtahLSM.
#
# This software is free and is distributed under the MIT License.
# See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
#
"""Compare a UtahLSM GABLS3 run against the Cabauw observations.

Plots the model's sensible, latent, and ground heat fluxes (plus u*
and the canopy soil/veg latent partition, when available) against the
merged CESAR surface-flux record in
``cases/gabls3/observations/gabls3_fluxes.nc``.

The model output is expected to be a NetCDF file produced by
``utahlsm_offline.py`` for the GABLS3 case. The observation record
spans 2006-07-01 00:00 UTC through 2006-07-03 06:00 UTC at 1 h cadence.

Example:
-------
::

    python scripts/compare_gabls3_observations.py \
        --model lsm_gabls3_py.nc \
        --obs ../cases/gabls3/observations/gabls3_fluxes.nc \
        --out /tmp/gabls3_compare.png
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

# GABLS3 simulation start time in UTC (matches lsm_namelist.json: utc_year
# 2006, julian_day 183 = 2006-07-02, utc_start = 0). The observation file
# uses ``hours since 2006-07-01 00:00:00`` so both axes are converted to
# datetimes for a shared plot.
MODEL_T0 = np.datetime64("2006-07-02T00:00:00")
OBS_T0 = np.datetime64("2006-07-01T00:00:00")


def _model_times(ds: nc.Dataset) -> np.ndarray:
    """Model time variable (seconds) converted to absolute datetimes."""
    t_sec = np.asarray(ds.variables["time"][:]).astype(float)
    return MODEL_T0 + (t_sec * 1000.0).astype("timedelta64[ms]")


def _obs_times(ds: nc.Dataset) -> np.ndarray:
    """Obs time variable (hours) converted to absolute datetimes."""
    t_hr = np.asarray(ds.variables["time"][:]).astype(float)
    return OBS_T0 + (t_hr * 3600.0 * 1000.0).astype("timedelta64[ms]")


def _overlap(t_model: np.ndarray, t_obs: np.ndarray) -> slice:
    """Return the slice of the obs time axis that overlaps the model run."""
    start = max(t_model[0], t_obs[0])
    end = min(t_model[-1], t_obs[-1])
    idx = np.where((t_obs >= start) & (t_obs <= end))[0]
    if idx.size == 0:
        return slice(0, 0)
    return slice(int(idx[0]), int(idx[-1]) + 1)


def _plot_panel(ax: Any, t_m: np.ndarray, y_m: np.ndarray | None, t_o: np.ndarray, y_o: np.ndarray | None, ylabel: str, title: str, model_label: str,
                partition: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None) -> None:
    """Plot a single flux panel (model vs obs, optionally partitioned)."""
    if y_m is None or y_o is None:
        return
    ax.plot(t_m, y_m, color="#d62728", lw=1.5, label=model_label)
    ax.plot(t_o, y_o, color="k", marker="o", ms=3, lw=0, label="CESAR obs")
    if partition is not None:
        (t_ms, y_soil, y_veg) = partition
        ax.plot(t_ms, y_soil, color="#ff7f0e", lw=1.0, ls="--",
                label="bare-soil evap")
        ax.plot(t_ms, y_veg, color="#2ca02c", lw=1.0, ls="--",
                label="transpiration")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.5)
    ax.axhline(0.0, color="gray", lw=0.5)
    ax.legend(loc="best", fontsize=8)


def _nan_fill(arr: np.ndarray) -> np.ndarray:
    """Convert masked-array fill values to NaN for clean plotting."""
    a = np.asarray(arr, dtype=float)
    mask = np.ma.getmaskarray(arr) if np.ma.isMaskedArray(arr) else None
    if mask is not None and mask.any():
        a = a.copy()
        a[mask] = np.nan
    # CESAR files use -9999 sentinels in some fields.
    a = np.where(a < -999.0, np.nan, a)
    return a


def compare(model_path: Path, obs_path: Path, out_path: Path | None,
            show: bool = False) -> None:
    """Create a multi-panel comparison figure."""
    with nc.Dataset(model_path) as mds, nc.Dataset(obs_path) as ods:
        t_model = _model_times(mds)
        t_obs = _obs_times(ods)
        sl = _overlap(t_model, t_obs)

        def mvar(name: str) -> np.ndarray | None:
            if name not in mds.variables:
                return None
            # Model vars may be (t, y, x) for 2x2 cases. We take column 0.
            v = np.asarray(mds.variables[name][:])
            if v.ndim >= 2:
                v = v.reshape(v.shape[0], -1)[:, 0]
            return _nan_fill(v)

        shf_m = mvar("shf")
        lhf_m = mvar("lhf")
        ghf_m = mvar("ghf")
        ust_m = mvar("ust")
        lhf_soil_m = mvar("lhf_soil")
        lhf_veg_m = mvar("lhf_veg")

        shf_o = _nan_fill(ods.variables["H"][sl])
        lhf_o = _nan_fill(ods.variables["LE"][sl])
        ghf_o = _nan_fill(ods.variables["G0"][sl])
        ust_o = _nan_fill(ods.variables["UST"][sl])
        t_obs_sl = t_obs[sl]

    fig, axes = plt.subplots(4, 1, figsize=(10, 11), sharex=True, squeeze=True)

    partition = None
    if lhf_soil_m is not None and lhf_veg_m is not None:
        partition = (t_model, lhf_soil_m, lhf_veg_m)

    _plot_panel(axes[0], t_model, shf_m, t_obs_sl, shf_o,
                "W m$^{-2}$", "Sensible heat flux", "UtahLSM H")
    _plot_panel(axes[1], t_model, lhf_m, t_obs_sl, lhf_o,
                "W m$^{-2}$", "Latent heat flux", "UtahLSM LE",
                partition=partition)
    _plot_panel(axes[2], t_model, ghf_m, t_obs_sl, ghf_o,
                "W m$^{-2}$", "Ground heat flux", "UtahLSM G0")
    _plot_panel(axes[3], t_model, ust_m, t_obs_sl, ust_o,
                "m s$^{-1}$", "Friction velocity", "UtahLSM u*")

    axes[-1].set_xlabel("Time (UTC)")
    axes[-1].xaxis.set_major_locator(mdates.HourLocator(interval=1))
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H%M"))
    fig.suptitle(f"GABLS3: UtahLSM vs Cabauw observations\n"
                 f"model: {model_path.name}   obs: {obs_path.name}",
                 fontsize=11)
    fig.autofmt_xdate()
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    if out_path is not None:
        fig.savefig(out_path, dpi=150)
        print(f"Saved figure → {out_path}")
    if show:
        plt.show()
    plt.close(fig)

    # Headline stats (overlap window).
    def _stats(label: str, m: np.ndarray | None, o: np.ndarray, t_m: np.ndarray, t_o: np.ndarray) -> None:
        if m is None:
            print(f"  {label:>8}: not in model output")
            return
        # Nearest-neighbor interpolate model onto obs times for bias/RMSE.
        mi: np.ndarray = np.interp(
            (t_o - t_o[0]).astype("timedelta64[s]").astype(float),
            (t_m - t_o[0]).astype("timedelta64[s]").astype(float),
            m,
        )
        diff_raw: np.ndarray = np.asarray(mi - o)
        diff: np.ndarray = np.asarray(diff_raw[np.isfinite(diff_raw)])
        if diff.size == 0:
            print(f"  {label:>8}: no overlap")
            return
        print(f"  {label:>8}: bias={diff.mean():+7.2f}  "
              f"rmse={np.sqrt((diff**2).mean()):6.2f}  "
              f"peak obs={np.nanmax(o):7.2f}  peak model={np.nanmax(np.asarray(mi)):7.2f}")

    print("\nSummary (overlap window only):")
    _stats("H", shf_m, shf_o, t_model, t_obs_sl)
    _stats("LE", lhf_m, lhf_o, t_model, t_obs_sl)
    _stats("G0", ghf_m, ghf_o, t_model, t_obs_sl)
    _stats("u*", ust_m, ust_o, t_model, t_obs_sl)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for GABLS3 comparison script."""
    here = Path(__file__).resolve().parent
    default_model = here.parent / "lsm_gabls3_py.nc"
    default_obs = (here.parent.parent / "cases" / "gabls3"
                   / "observations" / "gabls3_fluxes.nc")
    default_out = here.parent / "gabls3_compare.png"

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--model", type=Path, default=default_model,
                   help="UtahLSM NetCDF output (default: %(default)s)")
    p.add_argument("--obs", type=Path, default=default_obs,
                   help="CESAR merged flux NetCDF (default: %(default)s)")
    p.add_argument("--out", type=Path, default=default_out,
                   help="Output figure path (default: %(default)s)")
    p.add_argument("--show", action="store_true",
                   help="Display the figure interactively in addition "
                        "to saving it.")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    compare(args.model, args.obs, args.out, show=args.show)
