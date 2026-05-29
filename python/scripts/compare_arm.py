#!/usr/bin/env python
#
# UtahLSM
#
# Copyright (c) 2017-2026 Jeremy A. Gibbs
# Copyright (c) 2017-2026 Rob Stoll
# Copyright (c) 2017-2026 Eric Pardyjak
# Copyright (c) 2017-2026 Pete Willemsen
#
# This file is part of UtahLSM.
#
# This software is free and is distributed under the MIT License.
# See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
#
"""Compare a UtahLSM ARM SGP run against ARM observations.

Produces three figures from one command:

1. Surface fluxes: sensible heat, latent heat, and ground heat flux.
2. Soil temperature: timeseries by observation depth plus final profile.
3. Soil moisture: timeseries by observation depth plus final profile.

Example:
-------
::

    python scripts/compare_arm.py \
        --model lsm_arm_py.nc \
        --obs-dir ../cases/arm/observations \
        --out /tmp/arm_compare.png
"""

from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import matplotlib.axes
import matplotlib.dates as mdates
import matplotlib.figure
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

MODEL_T0 = np.datetime64("2017-06-17T18:00:00")

SOIL_DEPTHS: dict[str, float] = {
    "z05": 0.05,
    "z10": 0.10,
    "z20": 0.20,
    "z50": 0.50,
}
SOIL_TEMP_VARS = (
    "soil_temperature_west",
    "soil_temperature_east",
    "soil_temperature_south",
)
SOIL_MOIS_VARS = (
    "soil_specific_water_content_west",
    "soil_specific_water_content_east",
    "soil_specific_water_content_south",
)
TEMP_MIN_SPAN_K = 2.0
MOIS_MIN_SPAN = 0.05


@dataclass(frozen=True)
class ModelData:
    """ARM model output needed by all comparison figures."""

    path: Path
    time: np.ndarray
    soil_z: np.ndarray
    soil_t: np.ndarray | None
    soil_q: np.ndarray | None
    shf: np.ndarray | None
    lhf: np.ndarray | None
    ghf: np.ndarray | None
    lhf_soil: np.ndarray | None
    lhf_veg: np.ndarray | None


@dataclass(frozen=True)
class OutputPaths:
    """Destination paths for the three generated figures."""

    flux: Path
    soil_temperature: Path
    soil_moisture: Path


@dataclass(frozen=True)
class CliArgs:
    """Typed command-line arguments."""

    model: Path
    obs_dir: Path
    surface_flux_obs: Path
    surface_radn_obs: Path
    soil_obs: Path
    out: Path
    show: bool


def _model_times(ds: nc.Dataset) -> np.ndarray:
    """Return model seconds converted to absolute UTC datetimes."""
    t_sec = np.asarray(ds.variables["time"][:]).astype(float)
    return MODEL_T0 + (t_sec * 1000.0).astype("timedelta64[ms]")


def _obs_times(ds: nc.Dataset) -> np.ndarray:
    """Return ARM observation seconds converted to absolute UTC datetimes."""
    t_sec = np.asarray(ds.variables["time"][:]).astype(float)
    units = getattr(ds.variables["time"], "units", "")
    match = re.search(r"since\s+(\d{4}-\d{2}-\d{2})\s+(\d{1,2}:\d{2}:\d{2})", units)
    if match is None:
        raise RuntimeError(f"Cannot parse ARM time units: {units!r}")
    base = np.datetime64(f"{match.group(1)}T{match.group(2)}")
    return base + (t_sec * 1000.0).astype("timedelta64[ms]")


def _nan_fill(values: np.ndarray) -> np.ndarray:
    """Convert common NetCDF fill values and masks to NaN."""
    out = np.asarray(values, dtype=float)
    if np.ma.isMaskedArray(values):
        mask = np.ma.getmaskarray(values)
        if mask.any():
            out = out.copy()
            out[mask] = np.nan
    return np.where(out < -999.0, np.nan, out)


def _clean_variable(
    ds: nc.Dataset,
    name: str,
    obs_slice: slice,
    *,
    min_value: float | None = None,
    max_value: float | None = None,
) -> np.ndarray:
    """Return an ARM variable with missing, QC-flagged, and bad-range values set to NaN."""
    var = ds.variables[name]
    values = _nan_fill(var[obs_slice])
    valid = np.isfinite(values)

    for attr in ("missing_value", "_FillValue"):
        if hasattr(var, attr):
            valid &= values != float(getattr(var, attr))

    if min_value is not None:
        valid &= values >= min_value
    if max_value is not None:
        valid &= values <= max_value

    qc_name = f"qc_{name}"
    if qc_name in ds.variables:
        qc = np.asarray(ds.variables[qc_name][obs_slice])
        valid &= qc == 0

    cleaned = values.astype(float)
    cleaned[~valid] = np.nan
    return cleaned


def _overlap_slice(t_model: np.ndarray, t_obs: np.ndarray) -> slice:
    """Return the slice of observation times overlapping the model run."""
    start = max(t_model[0], t_obs[0])
    end = min(t_model[-1], t_obs[-1])
    idx = np.where((t_obs >= start) & (t_obs <= end))[0]
    if idx.size == 0:
        return slice(0, 0)
    return slice(int(idx[0]), int(idx[-1]) + 1)


def _column0(values: np.ndarray) -> np.ndarray:
    """Collapse optional horizontal grid dimensions by taking column zero."""
    out = np.asarray(values)
    if out.ndim >= 2:
        out = out.reshape(out.shape[0], -1)[:, 0]
    return out


def _profile_column0(values: np.ndarray) -> np.ndarray:
    """Collapse optional horizontal grid dimensions in a time/depth profile."""
    out = np.asarray(values)
    if out.ndim > 2:
        out = out.reshape(out.shape[0], out.shape[1], -1)[:, :, 0]
    return out


def _model_var(ds: nc.Dataset, name: str) -> np.ndarray | None:
    """Return a model variable as a 1-D timeseries when present."""
    if name not in ds.variables:
        return None
    return _nan_fill(_column0(np.asarray(ds.variables[name][:])))


def load_model(model_path: Path) -> ModelData:
    """Load model output once for all comparisons."""
    with nc.Dataset(model_path) as ds:
        soil_z = -np.asarray(ds.variables["soil_z"][:]).astype(float)
        if np.nanmean(soil_z) < 0.0:
            soil_z = -soil_z

        soil_t = None
        soil_q = None
        if "soil_T" in ds.variables:
            soil_t = _nan_fill(_profile_column0(np.asarray(ds.variables["soil_T"][:])))
        if "soil_q" in ds.variables:
            soil_q = _nan_fill(_profile_column0(np.asarray(ds.variables["soil_q"][:])))

        return ModelData(
            path=model_path,
            time=_model_times(ds),
            soil_z=soil_z,
            soil_t=soil_t,
            soil_q=soil_q,
            shf=_model_var(ds, "shf"),
            lhf=_model_var(ds, "lhf"),
            ghf=_model_var(ds, "ghf"),
            lhf_soil=_model_var(ds, "lhf_soil"),
            lhf_veg=_model_var(ds, "lhf_veg"),
        )


def _finite_values(*series: np.ndarray) -> np.ndarray:
    """Return all finite values from the input arrays."""
    finite: list[np.ndarray] = []
    for values in series:
        arr = np.asarray(values, dtype=float)
        mask = np.isfinite(arr)
        if mask.any():
            finite.append(arr[mask])
    if not finite:
        return np.array([], dtype=float)
    return np.concatenate(finite)


def _set_axis_limits(
    ax: matplotlib.axes.Axes,
    axis: str,
    values: np.ndarray,
    min_span: float,
) -> None:
    """Set padded axis limits while enforcing a minimum plotted span."""
    if values.size == 0:
        return
    lower = float(values.min())
    upper = float(values.max())
    span = upper - lower
    if span < min_span:
        center = 0.5 * (lower + upper)
        lower = center - 0.5 * min_span
        upper = center + 0.5 * min_span
        span = min_span
    pad = max(0.08 * span, 0.02 * min_span)
    getattr(ax, f"set_{axis}lim")(lower - pad, upper + pad)


def _format_time_axis(
    ax: matplotlib.axes.Axes,
    t_start: np.datetime64,
    t_end: np.datetime64,
) -> None:
    """Apply compact UTC date ticks to a timeseries axis."""
    locator = cast(Any, mdates.HourLocator)(interval=3)
    x_start = float(
        cast(Any, mdates.date2num)(t_start.astype("datetime64[ms]").astype(object))
    )
    x_end = float(
        cast(Any, mdates.date2num)(t_end.astype("datetime64[ms]").astype(object))
    )
    ax.set_xlim(x_start, x_end)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(cast(Any, mdates.ConciseDateFormatter)(locator))
    ax.tick_params(axis="x", labelrotation=20)


def _series_stats(
    t_model: np.ndarray,
    model_vals: np.ndarray,
    t_obs: np.ndarray,
    obs_vals: np.ndarray,
) -> tuple[float, float] | None:
    """Return bias and RMSE after interpolating model values to obs times."""
    finite = np.isfinite(obs_vals)
    if t_obs.size == 0 or not finite.any():
        return None
    t_obs_f = t_obs[finite]
    obs_f = obs_vals[finite]
    model_interp = np.asarray(
        np.interp(
            (t_obs_f - t_obs_f[0]).astype("timedelta64[s]").astype(float),
            (t_model - t_obs_f[0]).astype("timedelta64[s]").astype(float),
            model_vals,
        )
    )
    diff = np.asarray(model_interp - obs_f)
    diff = np.asarray(diff[np.isfinite(diff)])
    if diff.size == 0:
        return None
    return float(diff.mean()), float(np.sqrt((diff**2).mean()))


def _print_stats(
    label: str,
    model_values: np.ndarray | None,
    obs_values: np.ndarray,
    t_model: np.ndarray,
    t_obs: np.ndarray,
) -> None:
    """Print bias and RMSE summary for one observed quantity."""
    if model_values is None:
        print(f"  {label:>8}: not in model output")
        return
    finite = np.isfinite(obs_values)
    stats = _series_stats(t_model, model_values, t_obs, obs_values)
    if stats is None or not finite.any():
        print(f"  {label:>8}: no overlap")
        return
    t_obs_f = t_obs[finite]
    model_interp = np.asarray(
        np.interp(
            (t_obs_f - t_obs_f[0]).astype("timedelta64[s]").astype(float),
            (t_model - t_obs_f[0]).astype("timedelta64[s]").astype(float),
            model_values,
        )
    )
    obs_f = obs_values[finite]
    bias, rmse = stats
    print(
        f"  {label:>8}: bias={bias:+7.2f}  rmse={rmse:6.2f}  "
        f"min obs={np.nanmin(obs_f):7.2f}  "
        f"min model={np.nanmin(model_interp):7.2f}  "
        f"max obs={np.nanmax(obs_f):7.2f}  "
        f"max model={np.nanmax(model_interp):7.2f}"
    )


def _plot_flux_panel(
    ax: matplotlib.axes.Axes,
    t_model: np.ndarray,
    y_model: np.ndarray | None,
    t_obs: np.ndarray,
    y_obs: np.ndarray,
    ylabel: str,
    title: str,
    model_label: str,
    partition: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None,
) -> None:
    """Plot one surface-flux panel."""
    if y_model is None:
        ax.set_visible(False)
        return
    ax.plot(t_model, y_model, color="#d62728", lw=1.5, label=model_label)
    ax.plot(t_obs, y_obs, color="k", marker="o", ms=3, lw=0, label="ARM obs")
    if partition is not None:
        t_part, y_soil, y_veg = partition
        ax.plot(
            t_part, y_soil, color="#ff7f0e", lw=1.0, ls="--", label="bare-soil evap"
        )
        ax.plot(t_part, y_veg, color="#2ca02c", lw=1.0, ls="--", label="transpiration")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.5)
    ax.axhline(0.0, color="gray", lw=0.5)
    ax.legend(loc="best", fontsize=8)


def compare_fluxes(
    model: ModelData,
    surface_flux_obs_path: Path,
    surface_radn_obs_path: Path,
    out_path: Path | None,
    show: bool,
) -> None:
    """Create the surface-flux comparison figure."""
    with (
        nc.Dataset(surface_flux_obs_path) as fds,
        nc.Dataset(surface_radn_obs_path) as rds,
    ):
        t_flux = _obs_times(fds)
        t_radn = _obs_times(rds)
        flux_slice = _overlap_slice(model.time, t_flux)
        radn_slice = _overlap_slice(model.time, t_radn)
        t_flux_sl = t_flux[flux_slice]
        t_radn_sl = t_radn[radn_slice]
        shf_obs = _clean_variable(fds, "corrected_sensible_heat_flux", flux_slice)
        lhf_obs = _clean_variable(fds, "corrected_latent_heat_flux", flux_slice)
        # ARM surface_soil_heat_flux_avg is positive upward. UtahLSM ghf is positive downward.
        ghf_obs = -_clean_variable(rds, "surface_soil_heat_flux_avg", radn_slice)

    fig, axes = cast(
        tuple[matplotlib.figure.Figure, np.ndarray],
        cast(Any, plt.subplots)(3, 1, figsize=(10, 8.5), sharex=True, squeeze=True),
    )
    axes = np.asarray(axes, dtype=object)

    partition = None
    if model.lhf_soil is not None and model.lhf_veg is not None:
        partition = (model.time, model.lhf_soil, model.lhf_veg)

    _plot_flux_panel(
        axes[0],
        model.time,
        model.shf,
        t_flux_sl,
        shf_obs,
        "W m$^{-2}$",
        "Sensible heat flux",
        "UtahLSM H",
    )
    _plot_flux_panel(
        axes[1],
        model.time,
        model.lhf,
        t_flux_sl,
        lhf_obs,
        "W m$^{-2}$",
        "Latent heat flux",
        "UtahLSM LE",
        partition=partition,
    )
    _plot_flux_panel(
        axes[2],
        model.time,
        model.ghf,
        t_radn_sl,
        ghf_obs,
        "W m$^{-2}$",
        "Ground heat flux",
        "UtahLSM G0",
    )

    axes[-1].set_xlabel("Time (UTC)")
    t_ax_end = max(
        t_flux_sl[-1] if t_flux_sl.size > 0 else model.time[0],
        t_radn_sl[-1] if t_radn_sl.size > 0 else model.time[0],
    )
    _format_time_axis(axes[-1], model.time[0], t_ax_end)
    fig.suptitle(
        "ARM SGP: UtahLSM vs ARM observations\n"
        f"model: {model.path.name}   flux obs: {surface_flux_obs_path.name}   "
        f"radn obs: {surface_radn_obs_path.name}",
        fontsize=11,
    )
    fig.autofmt_xdate()
    fig.tight_layout(rect=(0, 0, 1, 0.95))

    if out_path is not None:
        fig.savefig(out_path, dpi=150)
        print(f"Saved figure -> {out_path}")
    if show:
        plt.show()
    plt.close(fig)

    print("\nSurface-flux summary (overlap window only):")
    _print_stats("H", model.shf, shf_obs, model.time, t_flux_sl)
    _print_stats("LE", model.lhf, lhf_obs, model.time, t_flux_sl)
    _print_stats("G0", model.ghf, ghf_obs, model.time, t_radn_sl)


def _interp_model_to_depth(
    profile: np.ndarray,
    z_model: np.ndarray,
    depth: float,
) -> np.ndarray:
    """Interpolate a model profile to a fixed positive-downward depth."""
    return np.asarray(
        [np.interp(depth, z_model, profile_at_time) for profile_at_time in profile],
        dtype=float,
    )


def _load_arm_soil_series(
    ds: nc.Dataset,
    variable_names: tuple[str, ...],
    obs_slice: slice,
    scale: float,
    offset: float,
) -> tuple[np.ndarray, dict[str, np.ndarray], list[str]]:
    """Average ARM west/east/south soil observations by depth."""
    depths_cm = np.asarray(ds.variables["depth"][:], dtype=float)
    depth_indices = [
        int(np.argmin(np.abs(depths_cm - depth * 100.0)))
        for depth in SOIL_DEPTHS.values()
    ]
    site_values: list[np.ndarray] = []

    for name in variable_names:
        if name not in ds.variables:
            continue
        values = _clean_variable(ds, name, obs_slice)
        site_values.append(values[:, depth_indices] * scale + offset)

    if not site_values:
        return depths_cm[depth_indices] / 100.0, {}, list(SOIL_DEPTHS)

    stacked = np.stack(site_values, axis=0)
    obs_by_label: dict[str, np.ndarray] = {}
    empty_labels: list[str] = []
    for i, label in enumerate(SOIL_DEPTHS):
        values = np.nanmean(stacked[:, :, i], axis=0)
        if np.isfinite(values).any():
            obs_by_label[label] = values
        else:
            empty_labels.append(label)
    return depths_cm[depth_indices] / 100.0, obs_by_label, empty_labels


def _make_panel_grid(
    n_timeseries: int,
) -> tuple[matplotlib.figure.Figure, np.ndarray, matplotlib.axes.Axes, int, int]:
    """Create timeseries panels plus one final-profile panel."""
    total_panels = n_timeseries + 1
    if total_panels <= 1:
        ncols = 1
    elif total_panels <= 4:
        ncols = 2
    elif total_panels <= 9:
        ncols = 3
    else:
        ncols = 4
    nrows = math.ceil(total_panels / ncols)
    fig, axes = cast(
        tuple[matplotlib.figure.Figure, np.ndarray],
        cast(Any, plt.subplots)(
            nrows,
            ncols,
            figsize=(5.1 * ncols, 2.8 * nrows + 0.4),
            squeeze=False,
        ),
    )
    axes = np.asarray(axes, dtype=object)
    axes_flat = axes.ravel()
    for ax in axes_flat[total_panels:]:
        ax.set_visible(False)
    ax_profile = cast(matplotlib.axes.Axes, axes_flat[n_timeseries])
    return fig, axes_flat[:n_timeseries], ax_profile, ncols, total_panels


def _is_bottom_row(panel_idx: int, total_panels: int, ncols: int) -> bool:
    """Return whether panel_idx appears on the bottom row."""
    return panel_idx // ncols == math.ceil(total_panels / ncols) - 1


def _plot_timeseries_panels(
    axes: np.ndarray,
    t_model: np.ndarray,
    model_vals_by_depth: dict[float, np.ndarray],
    t_obs: np.ndarray,
    obs_vals_by_label: dict[str, np.ndarray],
    depth_map: dict[str, float],
    ylabel: str,
    title_prefix: str,
    min_span: float,
    stats_fmt: str,
) -> None:
    """Plot one timeseries subplot per observation depth."""
    for ax, (label, depth) in zip(axes, depth_map.items()):
        if label not in obs_vals_by_label:
            ax.set_visible(False)
            continue
        model_vals = model_vals_by_depth[depth]
        obs_vals = obs_vals_by_label[label]
        ax.plot(
            t_model,
            model_vals,
            color="#d62728",
            lw=1.2,
            label=f"UtahLSM @ {depth * 100:.0f} cm",
        )
        ax.plot(t_obs, obs_vals, color="k", lw=0.8, alpha=0.8, label=f"ARM {label}")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{title_prefix} @ {depth * 100:.0f} cm")
        ax.grid(True, alpha=0.3)
        _set_axis_limits(ax, "y", _finite_values(model_vals, obs_vals), min_span)
        stats = _series_stats(t_model, model_vals, t_obs, obs_vals)
        if stats is not None:
            bias, rmse = stats
            ax.text(
                0.98,
                0.03,
                stats_fmt.format(bias=bias, rmse=rmse),
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                fontsize=8,
                bbox={
                    "facecolor": "white",
                    "alpha": 0.75,
                    "edgecolor": "none",
                    "boxstyle": "round,pad=0.2",
                },
            )
        ax.legend(loc="upper left", fontsize=8, framealpha=0.85)


def _plot_final_profile(
    ax: matplotlib.axes.Axes,
    z_model: np.ndarray,
    model_final: np.ndarray,
    obs_depths: np.ndarray,
    obs_vals_final: np.ndarray,
    xlabel: str,
    title: str,
    min_span: float,
) -> None:
    """Plot the final model and observed profile."""
    ax.plot(
        model_final, z_model, color="#d62728", marker="o", lw=1.5, label="UtahLSM final"
    )
    if obs_depths.size > 0:
        ax.plot(
            obs_vals_final,
            obs_depths,
            color="k",
            marker="s",
            ms=6,
            lw=0,
            label="ARM final",
        )
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    _set_axis_limits(ax, "x", _finite_values(model_final, obs_vals_final), min_span)
    ax.legend(loc="best", fontsize=8)


def compare_soil_quantity(
    model: ModelData,
    obs_path: Path,
    out_path: Path | None,
    show: bool,
    model_profile: np.ndarray | None,
    variable_names: tuple[str, ...],
    scale: float,
    obs_offset: float,
    ylabel: str,
    title_prefix: str,
    profile_title: str,
    min_span: float,
    stats_fmt: str,
    summary_label: str,
) -> None:
    """Create a soil profile and timeseries comparison figure."""
    if model_profile is None or model.soil_z.size == 0:
        print(f"Skipping {summary_label}: required model fields are missing")
        return

    with nc.Dataset(obs_path) as ods:
        t_obs = _obs_times(ods)
        obs_slice = _overlap_slice(model.time, t_obs)
        t_obs_sl = t_obs[obs_slice]
        if t_obs_sl.size == 0:
            print(f"No overlapping {summary_label} samples in {obs_path}")
            return
        obs_depths, obs_vals, empty_labels = _load_arm_soil_series(
            ods, variable_names, obs_slice, scale, obs_offset
        )

    if empty_labels:
        print(f"Skipping empty {summary_label} obs series: " + ", ".join(empty_labels))

    depths_present = {
        label: depth for label, depth in SOIL_DEPTHS.items() if label in obs_vals
    }
    if not depths_present:
        print(f"No finite {summary_label} observations found in {obs_path}")
        return

    model_by_depth = {
        depth: _interp_model_to_depth(model_profile, model.soil_z, depth)
        for depth in depths_present.values()
    }

    fig, axes_ts, ax_profile, ncols, total_panels = _make_panel_grid(
        len(depths_present)
    )
    _plot_timeseries_panels(
        axes_ts,
        model.time,
        model_by_depth,
        t_obs_sl,
        obs_vals,
        depths_present,
        ylabel,
        title_prefix,
        min_span,
        stats_fmt,
    )
    for i, ax in enumerate(axes_ts):
        _format_time_axis(ax, model.time[0], t_obs_sl[-1])
        if _is_bottom_row(i, total_panels, ncols):
            ax.set_xlabel("Time (UTC)")

    obs_final = np.array(
        [
            obs_vals[var][-1] if np.isfinite(obs_vals[var][-1]) else np.nan
            for var in depths_present
        ]
    )
    mask = np.isfinite(obs_final)
    _plot_final_profile(
        ax_profile,
        model.soil_z,
        model_profile[-1],
        obs_depths[mask],
        obs_final[mask],
        ylabel,
        profile_title,
        min_span,
    )

    fig.suptitle(
        f"ARM SGP {summary_label}: UtahLSM vs ARM observations\n"
        f"model: {model.path.name}   obs: {obs_path.name}",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))

    if out_path is not None:
        fig.savefig(out_path, dpi=150)
        print(f"Saved figure -> {out_path}")
    if show:
        plt.show()
    plt.close(fig)

    print(f"\n{summary_label.capitalize()} bias vs obs (overlap window):")
    for var, depth in depths_present.items():
        stats = _series_stats(
            model.time, model_by_depth[depth], t_obs_sl, obs_vals[var]
        )
        if stats is not None:
            bias, rmse = stats
            print(
                f"  {var:>5} ({depth * 100:4.0f} cm): "
                + stats_fmt.format(bias=bias, rmse=rmse).replace("\n", "  ")
            )


def compare_soil_temperature(
    model: ModelData,
    obs_path: Path,
    out_path: Path | None,
    show: bool,
) -> None:
    """Create the soil-temperature comparison figure."""
    compare_soil_quantity(
        model=model,
        obs_path=obs_path,
        out_path=out_path,
        show=show,
        model_profile=model.soil_t,
        variable_names=SOIL_TEMP_VARS,
        scale=1.0,
        obs_offset=273.15,
        ylabel="T [K]",
        title_prefix="Soil temperature",
        profile_title="Final soil temperature profile",
        min_span=TEMP_MIN_SPAN_K,
        stats_fmt="bias={bias:+.2f} K\nrmse={rmse:.2f} K",
        summary_label="soil temperature",
    )


def compare_soil_moisture(
    model: ModelData,
    obs_path: Path,
    out_path: Path | None,
    show: bool,
) -> None:
    """Create the soil-moisture comparison figure."""
    compare_soil_quantity(
        model=model,
        obs_path=obs_path,
        out_path=out_path,
        show=show,
        model_profile=model.soil_q,
        variable_names=SOIL_MOIS_VARS,
        scale=0.01,
        obs_offset=0.0,
        ylabel=r"$\theta$ [m$^3$/m$^3$]",
        title_prefix="Soil moisture",
        profile_title="Final soil moisture profile",
        min_span=MOIS_MIN_SPAN,
        stats_fmt="bias={bias:+.3f}\nrmse={rmse:.3f}",
        summary_label="soil moisture",
    )


def _output_paths(out_arg: Path) -> OutputPaths:
    """Return all figure paths from an output directory or file-like prefix."""
    if out_arg.suffix:
        out_dir = out_arg.parent if out_arg.parent != Path("") else Path(".")
        prefix = f"{out_arg.stem}_"
    else:
        out_dir = out_arg
        prefix = ""
    out_dir.mkdir(parents=True, exist_ok=True)
    return OutputPaths(
        flux=out_dir / f"{prefix}arm_surf_flux.png",
        soil_temperature=out_dir / f"{prefix}arm_soil_temp.png",
        soil_moisture=out_dir / f"{prefix}arm_soil_mois.png",
    )


def parse_args() -> CliArgs:
    """Parse command-line arguments."""
    here = Path(__file__).resolve().parent
    default_model = here.parent / "lsm_arm_py.nc"
    default_obs_dir = here.parent.parent / "cases" / "arm" / "observations"

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--model",
        type=Path,
        default=default_model,
        help="UtahLSM NetCDF output (default: %(default)s)",
    )
    parser.add_argument(
        "--obs-dir",
        type=Path,
        default=default_obs_dir,
        help="Directory with ARM observation files (default: %(default)s)",
    )
    parser.add_argument(
        "--surface-flux-obs",
        type=Path,
        default=None,
        help="ARM surface flux NetCDF (default: OBS_DIR/arm_surf_flux.nc)",
    )
    parser.add_argument(
        "--surface-radn-obs",
        type=Path,
        default=None,
        help="ARM surface radiation NetCDF (default: OBS_DIR/arm_surf_radn.nc)",
    )
    parser.add_argument(
        "--soil-obs",
        type=Path,
        default=None,
        help="ARM soil observation NetCDF (default: OBS_DIR/arm_soil_data.nc)",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=here.parent,
        help="Output directory, or file-like path used as a figure prefix (default: %(default)s)",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display figures interactively in addition to saving them.",
    )

    raw = parser.parse_args()
    obs_dir = cast(Path, raw.obs_dir)
    return CliArgs(
        model=cast(Path, raw.model),
        obs_dir=obs_dir,
        surface_flux_obs=cast(Path | None, raw.surface_flux_obs)
        or obs_dir / "arm_surf_flux.nc",
        surface_radn_obs=cast(Path | None, raw.surface_radn_obs)
        or obs_dir / "arm_surf_radn.nc",
        soil_obs=cast(Path | None, raw.soil_obs) or obs_dir / "arm_soil_data.nc",
        out=cast(Path, raw.out),
        show=cast(bool, raw.show),
    )


def main() -> None:
    """Run all requested ARM comparisons."""
    args = parse_args()
    outputs = _output_paths(args.out)
    model = load_model(args.model)

    if args.surface_flux_obs.exists() and args.surface_radn_obs.exists():
        compare_fluxes(
            model,
            args.surface_flux_obs,
            args.surface_radn_obs,
            outputs.flux,
            args.show,
        )
    else:
        print(
            "Skipping surface fluxes: "
            f"{args.surface_flux_obs} or {args.surface_radn_obs} not found"
        )

    if args.soil_obs.exists():
        compare_soil_temperature(
            model,
            args.soil_obs,
            outputs.soil_temperature,
            args.show,
        )
        compare_soil_moisture(
            model,
            args.soil_obs,
            outputs.soil_moisture,
            args.show,
        )
    else:
        print(f"Skipping soil comparisons: {args.soil_obs} not found")


if __name__ == "__main__":
    main()
