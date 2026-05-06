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
"""Compare a UtahLSM GABLS3 run against Cabauw observations.

Produces three figures from one command:

1. Surface fluxes: sensible heat, latent heat, ground heat flux, and u*.
2. Soil temperature: timeseries by observation depth plus final profile.
3. Soil moisture: timeseries by observation depth plus final profile.

Example:
-------
::

    python scripts/compare_gabls3.py \
        --model lsm_gabls3_py.nc \
        --obs-dir ../cases/gabls3/observations \
        --out /tmp/gabls3_compare.png
"""
from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import matplotlib.axes
import matplotlib.dates as mdates
import matplotlib.figure
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

MODEL_T0 = np.datetime64("2006-07-02T00:00:00")
OBS_T0 = np.datetime64("2006-07-01T00:00:00")

TEMP_DEPTHS: dict[str, float] = {
    "TS00": 0.00,
    "TS02": 0.02,
    "TS04": 0.04,
    "TS06": 0.06,
    "TS08": 0.08,
    "TS12": 0.12,
    "TS20": 0.20,
    "TS30": 0.30,
    "TS50": 0.50,
}

MOIS_DEPTHS_TDR: dict[str, float] = {
    "TH03": 0.03,
    "TH08": 0.08,
    "TH20": 0.20,
}
MOIS_DEPTHS_EB: dict[str, float] = {
    "TH05": 0.05,
    "TH19": 0.19,
    "TH33": 0.33,
    "TH40": 0.40,
    "TH56": 0.56,
}

TEMP_MIN_SPAN_K = 2.0
MOIS_MIN_SPAN = 0.05


@dataclass(frozen=True)
class ModelData:
    """GABLS3 model output needed by all comparison figures."""

    path: Path
    time: np.ndarray
    soil_z: np.ndarray
    soil_t: np.ndarray | None
    soil_q: np.ndarray | None
    shf: np.ndarray | None
    lhf: np.ndarray | None
    ghf: np.ndarray | None
    ust: np.ndarray | None
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
    soil_heat_obs: Path
    soil_moisture_obs: Path
    out: Path
    show: bool


def _model_times(ds: nc.Dataset) -> np.ndarray:
    """Return model seconds converted to absolute UTC datetimes."""
    t_sec = np.asarray(ds.variables["time"][:]).astype(float)
    return MODEL_T0 + (t_sec * 1000.0).astype("timedelta64[ms]")


def _obs_times(ds: nc.Dataset) -> np.ndarray:
    """Return observation hours converted to absolute UTC datetimes."""
    t_hr = np.asarray(ds.variables["time"][:]).astype(float)
    return OBS_T0 + (t_hr * 3600.0 * 1000.0).astype("timedelta64[ms]")


def _nan_fill(values: np.ndarray) -> np.ndarray:
    """Convert common NetCDF fill values and masks to NaN."""
    out = np.asarray(values, dtype=float)
    if np.ma.isMaskedArray(values):
        mask = np.ma.getmaskarray(values)
        if mask.any():
            out = out.copy()
            out[mask] = np.nan
    return np.where(out < -999.0, np.nan, out)


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
        soil_z: np.ndarray
        soil_t: np.ndarray | None = None
        soil_q: np.ndarray | None = None

        t_model = _model_times(ds)
        if "soil_z" in ds.variables:
            # Model soil_z is negative downward; comparison plots use positive down.
            soil_z = -np.asarray(ds.variables["soil_z"][:]).astype(float)
        else:
            soil_z = np.array([], dtype=float)

        if "soil_T" in ds.variables:
            soil_t = _nan_fill(_profile_column0(np.asarray(ds.variables["soil_T"][:])))
        if "soil_q" in ds.variables:
            soil_q = _nan_fill(_profile_column0(np.asarray(ds.variables["soil_q"][:])))

        return ModelData(
            path=model_path,
            time=t_model,
            soil_z=soil_z,
            soil_t=soil_t,
            soil_q=soil_q,
            shf=_model_var(ds, "shf"),
            lhf=_model_var(ds, "lhf"),
            ghf=_model_var(ds, "ghf"),
            ust=_model_var(ds, "ust"),
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
    locator = cast(Any, mdates.AutoDateLocator)(minticks=4, maxticks=7)
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
    if t_obs.size == 0:
        return None
    model_interp = np.asarray(
        np.interp(
            (t_obs - t_obs[0]).astype("timedelta64[s]").astype(float),
            (t_model - t_obs[0]).astype("timedelta64[s]").astype(float),
            model_vals,
        )
    )
    diff = np.asarray(model_interp - obs_vals)
    diff = np.asarray(diff[np.isfinite(diff)])
    if diff.size == 0:
        return None
    return float(diff.mean()), float(np.sqrt((diff**2).mean()))


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
    ax.plot(t_obs, y_obs, color="k", marker="o", ms=3, lw=0, label="CESAR obs")
    if partition is not None:
        t_part, y_soil, y_veg = partition
        ax.plot(t_part, y_soil, color="#ff7f0e", lw=1.0, ls="--", label="bare-soil evap")
        ax.plot(t_part, y_veg, color="#2ca02c", lw=1.0, ls="--", label="transpiration")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.5)
    ax.axhline(0.0, color="gray", lw=0.5)
    ax.legend(loc="best", fontsize=8)


def compare_fluxes(
    model: ModelData,
    surface_flux_obs_path: Path,
    soil_heat_obs_path: Path,
    out_path: Path | None,
    show: bool,
) -> None:
    """Create the surface-flux comparison figure."""
    with nc.Dataset(surface_flux_obs_path) as ods, nc.Dataset(soil_heat_obs_path) as sds:
        t_obs = _obs_times(ods)
        t_soil = _obs_times(sds)
        obs_slice = _overlap_slice(model.time, t_obs)
        soil_slice = _overlap_slice(model.time, t_soil)
        t_obs_sl = t_obs[obs_slice]
        t_soil_sl = t_soil[soil_slice]
        shf_obs = _nan_fill(ods.variables["H"][obs_slice])
        lhf_obs = _nan_fill(ods.variables["LE"][obs_slice])
        ust_obs = _nan_fill(ods.variables["UST"][obs_slice])
        ghf_obs = _nan_fill(sds.variables["FG0"][soil_slice])

    fig, axes = cast(
        tuple[matplotlib.figure.Figure, np.ndarray],
        cast(Any, plt.subplots)(4, 1, figsize=(10, 11), sharex=True, squeeze=True),
    )
    axes = np.asarray(axes, dtype=object)

    partition = None
    if model.lhf_soil is not None and model.lhf_veg is not None:
        partition = (model.time, model.lhf_soil, model.lhf_veg)

    _plot_flux_panel(
        axes[0], model.time, model.shf, t_obs_sl, shf_obs,
        "W m$^{-2}$", "Sensible heat flux", "UtahLSM H",
    )
    _plot_flux_panel(
        axes[1], model.time, model.lhf, t_obs_sl, lhf_obs,
        "W m$^{-2}$", "Latent heat flux", "UtahLSM LE", partition=partition,
    )
    _plot_flux_panel(
        axes[2], model.time, model.ghf, t_soil_sl, ghf_obs,
        "W m$^{-2}$", "Ground heat flux", "UtahLSM G0",
    )
    _plot_flux_panel(
        axes[3], model.time, model.ust, t_obs_sl, ust_obs,
        "m s$^{-1}$", "Friction velocity", "UtahLSM u*",
    )

    axes[-1].set_xlabel("Time (UTC)")
    axes[-1].xaxis.set_major_locator(cast(Any, mdates.HourLocator)(interval=1))
    axes[-1].xaxis.set_major_formatter(cast(Any, mdates.DateFormatter)("%m-%d %H%M"))
    fig.suptitle(
        "GABLS3: UtahLSM vs Cabauw observations\n"
        f"model: {model.path.name}   flux obs: {surface_flux_obs_path.name}   "
        f"soil heat obs: {soil_heat_obs_path.name}",
        fontsize=11,
    )
    fig.autofmt_xdate()
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    if out_path is not None:
        fig.savefig(out_path, dpi=150)
        print(f"Saved figure -> {out_path}")
    if show:
        plt.show()
    plt.close(fig)

    print("\nSurface-flux summary (overlap window only):")
    _print_stats("H", model.shf, shf_obs, model.time, t_obs_sl)
    _print_stats("LE", model.lhf, lhf_obs, model.time, t_obs_sl)
    _print_stats("FG0", model.ghf, ghf_obs, model.time, t_soil_sl)
    _print_stats("u*", model.ust, ust_obs, model.time, t_obs_sl)


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
    stats = _series_stats(t_model, model_values, t_obs, obs_values)
    if stats is None:
        print(f"  {label:>8}: no overlap")
        return
    model_interp = np.asarray(
        np.interp(
            (t_obs - t_obs[0]).astype("timedelta64[s]").astype(float),
            (t_model - t_obs[0]).astype("timedelta64[s]").astype(float),
            model_values,
        )
    )
    bias, rmse = stats
    print(
        f"  {label:>8}: bias={bias:+7.2f}  rmse={rmse:6.2f}  "
        f"min obs={np.nanmin(obs_values):7.2f}  "
        f"min model={np.nanmin(model_interp):7.2f}  "
        f"max obs={np.nanmax(obs_values):7.2f}  "
        f"max model={np.nanmax(model_interp):7.2f}"
    )


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


def _load_obs_series(
    ds: nc.Dataset,
    labels: dict[str, float],
    obs_slice: slice,
    offset: float = 0.0,
) -> tuple[dict[str, np.ndarray], list[str]]:
    """Load finite observation series by variable label."""
    obs_vals: dict[str, np.ndarray] = {}
    empty_labels: list[str] = []
    for label in labels:
        if label not in ds.variables:
            continue
        values = _nan_fill(ds.variables[label][obs_slice]) + offset
        if np.isfinite(values).any():
            obs_vals[label] = values
        else:
            empty_labels.append(label)
    return obs_vals, empty_labels


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
        ax.plot(t_obs, obs_vals, color="k", lw=0.8, alpha=0.8, label=f"obs {label}")
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
    ax.plot(model_final, z_model, color="#d62728", marker="o", lw=1.5, label="UtahLSM final")
    if obs_depths.size > 0:
        ax.plot(obs_vals_final, obs_depths, color="k", marker="s", ms=6, lw=0, label="obs final")
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
    depth_map: dict[str, float],
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
        obs_vals, empty_labels = _load_obs_series(ods, depth_map, obs_slice, obs_offset)

    if empty_labels:
        print(f"Skipping empty {summary_label} obs series: " + ", ".join(empty_labels))

    depths_present = {label: depth for label, depth in depth_map.items() if label in obs_vals}
    if not depths_present:
        print(f"No finite {summary_label} observations found in {obs_path}")
        return

    model_by_depth = {
        depth: _interp_model_to_depth(model_profile, model.soil_z, depth)
        for depth in depths_present.values()
    }

    fig, axes_ts, ax_profile, ncols, total_panels = _make_panel_grid(len(depths_present))
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
        _format_time_axis(ax, model.time[0], model.time[-1])
        if _is_bottom_row(i, total_panels, ncols):
            ax.set_xlabel("Time (UTC)")

    depths_arr = np.array(list(depths_present.values()))
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
        depths_arr[mask],
        obs_final[mask],
        ylabel,
        profile_title,
        min_span,
    )

    fig.suptitle(
        f"GABLS3 {summary_label}: UtahLSM vs Cabauw\n"
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
        stats = _series_stats(model.time, model_by_depth[depth], t_obs_sl, obs_vals[var])
        if stats is not None:
            bias, rmse = stats
            print(f"  {var:>5} ({depth * 100:4.0f} cm): " + stats_fmt.format(bias=bias, rmse=rmse).replace("\n", "  "))


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
        depth_map=TEMP_DEPTHS,
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
        depth_map={**MOIS_DEPTHS_TDR, **MOIS_DEPTHS_EB},
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
        flux=out_dir / f"{prefix}gabls3_surf_flux.png",
        soil_temperature=out_dir / f"{prefix}gabls3_soil_temp.png",
        soil_moisture=out_dir / f"{prefix}gabls3_soil_mois.png",
    )


def parse_args() -> CliArgs:
    """Parse command-line arguments."""
    here = Path(__file__).resolve().parent
    default_model = here.parent / "lsm_gabls3_py.nc"
    default_obs_dir = here.parent.parent / "cases" / "gabls3" / "observations"

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=default_model, help="UtahLSM NetCDF output (default: %(default)s)")
    parser.add_argument("--obs-dir", type=Path, default=default_obs_dir, help="Directory with CESAR observation files (default: %(default)s)")
    parser.add_argument("--surface-flux-obs", type=Path, default=None, help="CESAR merged flux NetCDF (default: OBS_DIR/gabls3_surf_flux.nc)")
    parser.add_argument("--soil-heat-obs", type=Path, default=None, help="CESAR soil heat NetCDF (default: OBS_DIR/gabls3_soil_heat.nc)")
    parser.add_argument("--soil-moisture-obs", type=Path, default=None, help="CESAR soil moisture NetCDF (default: OBS_DIR/gabls3_soil_mois_thc.nc)")
    parser.add_argument("--out", type=Path, default=here.parent, help="Output directory, or file-like path used as a figure prefix (default: %(default)s)")
    parser.add_argument("--show", action="store_true", help="Display figures interactively in addition to saving them.")

    raw = parser.parse_args()
    obs_dir = cast(Path, raw.obs_dir)
    return CliArgs(
        model=cast(Path, raw.model),
        obs_dir=obs_dir,
        surface_flux_obs=cast(Path | None, raw.surface_flux_obs) or obs_dir / "gabls3_surf_flux.nc",
        soil_heat_obs=cast(Path | None, raw.soil_heat_obs) or obs_dir / "gabls3_soil_heat.nc",
        soil_moisture_obs=cast(Path | None, raw.soil_moisture_obs) or obs_dir / "gabls3_soil_mois_thc.nc",
        out=cast(Path, raw.out),
        show=cast(bool, raw.show),
    )


def main() -> None:
    """Run all requested GABLS3 comparisons."""
    args = parse_args()
    outputs = _output_paths(args.out)
    model = load_model(args.model)

    if args.surface_flux_obs.exists() and args.soil_heat_obs.exists():
        compare_fluxes(
            model,
            args.surface_flux_obs,
            args.soil_heat_obs,
            outputs.flux,
            args.show,
        )
    else:
        print(
            "Skipping surface fluxes: "
            f"{args.surface_flux_obs} or {args.soil_heat_obs} not found"
        )

    if args.soil_heat_obs.exists():
        compare_soil_temperature(
            model,
            args.soil_heat_obs,
            outputs.soil_temperature,
            args.show,
        )
    else:
        print(f"Skipping soil temperature: {args.soil_heat_obs} not found")

    if args.soil_moisture_obs.exists():
        compare_soil_moisture(
            model,
            args.soil_moisture_obs,
            outputs.soil_moisture,
            args.show,
        )
    else:
        print(f"Skipping soil moisture: {args.soil_moisture_obs} not found")


if __name__ == "__main__":
    main()
