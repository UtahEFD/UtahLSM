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
"""Compare a UtahLSM GABLS3 run against the Cabauw soil observations.

Produces two figures:

1. Soil temperature — model profile vs CESAR thermistor chain
   (``TS00`` / ``TS02`` … ``TS50``). Plots timeseries at each obs depth
   with the model interpolated to that depth, plus a final-time
   profile comparison panel.

2. Soil moisture — model profile vs CESAR soil-water probes
   (Campbell-calibrated TDR ``TH03/TH08/TH20`` and the EB-field probes
   at 5/19/33/40/56 cm). Same layout as the temperature plot.

Example:
-------
::

    python scripts/compare_gabls3_soil.py \
        --model lsm_gabls3_py.nc \
        --obs-dir ../cases/gabls3/observations \
        --out-dir /tmp
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import cast

import matplotlib.axes
import matplotlib.dates as mdates
import matplotlib.figure
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

MODEL_T0 = np.datetime64("2006-07-02T00:00:00")
OBS_T0 = np.datetime64("2006-07-01T00:00:00")

# CESAR soil-temperature depths (metres, positive downward).
TEMP_DEPTHS = {
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

# CESAR soil-water depths. The Campbell-calibrated TDR set (TH03/08/20)
# is more reliable than the EB-field probes; we plot both for context.
MOIS_DEPTHS_TDR = {
    "TH03": 0.03,
    "TH08": 0.08,
    "TH20": 0.20,
}
MOIS_DEPTHS_EB = {
    "TH05": 0.05,
    "TH19": 0.19,
    "TH33": 0.33,
    "TH40": 0.40,
    "TH56": 0.56,
}

TEMP_MIN_SPAN_K = 2.0
MOIS_MIN_SPAN = 0.05


def _model_times(ds: nc.Dataset) -> np.ndarray:
    t_sec = np.asarray(ds.variables["time"][:]).astype(float)  # type: ignore[arg-type]
    return MODEL_T0 + (t_sec * 1000.0).astype("timedelta64[ms]")


def _obs_times(ds: nc.Dataset) -> np.ndarray:
    t_hr = np.asarray(ds.variables["time"][:]).astype(float)  # type: ignore[arg-type]
    return OBS_T0 + (t_hr * 3600.0 * 1000.0).astype("timedelta64[ms]")


def _nan_fill(arr: np.ndarray) -> np.ndarray:
    a = np.asarray(arr, dtype=float)
    if np.ma.isMaskedArray(arr):
        mask = np.ma.getmaskarray(arr)
        if mask.any():
            a = a.copy()
            a[mask] = np.nan
    return np.where(a < -9000.0, np.nan, a)


def _overlap_slice(t_model: np.ndarray, t_obs: np.ndarray) -> slice:
    start = max(t_model[0], t_obs[0])
    end = min(t_model[-1], t_obs[-1])
    idx = np.where((t_obs >= start) & (t_obs <= end))[0]
    if idx.size == 0:
        return slice(0, 0)
    return slice(int(idx[0]), int(idx[-1]) + 1)


def _interp_model_to_depth(
    profile: np.ndarray, z_model: np.ndarray, depth: float
) -> np.ndarray:
    """Linear interp of a (ntime, nz) model profile to a fixed depth [m].

    ``z_model`` is in metres, positive downward, z_model[0] = 0.
    """
    # np.interp requires ascending xp; z_model is already ascending.
    out = np.empty(profile.shape[0])
    for i in range(profile.shape[0]):
        out[i] = np.interp(depth, z_model, profile[i])
    return out


def _load_model(model_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with nc.Dataset(model_path) as mds:
        t_model = _model_times(mds)
        # soil_z in model is negative-downward (z[0]=0, z[1]<0, ...).
        z_neg = np.asarray(mds.variables["soil_z"][:]).astype(float)  # type: ignore[arg-type]
        z_model = -z_neg  # positive-downward metres
        soil_T = np.asarray(mds.variables["soil_T"][:]).astype(float)  # type: ignore[arg-type]
        soil_q = np.asarray(mds.variables["soil_q"][:]).astype(float)  # type: ignore[arg-type]
        # Collapse any (t, z, y, x) to (t, z) by taking column 0.
        if soil_T.ndim > 2:
            soil_T = soil_T.reshape(soil_T.shape[0], soil_T.shape[1], -1)[:, :, 0]
        if soil_q.ndim > 2:
            soil_q = soil_q.reshape(soil_q.shape[0], soil_q.shape[1], -1)[:, :, 0]
    return t_model, z_model, soil_T, soil_q


def _load_obs_series(
    ds: nc.Dataset, labels: dict[str, float], sl: slice, offset: float = 0.0
) -> tuple[dict[str, np.ndarray], list[str]]:
    obs_vals: dict[str, np.ndarray] = {}
    empty_labels: list[str] = []
    for label in labels:
        if label not in ds.variables:
            continue
        values = _nan_fill(ds.variables[label][sl]) + offset
        if np.isfinite(values).any():
            obs_vals[label] = values
        else:
            empty_labels.append(label)
    return obs_vals, empty_labels


def _finite_values(*series: np.ndarray) -> np.ndarray:
    finite: list[np.ndarray] = []
    for values in series:
        arr = np.asarray(values, dtype=float)
        mask = np.isfinite(arr)
        if mask.any():
            finite.append(arr[mask])
    if not finite:
        return np.array([], dtype=float)
    return np.concatenate(finite)


def _set_axis_limits(ax: matplotlib.axes.Axes, axis: str, values: np.ndarray, min_span: float) -> None:
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


def _make_panel_grid(n_timeseries: int) -> tuple[matplotlib.figure.Figure, np.ndarray, matplotlib.axes.Axes, int, int]:
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
    fig: matplotlib.figure.Figure
    axes: np.ndarray
    fig, axes = plt.subplots(  # type: ignore[misc]
        nrows,
        ncols,
        figsize=(5.1 * ncols, 2.8 * nrows + 0.4),
        squeeze=False,
    )
    axes_flat = axes.ravel()
    for ax in axes_flat[total_panels:]:
        ax.set_visible(False)
    ax_profile = cast(matplotlib.axes.Axes, axes_flat[n_timeseries])
    return fig, axes_flat[:n_timeseries], ax_profile, ncols, total_panels


def _format_time_axis(ax: matplotlib.axes.Axes, t_start: np.datetime64, t_end: np.datetime64) -> None:
    locator = mdates.AutoDateLocator(minticks=4, maxticks=7)
    ax.set_xlim(t_start, t_end)  # type: ignore[arg-type]
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
    ax.tick_params(axis="x", labelrotation=20)  # type: ignore[call-arg]


def _is_bottom_row(panel_idx: int, total_panels: int, ncols: int) -> bool:
    return panel_idx // ncols == math.ceil(total_panels / ncols) - 1


def _series_stats(
    t_model: np.ndarray,
    model_vals: np.ndarray,
    t_obs: np.ndarray,
    obs_vals: np.ndarray,
) -> tuple[float, float] | None:
    if t_obs.size == 0:
        return None
    model_interp: np.ndarray = np.asarray(np.interp(  # type: ignore[arg-type]
        (t_obs - t_obs[0]).astype("timedelta64[s]").astype(float),  # type: ignore[misc]
        (t_model - t_obs[0]).astype("timedelta64[s]").astype(float),  # type: ignore[misc]
        model_vals,
    ))
    diff: np.ndarray = np.asarray(model_interp - obs_vals)  # type: ignore[arg-type]
    diff = np.asarray(diff[np.isfinite(diff)])  # type: ignore[arg-type]
    if diff.size == 0:
        return None
    return float(diff.mean()), float(np.sqrt((diff**2).mean()))


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
    """One subplot per obs depth."""
    for ax, (label, depth) in zip(axes, depth_map.items()):
        if label not in obs_vals_by_label:
            ax.set_visible(False)
            continue
        model_vals = model_vals_by_depth[depth]
        obs_vals = obs_vals_by_label[label]
        ax.plot(t_model, model_vals, color="#d62728",  # type: ignore[arg-type]
                lw=1.2, label=f"UtahLSM @ {depth*100:.0f} cm")
        ax.plot(t_obs, obs_vals, color="k", lw=0.8,  # type: ignore[arg-type]
                alpha=0.8, label=f"obs {label}")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{title_prefix} @ {depth*100:.0f} cm")  # type: ignore[misc]
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
    ax.plot(model_final, z_model, color="#d62728", marker="o",  # type: ignore[arg-type]
            lw=1.5, label="UtahLSM final")
    if obs_depths.size > 0:
        ax.plot(obs_vals_final, obs_depths, color="k", marker="s",  # type: ignore[arg-type]
                ms=6, lw=0, label="obs final")
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)  # type: ignore[misc]
    ax.set_title(title)  # type: ignore[misc]
    ax.grid(True, alpha=0.3) # type: ignore[misc]
    _set_axis_limits(ax, "x", _finite_values(model_final, obs_vals_final), min_span)
    ax.legend(loc="best", fontsize=8) # type: ignore[misc]


def compare_soil_temperature(
    model_path: Path, obs_path: Path, out_path: Path | None, show: bool
) -> None:
    """Compare modelled and observed soil temperature profiles."""
    t_model, z_model, soil_T, _ = _load_model(model_path)

    with nc.Dataset(obs_path) as ods:
        t_obs = _obs_times(ods)
        sl = _overlap_slice(t_model, t_obs)
        t_obs_sl = t_obs[sl]
        if t_obs_sl.size == 0:
            print(f"No overlapping soil temperature samples in {obs_path}")
            return
        obs_vals, empty_labels = _load_obs_series(
            ods, TEMP_DEPTHS, sl, offset=273.15
        )

    if empty_labels:
        print("Skipping empty temperature obs series: "
              + ", ".join(empty_labels))

    # Interpolate model to each obs depth.
    depths_present = {k: v for k, v in TEMP_DEPTHS.items() if k in obs_vals}
    if not depths_present:
        print(f"No finite soil temperature observations found in {obs_path}")
        return
    model_by_depth = {d: _interp_model_to_depth(soil_T, z_model, d)
                      for d in depths_present.values()}

    n = len(depths_present)
    fig, axes_ts, ax_profile, ncols, total_panels = _make_panel_grid(n)
    _plot_timeseries_panels(
        axes_ts, t_model, model_by_depth, t_obs_sl, obs_vals,
        depths_present, "T [K]", "Soil temperature", TEMP_MIN_SPAN_K,
        "bias={bias:+.2f} K\nrmse={rmse:.2f} K",
    )
    for i, ax in enumerate(axes_ts):
        _format_time_axis(ax, t_model[0], t_model[-1])
        if _is_bottom_row(i, total_panels, ncols):
            ax.set_xlabel("Time (UTC)")  # type: ignore[misc]

    # Final-time profile comparison (use the last overlapping obs sample).
    depths_arr = np.array(list(depths_present.values()))
    obs_final = np.array([
        obs_vals[var][-1] if np.isfinite(obs_vals[var][-1])
        else np.nan
        for var in depths_present
    ])
    mask = np.isfinite(obs_final)
    _plot_final_profile(
        ax_profile, z_model, soil_T[-1], depths_arr[mask], obs_final[mask],
        "T [K]", "Final soil temperature profile", TEMP_MIN_SPAN_K,
    )

    fig.suptitle(  # type: ignore[misc]
        f"GABLS3 soil temperature: UtahLSM vs Cabauw\n"
        f"model: {model_path.name}   obs: {obs_path.name}",
        fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))

    if out_path is not None:
        fig.savefig(out_path, dpi=150)  # type: ignore[misc]
        print(f"Saved figure → {out_path}")
    if show:
        plt.show()  # type: ignore[misc]
    plt.close(fig)

    print("\nSoil temperature bias vs obs (K, overlap window):")
    for var, depth in depths_present.items():
        stats = _series_stats(t_model, model_by_depth[depth], t_obs_sl, obs_vals[var])
        if stats is not None:
            bias, rmse = stats
            print(f"  {var:>5} ({depth*100:4.0f} cm): "
                  f"bias={bias:+6.2f}  rmse={rmse:5.2f}")


def compare_soil_moisture(
    model_path: Path, obs_path: Path, out_path: Path | None, show: bool
) -> None:
    """Compare modelled and observed soil moisture profiles."""
    t_model, z_model, _, soil_q = _load_model(model_path)

    with nc.Dataset(obs_path) as ods:
        t_obs = _obs_times(ods)
        sl = _overlap_slice(t_model, t_obs)
        t_obs_sl = t_obs[sl]
        if t_obs_sl.size == 0:
            print(f"No overlapping soil moisture samples in {obs_path}")
            return
        obs_vals, empty_labels = _load_obs_series(
            ods, {**MOIS_DEPTHS_TDR, **MOIS_DEPTHS_EB}, sl
        )

    if empty_labels:
        print("Skipping empty moisture obs series: " + ", ".join(empty_labels))

    depths_present = {**MOIS_DEPTHS_TDR, **MOIS_DEPTHS_EB}
    depths_present = {k: v for k, v in depths_present.items()
                      if k in obs_vals}
    if not depths_present:
        print(f"No finite soil moisture observations found in {obs_path}")
        return
    model_by_depth = {d: _interp_model_to_depth(soil_q, z_model, d)
                      for d in depths_present.values()}

    n = len(depths_present)
    fig, axes_ts, ax_profile, ncols, total_panels = _make_panel_grid(n)
    _plot_timeseries_panels(
        axes_ts, t_model, model_by_depth, t_obs_sl, obs_vals,
        depths_present, r"$\theta$ [m$^3$/m$^3$]", "Soil moisture",
        MOIS_MIN_SPAN, "bias={bias:+.3f}\nrmse={rmse:.3f}",
    )
    for i, ax in enumerate(axes_ts):
        _format_time_axis(ax, t_model[0], t_model[-1])
        if _is_bottom_row(i, total_panels, ncols):
            ax.set_xlabel("Time (UTC)")  # type: ignore[misc]

    depths_arr = np.array(list(depths_present.values()))
    obs_final = np.array([
        obs_vals[var][-1] if np.isfinite(obs_vals[var][-1])
        else np.nan
        for var in depths_present
    ])
    mask = np.isfinite(obs_final)
    _plot_final_profile(
        ax_profile, z_model, soil_q[-1], depths_arr[mask], obs_final[mask],
        r"$\theta$ [m$^3$/m$^3$]", "Final soil moisture profile",
        MOIS_MIN_SPAN,
    )

    fig.suptitle(  # type: ignore[misc]
        f"GABLS3 soil moisture: UtahLSM vs Cabauw\n"
        f"model: {model_path.name}   obs: {obs_path.name}",
        fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))

    if out_path is not None:
        fig.savefig(out_path, dpi=150)  # type: ignore[misc]
        print(f"Saved figure → {out_path}")
    if show:
        plt.show()  # type: ignore[misc]
    plt.close(fig)

    print("\nSoil moisture bias vs obs (m3/m3, overlap window):")
    for var, depth in depths_present.items():
        stats = _series_stats(t_model, model_by_depth[depth], t_obs_sl, obs_vals[var])
        if stats is not None:
            bias, rmse = stats
            print(f"  {var:>5} ({depth*100:4.0f} cm): "
                  f"bias={bias:+7.4f}  rmse={rmse:7.4f}")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for GABLS3 soil comparison script."""
    here = Path(__file__).resolve().parent
    default_model = here.parent / "lsm_gabls3_py.nc"
    default_obs_dir = (here.parent.parent / "cases" / "gabls3"
                      / "observations")
    default_out_dir = here.parent

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--model", type=Path, default=default_model,
                   help="UtahLSM NetCDF output (default: %(default)s)")
    p.add_argument("--obs-dir", type=Path, default=default_obs_dir,
                   help="Directory with CESAR soil observation files "
                        "(default: %(default)s)")
    p.add_argument("--out-dir", type=Path, default=default_out_dir,
                   help="Directory for output figures (default: %(default)s)")
    p.add_argument("--show", action="store_true",
                   help="Display figures interactively in addition to "
                        "saving them.")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    temp_obs = args.obs_dir / "cesar_soil_heat_lb1_t10_v1.0_200607.nc"
    mois_obs = args.obs_dir / "cesar_soil_water_lb1_t10_v1.1_200607.nc"

    if temp_obs.exists():
        compare_soil_temperature(
            args.model, temp_obs,
            args.out_dir / "gabls3_soil_temperature.png", args.show,
        )
    else:
        print(f"Skipping temperature: {temp_obs} not found")

    if mois_obs.exists():
        compare_soil_moisture(
            args.model, mois_obs,
            args.out_dir / "gabls3_soil_moisture.png", args.show,
        )
    else:
        print(f"Skipping moisture: {mois_obs} not found")
