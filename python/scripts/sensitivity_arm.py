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
"""Run ARM SGP sensitivity experiments and write comparison metrics.

The suite is intentionally structured as a screening design rather than a
large factorial search. Each named variant changes one physical block from
the current ARM namelist so the first pass can identify whether the dominant
error is hydraulic, thermal, canopy, ground-coupling, or surface-layer related.

Examples
--------
::

    python scripts/sensitivity_arm.py --dry-run --suite screen
    python scripts/sensitivity_arm.py --suite screen --rebuild-input
    python scripts/sensitivity_arm.py --suite full --groups soil,ground
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence, cast

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "utahlsm-matplotlib")
)

import netCDF4 as nc
import numpy as np

import compare_arm

WINDOWS: dict[str, tuple[np.datetime64 | None, np.datetime64 | None]] = {
    "all": (None, None),
    "evening": (np.datetime64("2017-06-17T18:00"), np.datetime64("2017-06-18T00:00")),
    "night": (np.datetime64("2017-06-18T00:00"), np.datetime64("2017-06-18T11:00")),
    "morning": (np.datetime64("2017-06-18T11:00"), np.datetime64("2017-06-18T15:00")),
    "afternoon": (np.datetime64("2017-06-18T15:00"), np.datetime64("2017-06-18T18:00")),
}

SOIL_EVENT_PERIODS: dict[str, tuple[np.datetime64, np.datetime64]] = {
    "storm_step": (
        np.datetime64("2017-06-18T06:00"),
        np.datetime64("2017-06-18T06:30"),
    ),
    "storm_total": (
        np.datetime64("2017-06-18T06:00"),
        np.datetime64("2017-06-18T12:00"),
    ),
    "day_total": (np.datetime64("2017-06-18T06:00"), np.datetime64("2017-06-18T18:00")),
}

SCALAR_DIAGNOSTICS = (
    "shf",
    "lhf",
    "ghf",
    "lhf_soil",
    "lhf_veg",
    "lhf_wet",
    "r_s",
    "theta_root",
    "ust",
    "obl",
    "seb_res",
    "seb_storage",
    "precip",
    "runoff",
    "canopy_water",
)


@dataclass(frozen=True)
class Variant:
    """One named ARM sensitivity experiment."""

    name: str
    group: str
    description: str
    patches: Mapping[str, Any]


VARIANTS: tuple[Variant, ...] = (
    Variant("baseline", "control", "current ARM namelist", {}),
    Variant(
        "soil_ch_campbell",
        "soil",
        "Clapp-Hornberger table with Campbell hydraulics",
        {"soil.model": "campbell"},
    ),
    Variant(
        "soil_ch_brooks_corey",
        "soil",
        "Clapp-Hornberger table with Brooks-Corey hydraulics",
        {"soil.model": "brooks-corey"},
    ),
    Variant(
        "soil_cosby_vg",
        "soil",
        "Cosby property table with van Genuchten hydraulics",
        {"soil.properties": "cosby", "soil.model": "van-genuchten"},
    ),
    Variant(
        "soil_rawls_vg",
        "soil",
        "Rawls-Brakensiek property table with van Genuchten hydraulics",
        {"soil.properties": "rawls-brakensiek", "soil.model": "van-genuchten"},
    ),
    Variant(
        "thermal_mccumber",
        "thermal",
        "McCumber-Pielke thermal conductivity",
        {"soil.thermal_conductivity_model": "mccumber-pielke"},
    ),
    Variant(
        "soil_ch_campbell_mccumber",
        "thermal",
        "Campbell hydraulics plus McCumber-Pielke thermal conductivity",
        {
            "soil.model": "campbell",
            "soil.thermal_conductivity_model": "mccumber-pielke",
        },
    ),
    Variant(
        "rground_000",
        "ground",
        "no in-canopy ground thermal resistance",
        {"canopy.r_ground": 0.0},
    ),
    Variant(
        "rground_100",
        "ground",
        "weaker in-canopy ground thermal resistance",
        {"canopy.r_ground": 100.0},
    ),
    Variant(
        "rground_400",
        "ground",
        "stronger in-canopy ground thermal resistance",
        {"canopy.r_ground": 400.0},
    ),
    Variant(
        "rground_800",
        "ground",
        "very strong in-canopy ground thermal resistance",
        {"canopy.r_ground": 800.0},
    ),
    Variant(
        "root_030",
        "canopy",
        "shallower root access",
        {"canopy.rooting_depth": 0.3},
    ),
    Variant(
        "root_100",
        "canopy",
        "deeper root access into clay loam",
        {"canopy.rooting_depth": 1.0},
    ),
    Variant(
        "root_150",
        "canopy",
        "deep root access through most of the active column",
        {"canopy.rooting_depth": 1.5},
    ),
    Variant(
        "beta_092",
        "canopy",
        "more top-heavy Jackson root profile",
        {"canopy.beta": 0.92},
    ),
    Variant(
        "beta_970",
        "canopy",
        "deeper Jackson root profile for grassland",
        {"canopy.beta": 0.97},
    ),
    Variant(
        "lai_150",
        "canopy",
        "sparser active canopy",
        {"canopy.lai": 1.5},
    ),
    Variant(
        "lai_350",
        "canopy",
        "denser active canopy",
        {"canopy.lai": 3.5},
    ),
    Variant(
        "veg_070",
        "canopy",
        "more exposed bare soil fraction",
        {"canopy.veg_fraction": 0.7},
    ),
    Variant(
        "veg_100",
        "canopy",
        "closed canopy",
        {"canopy.veg_fraction": 1.0},
    ),
    Variant(
        "rsmin_070",
        "canopy",
        "lower minimum stomatal resistance",
        {"canopy.rs_min": 70.0},
    ),
    Variant(
        "rsmin_200",
        "canopy",
        "higher minimum stomatal resistance",
        {"canopy.rs_min": 200.0},
    ),
    Variant(
        "vpdcoef_050",
        "canopy",
        "weaker VPD stress",
        {"canopy.vpd_coef": 5.0e-5},
    ),
    Variant(
        "vpdcoef_200",
        "canopy",
        "stronger VPD stress",
        {"canopy.vpd_coef": 2.0e-4},
    ),
    Variant(
        "rough_z0_030",
        "surface",
        "lower aerodynamic roughness",
        {"surface.z_o": 0.03},
    ),
    Variant(
        "rough_z0_100",
        "surface",
        "higher aerodynamic roughness",
        {"surface.z_o": 0.10},
    ),
    Variant(
        "rough_zt_00005",
        "surface",
        "smaller scalar roughness",
        {"surface.z_t": 0.00005},
    ),
    Variant(
        "rough_zt_005",
        "surface",
        "larger scalar roughness",
        {"surface.z_t": 0.005},
    ),
    Variant(
        "gust_000",
        "surface",
        "remove stable gustiness floor",
        {"surface.gustiness": 0.0},
    ),
    Variant(
        "gust_200",
        "surface",
        "stronger stable gustiness floor",
        {"surface.gustiness": 2.0},
    ),
    Variant(
        "zeta_500",
        "surface",
        "looser stable MOST zeta clamp",
        {"surface.zeta_max": 5.0},
    ),
    Variant(
        "theta_implicit",
        "numerics",
        "fully implicit heat diffusion",
        {"numerics.heat_diffusion_back_weight": 1.0},
    ),
    Variant(
        "combo_campbell_rground400",
        "combo",
        "Campbell hydraulics plus stronger ground resistance",
        {"soil.model": "campbell", "canopy.r_ground": 400.0},
    ),
    Variant(
        "combo_campbell_rground800",
        "combo",
        "Campbell hydraulics plus very strong ground resistance",
        {"soil.model": "campbell", "canopy.r_ground": 800.0},
    ),
    Variant(
        "combo_rsmin070_rground400",
        "combo",
        "lower minimum stomatal resistance plus stronger ground resistance",
        {"canopy.rs_min": 70.0, "canopy.r_ground": 400.0},
    ),
    Variant(
        "combo_lai350_rground400",
        "combo",
        "denser canopy plus stronger ground resistance",
        {"canopy.lai": 3.5, "canopy.r_ground": 400.0},
    ),
    Variant(
        "combo_campbell_rsmin070",
        "combo",
        "Campbell hydraulics plus lower minimum stomatal resistance",
        {"soil.model": "campbell", "canopy.rs_min": 70.0},
    ),
    Variant(
        "combo_campbell_rsmin070_rground400",
        "combo",
        "Campbell hydraulics, lower rs_min, and stronger ground resistance",
        {"soil.model": "campbell", "canopy.rs_min": 70.0, "canopy.r_ground": 400.0},
    ),
)

SUITES: dict[str, tuple[str, ...]] = {
    "smoke": ("baseline",),
    "screen": (
        "baseline",
        "soil_ch_campbell",
        "soil_cosby_vg",
        "soil_rawls_vg",
        "thermal_mccumber",
        "soil_ch_campbell_mccumber",
        "rground_000",
        "rground_400",
        "rground_800",
        "root_030",
        "root_100",
        "beta_970",
        "lai_150",
        "lai_350",
        "veg_070",
        "veg_100",
        "rsmin_070",
        "rsmin_200",
        "vpdcoef_050",
        "rough_z0_030",
        "rough_z0_100",
        "rough_zt_005",
        "gust_000",
        "gust_200",
    ),
    "interaction": (
        "baseline",
        "rground_400",
        "rground_800",
        "soil_ch_campbell",
        "rsmin_070",
        "lai_350",
        "combo_campbell_rground400",
        "combo_campbell_rground800",
        "combo_rsmin070_rground400",
        "combo_lai350_rground400",
        "combo_campbell_rsmin070",
        "combo_campbell_rsmin070_rground400",
    ),
    "full": tuple(variant.name for variant in VARIANTS),
}


def _repo_root() -> Path:
    """Return the repository root based on this script location."""
    return Path(__file__).resolve().parents[2]


def _variant_by_name() -> dict[str, Variant]:
    return {variant.name: variant for variant in VARIANTS}


def _select_variants(
    suite: str,
    names: Sequence[str],
    groups: Sequence[str],
    limit: int | None,
) -> list[Variant]:
    variants = _variant_by_name()
    requested = list(names) if names else list(SUITES[suite])
    unknown = sorted(set(requested) - set(variants))
    if unknown:
        raise ValueError("Unknown ARM sensitivity variant(s): " + ", ".join(unknown))

    selected = [variants[name] for name in requested]
    if groups:
        wanted_groups = set(groups)
        selected = [
            variant
            for variant in selected
            if variant.group in wanted_groups or variant.name == "baseline"
        ]
    if limit is not None:
        selected = selected[:limit]
    return selected


def _split_csv_arg(raw: str | None) -> tuple[str, ...]:
    if raw is None or raw.strip() == "":
        return ()
    return tuple(item.strip() for item in raw.split(",") if item.strip())


def _set_path(config: dict[str, Any], path: str, value: Any) -> None:
    current: dict[str, Any] = config
    keys = path.split(".")
    for key in keys[:-1]:
        current = cast(dict[str, Any], current[key])
    current[keys[-1]] = value


def _apply_variant(base: Mapping[str, Any], variant: Variant) -> dict[str, Any]:
    config = json.loads(json.dumps(base))
    for path, value in variant.patches.items():
        _set_path(config, path, value)
    return cast(dict[str, Any], config)


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=4)
        f.write("\n")


def _run_command(
    command: Sequence[str],
    cwd: Path,
    log_path: Path,
) -> None:
    result = subprocess.run(
        command,
        cwd=cwd,
        check=False,
        text=True,
        capture_output=True,
    )
    log_path.write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit code {result.returncode}: "
            + " ".join(command)
            + f"\nSee log: {log_path}"
        )


def _window_mask(
    times: np.ndarray,
    window: tuple[np.datetime64 | None, np.datetime64 | None],
) -> np.ndarray:
    start, end = window
    mask = np.ones(times.shape, dtype=bool)
    if start is not None:
        mask &= times >= start
    if end is not None:
        mask &= times <= end
    return mask


def _interp_model(
    t_model: np.ndarray,
    values: np.ndarray,
    t_obs: np.ndarray,
) -> np.ndarray:
    return np.asarray(
        np.interp(
            (t_obs - t_obs[0]).astype("timedelta64[s]").astype(float),
            (t_model - t_obs[0]).astype("timedelta64[s]").astype(float),
            values,
        )
    )


def _series_metrics(
    t_model: np.ndarray,
    model_values: np.ndarray | None,
    t_obs: np.ndarray,
    obs_values: np.ndarray,
    window: tuple[np.datetime64 | None, np.datetime64 | None],
) -> dict[str, float] | None:
    if model_values is None:
        return None
    mask = _window_mask(t_obs, window) & np.isfinite(obs_values)
    if not mask.any():
        return None
    t_obs_w = t_obs[mask]
    obs_w = obs_values[mask]
    model_w = _interp_model(t_model, model_values, t_obs_w)
    good = np.isfinite(model_w) & np.isfinite(obs_w)
    if not good.any():
        return None
    diff = model_w[good] - obs_w[good]
    return {
        "bias": float(np.mean(diff)),
        "rmse": float(np.sqrt(np.mean(diff**2))),
        "obs_mean": float(np.mean(obs_w[good])),
        "model_mean": float(np.mean(model_w[good])),
        "n": float(np.count_nonzero(good)),
    }


def _nearest_time_index(times: np.ndarray, target: np.datetime64) -> int:
    offsets = np.abs((times - target).astype("timedelta64[s]").astype(float))
    return int(np.argmin(offsets))


def _load_flux_obs(
    model: compare_arm.ModelData,
    obs_dir: Path,
) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray | None]]:
    with (
        nc.Dataset(obs_dir / "arm_surf_flux.nc") as fds,
        nc.Dataset(obs_dir / "arm_surf_radn.nc") as rds,
    ):
        t_flux = compare_arm._obs_times(fds)
        t_radn = compare_arm._obs_times(rds)
        flux_slice = compare_arm._overlap_slice(model.time, t_flux)
        radn_slice = compare_arm._overlap_slice(model.time, t_radn)
        t_flux_sl = t_flux[flux_slice]
        t_radn_sl = t_radn[radn_slice]
        shf_obs = compare_arm._clean_variable(
            fds, "corrected_sensible_heat_flux", flux_slice
        )
        lhf_obs = compare_arm._clean_variable(
            fds, "corrected_latent_heat_flux", flux_slice
        )
        ghf_obs = -compare_arm._clean_variable(
            rds, "surface_soil_heat_flux_avg", radn_slice
        )
    return {
        "H": (t_flux_sl, shf_obs, model.shf),
        "LE": (t_flux_sl, lhf_obs, model.lhf),
        "G0": (t_radn_sl, ghf_obs, model.ghf),
    }


def _flux_rows(
    variant: Variant,
    model: compare_arm.ModelData,
    obs_dir: Path,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    obs_by_quantity = _load_flux_obs(model, obs_dir)
    for window_name, window in WINDOWS.items():
        for quantity, (t_obs, obs_values, model_values) in obs_by_quantity.items():
            metrics = _series_metrics(
                model.time, model_values, t_obs, obs_values, window
            )
            if metrics is None:
                continue
            rows.append(
                {
                    "variant": variant.name,
                    "group": variant.group,
                    "window": window_name,
                    "quantity": quantity,
                    **metrics,
                }
            )
    return rows


def _soil_rows(
    variant: Variant,
    model: compare_arm.ModelData,
    obs_dir: Path,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with nc.Dataset(obs_dir / "arm_soil_data.nc") as ods:
        t_obs = compare_arm._obs_times(ods)
        obs_slice = compare_arm._overlap_slice(model.time, t_obs)
        t_obs_sl = t_obs[obs_slice]

        _, temp_obs, _ = compare_arm._load_arm_soil_series(
            ods,
            compare_arm.SOIL_TEMP_VARS,
            obs_slice,
            scale=1.0,
            offset=273.15,
        )
        _, mois_obs, _ = compare_arm._load_arm_soil_series(
            ods,
            compare_arm.SOIL_MOIS_VARS,
            obs_slice,
            scale=0.01,
            offset=0.0,
        )

    profiles = {
        "soil_temperature": model.soil_t,
        "soil_moisture": model.soil_q,
    }
    obs_sets = {
        "soil_temperature": temp_obs,
        "soil_moisture": mois_obs,
    }
    for quantity, profile in profiles.items():
        if profile is None:
            continue
        for label, depth in compare_arm.SOIL_DEPTHS.items():
            if label not in obs_sets[quantity]:
                continue
            model_at_depth = compare_arm._interp_model_to_depth(
                profile, model.soil_z, depth
            )
            metrics = _series_metrics(
                model.time,
                model_at_depth,
                t_obs_sl,
                obs_sets[quantity][label],
                WINDOWS["all"],
            )
            if metrics is None:
                continue
            rows.append(
                {
                    "variant": variant.name,
                    "group": variant.group,
                    "quantity": quantity,
                    "depth_label": label,
                    "depth_m": depth,
                    **metrics,
                }
            )
    return rows


def _soil_event_rows(
    variant: Variant,
    model: compare_arm.ModelData,
    obs_dir: Path,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with nc.Dataset(obs_dir / "arm_soil_data.nc") as ods:
        t_obs = compare_arm._obs_times(ods)
        obs_slice = compare_arm._overlap_slice(model.time, t_obs)
        t_obs_sl = t_obs[obs_slice]

        _, temp_obs, _ = compare_arm._load_arm_soil_series(
            ods,
            compare_arm.SOIL_TEMP_VARS,
            obs_slice,
            scale=1.0,
            offset=273.15,
        )
        _, mois_obs, _ = compare_arm._load_arm_soil_series(
            ods,
            compare_arm.SOIL_MOIS_VARS,
            obs_slice,
            scale=0.01,
            offset=0.0,
        )

    profiles = {
        "soil_temperature": model.soil_t,
        "soil_moisture": model.soil_q,
    }
    obs_sets = {
        "soil_temperature": temp_obs,
        "soil_moisture": mois_obs,
    }

    for quantity, profile in profiles.items():
        if profile is None:
            continue
        for label, depth in compare_arm.SOIL_DEPTHS.items():
            if label not in obs_sets[quantity]:
                continue
            model_at_depth = compare_arm._interp_model_to_depth(
                profile, model.soil_z, depth
            )
            obs_at_depth = obs_sets[quantity][label]
            for period, (start, end) in SOIL_EVENT_PERIODS.items():
                obs_i0 = _nearest_time_index(t_obs_sl, start)
                obs_i1 = _nearest_time_index(t_obs_sl, end)
                model_i0 = _nearest_time_index(model.time, start)
                model_i1 = _nearest_time_index(model.time, end)

                obs_before = float(obs_at_depth[obs_i0])
                obs_after = float(obs_at_depth[obs_i1])
                model_before = float(model_at_depth[model_i0])
                model_after = float(model_at_depth[model_i1])
                if not np.all(
                    np.isfinite([obs_before, obs_after, model_before, model_after])
                ):
                    continue

                obs_delta = obs_after - obs_before
                model_delta = model_after - model_before
                rows.append(
                    {
                        "variant": variant.name,
                        "group": variant.group,
                        "period": period,
                        "quantity": quantity,
                        "depth_label": label,
                        "depth_m": depth,
                        "obs_before": obs_before,
                        "obs_after": obs_after,
                        "obs_delta": obs_delta,
                        "model_before": model_before,
                        "model_after": model_after,
                        "model_delta": model_delta,
                        "delta_error": model_delta - obs_delta,
                    }
                )
    return rows


def _diagnostic_rows(
    variant: Variant,
    model: compare_arm.ModelData,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with nc.Dataset(model.path) as ds:
        for name in SCALAR_DIAGNOSTICS:
            values = compare_arm._model_var(ds, name)
            if values is None:
                continue
            for window_name, window in WINDOWS.items():
                mask = _window_mask(model.time, window) & np.isfinite(values)
                if not mask.any():
                    continue
                values_w = values[mask]
                row: dict[str, Any] = {
                    "variant": variant.name,
                    "group": variant.group,
                    "window": window_name,
                    "quantity": name,
                    "mean": float(np.mean(values_w)),
                    "min": float(np.min(values_w)),
                    "max": float(np.max(values_w)),
                }
                if name in {"precip", "runoff"} and model.time.size > 1:
                    dt = float(
                        np.median(
                            np.diff(model.time).astype("timedelta64[ms]").astype(float)
                        )
                        / 1000.0
                    )
                    row["total"] = float(np.sum(values_w) * dt)
                rows.append(row)
    return rows


def _score_rows(
    variants: Iterable[Variant],
    flux_rows: Sequence[Mapping[str, Any]],
    soil_rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for variant in variants:
        all_flux = [
            float(row["rmse"])
            for row in flux_rows
            if row["variant"] == variant.name and row["window"] == "all"
        ]
        night_flux = [
            float(row["rmse"])
            for row in flux_rows
            if row["variant"] == variant.name
            and row["window"] in {"evening", "night", "morning"}
        ]
        temp = [
            float(row["rmse"])
            for row in soil_rows
            if row["variant"] == variant.name and row["quantity"] == "soil_temperature"
        ]
        mois = [
            float(row["rmse"])
            for row in soil_rows
            if row["variant"] == variant.name and row["quantity"] == "soil_moisture"
        ]
        rows.append(
            {
                "variant": variant.name,
                "group": variant.group,
                "flux_rmse_score": _rms(all_flux),
                "transition_flux_rmse_score": _rms(night_flux),
                "soil_temp_rmse_mean": float(np.mean(temp)) if temp else np.nan,
                "soil_mois_rmse_mean": float(np.mean(mois)) if mois else np.nan,
            }
        )
    return rows


def _rms(values: Sequence[float]) -> float:
    if not values:
        return np.nan
    arr = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(arr**2)))


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _print_plan(variants: Sequence[Variant]) -> None:
    print("ARM sensitivity variants:")
    for variant in variants:
        patches = ", ".join(f"{key}={value}" for key, value in variant.patches.items())
        patch_text = patches if patches else "no changes"
        print(f"  {variant.name:28s} [{variant.group:8s}] {patch_text}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--suite",
        choices=sorted(SUITES),
        default="screen",
        help="Named variant suite to run. Default: %(default)s",
    )
    parser.add_argument(
        "--only",
        default=None,
        help="Comma-separated variant names. Overrides --suite selection.",
    )
    parser.add_argument(
        "--groups",
        default=None,
        help="Comma-separated variant groups to keep. Baseline is always included.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Run only the first N selected variants.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory. Default: REPO/python/arm_sensitivity",
    )
    parser.add_argument(
        "--obs-dir",
        type=Path,
        default=None,
        help="ARM observation directory. Default: REPO/cases/arm/observations",
    )
    parser.add_argument(
        "--rebuild-input",
        action="store_true",
        help="Run cases/arm/make_input.py before the suite.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Reuse existing NetCDF outputs when present.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print selected variants without running UtahLSM.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    repo = _repo_root()
    python_dir = repo / "python"
    arm_dir = repo / "cases" / "arm"
    obs_dir = args.obs_dir or arm_dir / "observations"
    out_dir = args.out_dir or python_dir / "arm_sensitivity"
    log_dir = out_dir / "logs"
    namelist_path = arm_dir / "lsm_namelist.json"

    selected = _select_variants(
        suite=args.suite,
        names=_split_csv_arg(args.only),
        groups=_split_csv_arg(args.groups),
        limit=args.limit,
    )

    _print_plan(selected)
    if args.dry_run:
        return

    out_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    if args.rebuild_input:
        _run_command(
            [sys.executable, "make_input.py"],
            cwd=arm_dir,
            log_path=log_dir / "make_input.log",
        )

    if not namelist_path.exists():
        raise FileNotFoundError(
            f"{namelist_path} does not exist. Run cases/arm/make_input.py "
            "or pass --rebuild-input."
        )

    base_config = cast(
        dict[str, Any],
        json.loads(namelist_path.read_text(encoding="utf-8")),
    )
    original_text = namelist_path.read_text(encoding="utf-8")

    all_flux_rows: list[dict[str, Any]] = []
    all_soil_rows: list[dict[str, Any]] = []
    all_soil_event_rows: list[dict[str, Any]] = []
    all_diag_rows: list[dict[str, Any]] = []

    try:
        for idx, variant in enumerate(selected, start=1):
            print(f"\n[{idx}/{len(selected)}] {variant.name}: {variant.description}")
            output_nc = out_dir / f"{variant.name}.nc"
            if args.resume and output_nc.exists():
                print(f"  using existing {output_nc}")
            else:
                config = _apply_variant(base_config, variant)
                _write_json(namelist_path, config)
                _run_command(
                    [
                        sys.executable,
                        "utahlsm_offline.py",
                        "-c",
                        "arm",
                        "-o",
                        str(output_nc),
                    ],
                    cwd=python_dir,
                    log_path=log_dir / f"{variant.name}.log",
                )

            model = compare_arm.load_model(output_nc)
            all_flux_rows.extend(_flux_rows(variant, model, obs_dir))
            all_soil_rows.extend(_soil_rows(variant, model, obs_dir))
            all_soil_event_rows.extend(_soil_event_rows(variant, model, obs_dir))
            all_diag_rows.extend(_diagnostic_rows(variant, model))
    finally:
        namelist_path.write_text(original_text, encoding="utf-8")

    score_rows = _score_rows(selected, all_flux_rows, all_soil_rows)
    _write_csv(out_dir / "flux_metrics.csv", all_flux_rows)
    _write_csv(out_dir / "soil_metrics.csv", all_soil_rows)
    _write_csv(out_dir / "soil_event_metrics.csv", all_soil_event_rows)
    _write_csv(out_dir / "diagnostics.csv", all_diag_rows)
    _write_csv(out_dir / "scores.csv", score_rows)

    ranked = sorted(
        score_rows,
        key=lambda row: (
            float(row["transition_flux_rmse_score"])
            if np.isfinite(float(row["transition_flux_rmse_score"]))
            else np.inf
        ),
    )
    print("\nTop variants by transition flux RMSE score:")
    for row in ranked[:8]:
        print(
            f"  {row['variant']:28s} "
            f"transition={float(row['transition_flux_rmse_score']):6.2f} "
            f"all_flux={float(row['flux_rmse_score']):6.2f} "
            f"soilT={float(row['soil_temp_rmse_mean']):5.2f}"
        )
    print(f"\nWrote sensitivity outputs -> {out_dir}")


if __name__ == "__main__":
    main()
