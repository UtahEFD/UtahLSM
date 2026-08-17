# Namelist Reference

UtahLSM uses a schema-validated JSON namelist. The schema lives in `utahlsm/util/io/schema_namelist.json`, and the configuration dataclasses in `utahlsm.data_models` are the runtime source of truth.

## Complete Example

```json
{
  "general": {
    "log_level": "info"
  },
  "numerics": {
    "iterations": {
      "sfc_flux": 100,
      "seb_bracket": 100,
      "seb_root": 100,
      "smb_flux": 100,
      "moisture_picard": 25,
      "coupling": 50
    },
    "tolerances": {
      "sfc_flux": 0.001,
      "seb_root": 1e-06,
      "smb_flux": 1e-06,
      "moisture_picard": 1e-08,
      "moisture_bounds": 1e-06,
      "coupling_temp": 0.01,
      "coupling_mois": 1e-05
    },
    "heat_diffusion_back_weight": 0.5
  },
  "time": {
    "utc_start": 0,
    "utc_year": 2006,
    "julian_day": 183
  },
  "grid": {
    "nx": 1,
    "ny": 1,
    "nz": 11
  },
  "surface": {
    "z_o": 0.03,
    "z_t": 0.0046,
    "z_m": 10.0,
    "z_s": 2.0,
    "albedo": 0.33,
    "emissivity": 0.99,
    "model": "most",
    "psi_stable": "beljaars-holtslag",
    "zeta_max": 0.5,
    "gustiness": 0.5,
    "gustiness_stable_only": true
  },
  "soil": {
    "properties": "rawls-brakensiek",
    "model": "campbell"
  },
  "radiation": {
    "model": "forcing",
    "latitude": 51.9711,
    "longitude": -4.9267
  },
  "canopy": {
    "model": "jarvis",
    "lai": 3.0,
    "veg_fraction": 0.95,
    "rooting_depth": 0.4,
    "beta": 0.965,
    "rs_min": 40.0,
    "rs_max": 5000.0,
    "rg_half": 100.0,
    "vpd_coef": 1.0e-4,
    "t_opt": 298.0,
    "t_coef": 1.6e-3
  },
  "output": {
    "save": true,
    "fields": ["all"]
  }
}
```

## Section Summary

| Section | Purpose |
| --- | --- |
| `general` | Logging verbosity |
| `numerics` | Iteration caps, tolerances, and heat-solver controls |
| `time` | UTC and Julian-day context for radiation timing |
| `grid` | Horizontal column count and number of soil levels |
| `surface` | Roughness lengths, measurement heights, and MOST options |
| `soil` | Soil dataset and constitutive model selection |
| `radiation` | Radiation model selection plus site latitude and longitude |
| `canopy` | Optional vegetation parameters |
| `output` | NetCDF output toggle and field list |

## Section Details

### `general`

`log_level`
: `info` or `debug`.

### `numerics`

`heat_diffusion_back_weight`
: Theta-scheme backward weight for the soil heat diffusion solve. `0.5` gives a Crank-Nicolson style weighting.

`iterations`
: Integer caps for the nonlinear iterations used by surface fluxes, SEB bracketing, SEB root finding, SMB root finding, the mixed-form soil moisture Picard solve, and the outer coupled solve.

`tolerances`
: Strictly positive floating-point convergence criteria for the same iteration families, plus the admissible post-solve soil moisture bounds overshoot before clipping or failure.

### `time`

`utc_start`
: Seconds from midnight UTC at the start of the simulation.

`utc_year`
: Calendar year used for leap-year handling and solar calculations.

`julian_day`
: Day of year used by the radiation model.

### `grid`

`nx`, `ny`
: Horizontal column counts. The model stores runtime state as flattened `ncol = nx × ny`.

`nz`
: Number of soil layers. At least three layers are required by the diffusion solvers.

### `surface`

`model`
: Currently only `"most"`, which selects `SurfaceMOST`.

`z_o`, `z_t`
: Positive momentum and thermal roughness lengths in meters.

`z_m`, `z_s`
: Positive reference heights for wind and scalar forcing in meters. The wind height must exceed `z_o`, and the scalar height must exceed `z_t`.

`albedo`, `emissivity`
: Surface radiative properties.

`psi_stable`
: Stable integrated MOST correction family. Supported values are `dyer-hicks`, `beljaars-holtslag`, and `cheng-brutsaert`.

`zeta_max`, `gustiness`, `gustiness_stable_only`
: Stability limiting and low-wind handling controls used by the surface flux solver.

### `soil`

`properties`
: Name of the bundled soil property table, or a path to a custom JSON file. Public bundled datasets are `clapp-hornberger`, `cosby`, `rawls-brakensiek`, and `carsel-parrish`. These are the base hydraulic-property tables exposed in the namelist; internally, each bundled dataset is supplemented with Peters-Lidard et al. (1998) texture-class quartz fractions for Johansen thermal conductivity and Letts et al. (2000) peat tiers (`peat_fibric`, `peat_hemic`, `peat_sapric`) for organic layers.

When `thermal_conductivity_model` is `"johansen"`, custom mineral soil types must provide `quartz_fraction`. Missing quartz is allowed only for named `peat_*` organic layers; otherwise the model raises a namelist error instead of guessing.

Soil-type entries define the retention curve either with Campbell-style `b`/`psi_sat` or with native van Genuchten `alpha` [1/m] and `n` — exactly one of the two pairs. Native van Genuchten entries (including all of `carsel-parrish`, from Carsel & Parrish 1988) are valid only with the `van-genuchten` soil model; `campbell` and `brooks-corey` reject them because no measured `b` exponent exists.

A dataset JSON may declare `"includes": ["a", "b", ...]` to merge soil types from other public bundled datasets or file paths. The current dataset's own `soil_types` are merged last and take precedence. Conflicts are governed by `"on_conflict"`: `"error"` (default), `"prefer_first"`, or `"prefer_last"`.

`model`
: Constitutive soil model selector:

- `"brooks-corey"`: Brooks-Corey hydraulics
- `"campbell"`: Campbell hydraulics
- `"van-genuchten"`: van Genuchten hydraulics

When `macropore_fraction` is greater than zero, both `macropore_z_top` and `macropore_z_bottom` are required. The bottom must be deeper than the top, remain inside the modeled soil domain, and enclose at least one prognostic soil-layer node.

### `radiation`

`model`
: `"forcing"` passes through radiation components from the forcing file. `"basic"` computes clear-sky radiation from site coordinates.

`latitude`, `longitude`
: Site coordinates in degrees.

### `canopy`

`model`
: `"none"` for bare soil or `"jarvis"` for the Jarvis resistance parameterization.

The remaining canopy fields are scalar-or-array parameters that can be supplied once or per horizontal column. Common fields are `lai`, `veg_fraction`, `rooting_depth`, `beta`, `rs_min`, `rs_max`, `rg_half`, `vpd_coef`, `t_opt`, `t_coef`, `r_ground`, `water_capacity_lai`, and `wet_cooling_max`.

`water_capacity_lai`
: Wet-canopy water holding capacity per LAI [kg m-2 per LAI]. The model uses `veg_fraction * lai * water_capacity_lai` as the per-ground-area dew/interception storage capacity.

`wet_cooling_max`
: Maximum diagnostic nighttime wet-canopy cooling below the soil/radiative skin [K] used by the dewfall path.

### `output`

`save`
: Enables or disables NetCDF writing.

`fields`
: Requested output variables. Use `["all"]` to write every available field for the active configuration, or provide an explicit list such as:

- `soil_z`
- `soil_type`
- `soil_T`
- `soil_q`
- `ust`
- `obl`
- `shf`
- `lhf`
- `ghf`
- `r_s`
- `theta_root`
- `lhf_soil`
- `lhf_veg`
- `lhf_wet`
- `canopy_water`

`time` is written automatically.

## Source of Truth

For the exact runtime types and field descriptions, see the docstring-generated entries for the configuration dataclasses in the [API Reference](reference.md).
