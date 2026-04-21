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
      "coupling": 50
    },
    "tolerances": {
      "sfc_flux": 0.001,
      "seb_root": 1e-06,
      "smb_flux": 1e-06,
      "coupling_temp": 0.01,
      "coupling_mois": 1e-05
    },
    "diffusion_back_weight": 0.5,
    "warm_start_turbulence": true,
    "initialize_surface_temperature_from_seb": true
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
    "model": 1,
    "psi_stable": "beljaars-holtslag",
    "zeta_max": 0.5,
    "gustiness": 0.5,
    "gustiness_stable_only": true
  },
  "soil": {
    "properties": "rawls-brakensiek",
    "model": 2
  },
  "radiation": {
    "model": 0,
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
| `numerics` | Iteration caps, tolerances, and diffusion scheme controls |
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

`diffusion_back_weight`
: Theta-scheme backward weight for the soil diffusion solves. `0.5` gives a Crank-Nicolson style weighting.

`warm_start_turbulence`
: If `true`, initialize the turbulence state from the first forcing record before the first coupled solve.

`initialize_surface_temperature_from_seb`
: If `true`, perform a standalone surface energy balance initialization before the normal timestep loop.

`iterations`
: Integer caps for the nonlinear iterations used by surface fluxes, SEB bracketing, SEB root finding, SMB root finding, and the outer coupled solve.

`tolerances`
: Floating-point convergence criteria for the same iteration families.

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
: Number of soil layers.

### `surface`

`model`
: Currently only `1`, which selects `SurfaceMOST`.

`z_o`, `z_t`
: Momentum and thermal roughness lengths in meters.

`z_m`, `z_s`
: Reference heights for wind and scalar forcing in meters.

`albedo`, `emissivity`
: Surface radiative properties.

`psi_stable`
: Stable integrated MOST correction family. Supported values are `dyer-hicks`, `beljaars-holtslag`, and `cheng-brutsaert`.

`zeta_max`, `gustiness`, `gustiness_stable_only`
: Stability limiting and low-wind handling controls used by the surface flux solver.

### `soil`

`properties`
: Name of the bundled soil property table. The repository currently ships `cabauw-heinen`, `clapp-hornberger`, `cosby`, and `rawls-brakensiek`.

`model`
: Integer selector for the constitutive soil model:

- `1`: `BrooksCorey`
- `2`: `Campbell`
- `3`: `VanGenuchten`

### `radiation`

`model`
: `0` disables online radiation and expects net radiation from forcing. `1` selects `RadBasic`.

`latitude`, `longitude`
: Site coordinates in degrees.

### `canopy`

`model`
: `"none"` for bare soil or `"jarvis"` for the Jarvis resistance parameterization.

The remaining canopy fields are scalar-or-array parameters that can be supplied once or per horizontal column. Common fields are `lai`, `veg_fraction`, `rooting_depth`, `beta`, `rs_min`, `rs_max`, `rg_half`, `vpd_coef`, `t_opt`, and `t_coef`.

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

`time` is written automatically.

## Source of Truth

For the exact runtime types and field descriptions, see the docstring-generated entries for the configuration dataclasses in the [API Reference](reference.md).
