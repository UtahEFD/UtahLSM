# UtahLSM
University of Utah Land Surface Model (UtahLSM), created at the University of Utah and the NOAA National Severe Storms Laboratory

## Overview
UtahLSM is a lightweight, fast-response land-surface model suitable for use in large-eddy simulation (LES) codes.

There is a C++ version of the model in the `/cpp` folder and a Python version of the model in the `/python` folder.

More documentation to come

## Compiling C++ Version

Compiling should be relatively quick and painless. To compile:

````shell
mkdir build
cd build
cmake -DNETCDF_DIR=/path/to/netcdf-c/include -DNETCDF_CXX_DIR=/path/to/netcdf-cxx/include ..
make
````
The executables will be placed in the `/bin` folder at the top level of `/cpp` file hierarchy.

## Documentation

### C++ Version

UtahLSM supports Doxygen documentation generation. After the `cmake` command above, optionally type:
````shell
make doc
````

The generated html and latex documentation will be placed in the `/doc` folder.

## Running

A series of offline test cases exist in the top-level `/cases` folder. Within each test case, run 
````shell
python makeInput.py
````
to generate the appropriate settings files and input data.

### C++

Running an offline test is easy. From the `/cpp` directory, run
````shell
./bin/utahlsm_offline -c [case name] -o [output filename, optional]
````
where `[case name]` is the folder name of a given case under `/cases`.

### Python

#### Installation & Setup

**Requirements:**
- Python 3.8 or later
- Dependencies: `numpy`, `netCDF4`

**Setup:**

A virtual environment is recommended to avoid conflicts with other projects:
```bash
cd python
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

Install dependencies:
```bash
pip install numpy netCDF4
```

#### Quick Start

The easiest way to get started is with the included GABLS3 test case:

```bash
cd python

# Run the included GABLS3 case (Cabauw, June 2006)
python utahlsm_offline.py -c GABLS3 -o my_output.nc
```

This will:
1. Load configuration from `../cases/GABLS3/lsm_namelist.json`
2. Read initial conditions from `../cases/GABLS3/lsm_init.nc`
3. Apply atmospheric forcing from `../cases/GABLS3/lsm_offline.nc`
4. Run a 9-hour simulation and save results to `my_output.nc`

Expected output: timing statistics printed to console, and a NetCDF file containing soil/surface variables.

#### Project Structure

```
python/
├── utahlsm_offline.py          # Main entry point for offline simulations
├── utahlsm/                    # Core package
│   ├── core.py                 # UtahLSM orchestrator class
│   ├── data_models.py          # Configuration and state dataclasses
│   ├── physics/                # Physics modules (soil, radiation, surface)
│   │   ├── soil/               # Soil models: Brooks-Corey, Campbell, Van Genuchten
│   │   ├── radiation/          # Radiation models
│   │   └── surface/            # Surface layer (MOST)
│   └── util/                   # Utilities
│       ├── solvers.py          # Numerical solvers (tridiagonal, root-finding)
│       ├── constants.py        # Physical constants
│       └── io/                 # Input/output
│           ├── input.py        # Reads namelist and NetCDF files
│           └── output.py       # Writes NetCDF output
├── data/                       # Soil property databases (JSON)
└── README.md
```

Key design patterns: **Strategy pattern** for physics modules (selected via configuration), **dataclasses** for type-safe state management, and **dependency injection** for modularity.

#### Configuration: The Namelist

The `lsm_namelist.json` file controls all simulation settings. Key sections:

**Numerics:**
- `iterations`: Max iterations for different solvers (recommend 20-50 for typical cases)
- `tolerances`: Convergence criteria (energy balance tolerance typically 0.1-1 W/m²)
- `heat_diffusion_back_weight`: Theta scheme parameter for the soil heat solver
  - `0.0` = Forward-Time Centered-Space (explicit, faster but less stable)
  - `0.5` = Crank-Nicolson (balanced)
  - `1.0` = Backward-Time Centered-Space (implicit, slower but very stable)
- `iterations.moisture_picard`: Iteration cap for the mixed-form soil moisture Picard solve
- `tolerances.moisture_picard`: Convergence tolerance for the mixed-form soil moisture Picard solve
- `tolerances.moisture_bounds`: Allowed post-solve overshoot before the model clips or raises

**Soil:**
- `model`: 1 (Brooks-Corey), 2 (Campbell), or 3 (Van Genuchten)
- `properties`: Public bundled soil database name ("rawls-brakensiek", "cosby", "clapp-hornberger", "carsel-parrish") or a custom JSON path. Bundled datasets include Peters-Lidard quartz fractions and Letts peat tiers internally. "carsel-parrish" provides native van Genuchten parameters (alpha/n) and is valid only with the Van Genuchten soil model.

**Surface:**
- `z_o`: Roughness length [m] (0.01-0.1 m typical for grass)
- `albedo`, `emissivity`: Surface properties [0-1]
- Model 1: Monin-Obukhov Similarity Theory (MOST)

**Radiation:**
- `model`: 0 (BasicRadiation—uses provided net radiation directly)
- `latitude`, `longitude`: For future solar calculations

**Time & Grid:**
- `utc_year`, `julian_day`, `utc_start`: Simulation timing
- `nx`, `ny`, `nz`: Grid dimensions (typically 1×1×N for offline column mode)

See `../cases/GABLS3/lsm_namelist.json` for a complete example.

#### Input & Output Files

**Input Files:**

1. **lsm_namelist.json** — JSON configuration file with all simulation settings
2. **lsm_init.nc** — NetCDF file with initial conditions:
   - `soil_z` [m]: Depth of soil layers
   - `soil_T` [K]: Initial soil temperature profile
   - `soil_q` [m³/m³]: Initial soil moisture profile
   - `soil_type` [string]: Soil type name (e.g., "clay", "sand")
3. **lsm_offline.nc** (optional) — NetCDF file with atmospheric forcing:
   - `tstep` [s]: Time step duration
   - `atm_U` [m/s]: Wind speed
   - `atm_T` [K]: Air temperature
   - `atm_q` [kg/kg]: Specific humidity
   - `atm_p` [Pa]: Atmospheric pressure
   - `R_net` [W/m²]: Net radiation

**Output File (lsm_*.nc):**
- Time series of soil temperature/moisture profiles, surface fluxes, and turbulent exchange
- Dimensions: `t` (time), `z` (soil depth)
- Analyze with `xarray`, `netCDF4`, or climate data tools

#### Workflow & Best Practices

**Complete simulation workflow:**

```
1. Prepare case data (or use GABLS3):
   cd ../cases/GABLS3
   python makeInput.py  # Generates input files from observations

2. Configure simulation (optional):
   Edit lsm_namelist.json to adjust:
   - Soil model and properties database
   - Surface parameters (roughness, albedo)
   - Solver tolerances and iteration limits
   - Output variables

3. Run simulation:
   python utahlsm_offline.py -c GABLS3 -o results.nc

4. Analyze output:
   import xarray as xr
   ds = xr.open_dataset('results.nc')
   # Plot soil temperature, surface fluxes, etc.
```

**Important conventions:**

- All physical quantities use **SI units** throughout: temperature [K], depth [m], fluxes [W/m²], etc.
- Soil layer indexing: 0 = top (surface), nz-1 = bottom
- Kinematic fluxes (in intermediate calculations) are [kg/(m²·s)] or [m/s]; these are converted to [W/m²] using air density and specific heat
- Stability parameter: Negative = unstable (heating), Positive = stable (cooling)

**Numerical stability considerations:**

- **Heat solver theta parameter**: Use `heat_diffusion_back_weight=1.0` (implicit) for maximum damping in the soil heat solve
- **Moisture Picard controls**: Increase `iterations.moisture_picard` or relax `tolerances.moisture_picard` if the Richards solve is not converging
- **SEB solver tolerance**: Tighter tolerance (e.g., 0.1 W/m²) improves surface energy balance closure but increases computation
- **Time step size**: Determined by forcing data; smaller steps improve accuracy but increase computation

**Troubleshooting:**

- **Heat solver oscillates**: Increase `heat_diffusion_back_weight` toward `1.0`; check soil temperature initial conditions for physical realism
- **Moisture solve fails bounds or convergence**: Increase `iterations.moisture_picard`, relax `tolerances.moisture_picard`, or inspect the forcing and soil hydraulic parameters
- **Slow convergence**: Increase `iterations` in numerics config; check atmospheric forcing data for realistic values
- **Memory issues with large grids**: Currently optimized for column simulations (nx=ny=1); 3D simulations require vectorization work
- **Unrealistic results**: Verify soil type matches your study domain; check that soil moisture is between residual and saturation

#### Extending the Model

UtahLSM uses a **strategy pattern** for physics modules, allowing easy addition of new schemes:

1. **Add a new soil model:**
   - Create class in `utahlsm/physics/soil/` inheriting from `Soil` abstract base
   - Implement required methods: `conductivity()`, `diffusivity()`, `water_potential()`, etc.
   - Register in `Soil.get_model()` factory with a new integer ID
   - Add ID to the namelist schema

2. **Add a new radiation model:**
   - Create class in `utahlsm/physics/radiation/` inheriting from `Radiation`
   - Implement `compute_net()` method
   - Register in `Radiation.get_model()` factory

#### Testing & Validation

The model validates input through:
- **JSON schema validation**: Namelist and soil property files are checked against defined schemas
- **Physical bounds checking**: Soil moisture, temperature, and other quantities are verified to be within realistic ranges
- **Case study validation**: Run against observations (GABLS3 includes observational data for comparison)

Example validation with GABLS3:
```python
import xarray as xr
import numpy as np

# Load model output and observations
model = xr.open_dataset('my_output.nc')
obs = xr.open_dataset('../cases/GABLS3/observations/obs_data.nc')

# Compare surface temperature
rmse = np.sqrt(((model.T_sfc - obs.T_obs)**2).mean())
print(f"Surface temperature RMSE: {rmse:.2f} K")
```
