# UtahLSM Python Version

Utah Land-Surface Model (UtahLSM) - A fast-response land-surface model for simulating heat, moisture, and momentum exchange between land and atmosphere.

## Quick Start

### Installation for Development

If you're developing UtahLSM or want to couple it with other models, install it in editable mode:

```bash
# Navigate to the python directory
cd python

# Create a virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in editable mode with development dependencies
pip install -e ".[dev]"
```

This allows you to:
- Import `utahlsm` from anywhere
- Make changes to the code and see them immediately
- Run tests easily

### Installation for Using as a Package

If you just want to use UtahLSM (e.g., coupling with another model):

```bash
pip install -e .
```

Or once it's on PyPI:
```bash
pip install utahlsm
```

### Running the Offline Driver

Run a test case with the included offline driver:

```bash
# From the python directory
python utahlsm_offline.py -c GABLS3 -o output.nc
```

Options:
- `-c case_name`: Case name (must exist in `../cases/` directory)
- `-o output_file`: Output NetCDF filename

### Running Tests

```bash
# Run all tests
pytest tests/ -v

# Run specific test suite
pytest tests/test_solvers.py -v

# Run tests by marker
pytest -m solver              # Numerical solver tests
pytest -m datamodel          # Data model tests
pytest -m soil               # Soil physics tests

# Run with coverage
pytest tests/ --cov=utahlsm
```

For more information, see [TESTING.md](TESTING.md).

## Project Structure

```
python/
├── pyproject.toml              # Package configuration and metadata
├── pytest.ini                  # Pytest configuration
├── TESTING.md                  # Testing guide and reference
├── CLAUDE.md                   # Architecture and development guide
├── utahlsm_offline.py          # Offline driver (main entry point)
├── utahlsm/                    # Main package
│   ├── __init__.py             # Package exports
│   ├── core.py                 # UtahLSM orchestrator class
│   ├── data_models.py          # State and configuration dataclasses
│   ├── exceptions.py           # Custom exceptions
│   ├── physics/                # Physics modules
│   │   ├── radiation/          # Radiation models
│   │   ├── soil/               # Soil physics models
│   │   └── surface/            # Surface layer models (MOST)
│   └── util/                   # Utilities
│       ├── constants.py        # Physical constants
│       ├── solvers.py          # Numerical solvers
│       └── io/                 # Input/output handling
└── tests/                      # Test suite
    ├── conftest.py             # Pytest fixtures
    ├── test_solvers.py         # Numerical solver tests
    ├── test_data_models.py     # Data model tests
    ├── test_soil_models.py     # Soil physics tests
    └── ...
```

## Using UtahLSM with External Models

### Basic Usage

```python
from utahlsm import UtahLSM
from utahlsm.util.io import Input, Output
from utahlsm.data_models import AtmosphericState

# Load configuration and initialize
input_obj = Input('path/to/lsm_namelist.json')
output_obj = Output('output.nc', input_obj)

# Create LSM instance
lsm = UtahLSM(input_obj, output_obj)

# Time-stepping loop
for step, runtime in enumerate(forcing_times):
    # Update with atmospheric forcing
    atm_state = AtmosphericState(
        wind_speed=wind_speed,
        temperature=air_temp,
        specific_humidity=humidity,
        pressure=pressure,
        radiation_net=net_radiation
    )

    lsm.update(dt, runtime, atm_state)

    # Run physics solvers
    lsm.run(step, runtime)

    # Save results
    lsm.save(step, runtime)
```

### Coupling with an Atmospheric Model

Example of coupling UtahLSM with an external atmospheric model:

```python
from atmospheric_model import AtmosphericModel
from utahlsm import UtahLSM

# Initialize both models
atm_model = AtmosphericModel('config.yaml')
lsm = UtahLSM(input_obj, output_obj)

# Coupling loop
for step in range(num_steps):
    # Get atmospheric state from atm model
    atm_state = atm_model.get_atmospheric_state()

    # Update LSM with forcing
    lsm.update(dt, time, atm_state)
    lsm.run(step, time)

    # Get surface fluxes from LSM
    surface_fluxes = lsm.get_surface_fluxes()

    # Feed back to atmospheric model
    atm_model.set_surface_fluxes(surface_fluxes)
    atm_model.step(dt)
```

## Configuration

UtahLSM is configured through a JSON namelist file. See `../cases/gabls3/lsm_namelist.json` for an example.

Key sections:
- `time`: Simulation time parameters
- `grid`: Soil grid configuration
- `surface`: Surface properties (roughness, albedo, etc.)
- `soil`: Soil model selection and properties
- `radiation`: Radiation model selection
- `numerics`: Solver parameters and tolerances

## Dependencies

### Required
- **numpy** (≥1.20): Numerical computations
- **netCDF4** (≥1.5): NetCDF I/O

### Optional
- **pytest** (≥8.0): Testing
- **ruff**: Code linting
- **mypy**: Static type checking
- **matplotlib**, **xarray**: Data analysis and plotting

## Performance Notes

UtahLSM is designed to be lightweight and fast for:
- Large-Eddy Simulation (LES) applications
- Coupled atmosphere-land simulations
- High-frequency updates (small time steps)

Typical timing:
- Single time step: < 1 ms
- GABLS3 full run (10,000 steps): ~5-10 seconds

## Documentation

- **CLAUDE.md**: High-level architecture, design patterns, extending the model
- **TESTING.md**: Comprehensive testing guide
- **docstrings**: Inline documentation in source code
- **data_models.py**: Well-documented configuration and state structures

## Contributing

Guidelines for contributing:
1. Follow PEP 8 style guidelines
2. Use type hints
3. Add docstrings to all classes and methods
4. Write tests for new features
5. Ensure all tests pass before submitting changes

Run code quality checks:
```bash
ruff check .
mypy utahlsm
```

## License

MIT License - See LICENSE file for details

## Authors

- Jeremy A. Gibbs (University of Utah)
- Rob Stoll (University of Utah)
- Eric Pardyjak (University of Utah)
- Pete Willemsen (NCAR)

## Citation

If you use UtahLSM in your research, please cite:

```bibtex
@software{gibbs2024utahlsm,
  title={UtahLSM: Utah Land-Surface Model},
  author={Gibbs, Jeremy A. and Stoll, Rob and Pardyjak, Eric and Willemsen, Pete},
  year={2024},
  url={https://github.com/yourusername/UtahLSM}
}
```

## Support

For issues, questions, or suggestions:
- Create an issue on GitHub
- Contact: jeremy.gibbs@utah.edu

## Changelog

### Version 0.1.0 (Initial Release)
- Core LSM functionality
- Three soil models (Brooks-Corey, Campbell, Van Genuchten)
- Monin-Obukhov similarity theory surface layer
- Offline and coupled driver support
- Comprehensive test suite (63+ tests)
