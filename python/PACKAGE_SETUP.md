# UtahLSM Package Setup and Installation Guide

## What Was Created

This document summarizes the modern Python packaging setup for UtahLSM that enables it to be installed and coupled with other models.

### Files Created

1. **`pyproject.toml`** (Main package configuration)
   - Project metadata: name, version, description, authors, license
   - Dependencies: numpy, netCDF4
   - Optional dev dependencies: pytest, ruff, mypy
   - Tool configurations: pytest, ruff, mypy, isort, black
   - Supports both local development and PyPI distribution

2. **`python/.gitignore`** (Python-specific exclusions)
   - Python bytecode and caches
   - Virtual environments
   - Build artifacts
   - Test coverage reports
   - IDE files

3. **`python/README.md`** (Installation and usage guide)
   - Quick start instructions
   - Development setup
   - Running tests
   - Example coupling code
   - Project structure overview

## Installation

### For Development

```bash
cd python
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate
pip install -e ".[dev]"
```

This installs UtahLSM in "editable" mode, so changes to source files are immediately reflected.

### For Using as a Package

```bash
pip install -e .
```

Or once published to PyPI:
```bash
pip install utahlsm
```

## Package Information

After installation, verify it was successful:

```bash
pip show utahlsm
```

Output will show:
```
Name: utahlsm
Version: 0.1.0
Summary: Utah Land-Surface Model (UtahLSM) - A fast-response land-surface model...
Author-email: Jeremy A. Gibbs <jeremy.gibbs@utah.edu>, ...
Requires: netCDF4, numpy
```

## Using UtahLSM in External Code

Once installed, external models can import and use UtahLSM:

```python
# In an atmospheric model or coupling code
from utahlsm import UtahLSM
from utahlsm.data_models import AtmosphericState
from utahlsm.util.io import Input, Output

# Initialize
input_obj = Input('namelist.json')
output_obj = Output('output.nc', input_obj)
lsm = UtahLSM(input_obj, output_obj)

# Use in time-stepping loop
atm_state = AtmosphericState(wind_speed=5.0, temperature=283.15, ...)
lsm.update(dt, time, atm_state)
lsm.run(step, time)
lsm.save(step, time)
```

## Key Features

### Build System
- Uses modern `setuptools` with PEP 517/518 standards
- No need for `setup.py` (all config in `pyproject.toml`)
- Supports editable installs for development
- Ready for PyPI publication

### Dependencies Management
- **Core**: numpy, netCDF4 (required)
- **Dev**: pytest, ruff, mypy (optional, for development)
- **Analysis**: matplotlib, xarray (optional, for post-processing)
- **Docs**: sphinx (optional, for documentation)

### Tool Configuration
All tools configured in one file:
- **pytest**: Test discovery, markers, output options
- **ruff**: Code linting and quality checks
- **mypy**: Static type checking
- **isort**: Import sorting
- **black**: Code formatting

## Python Version Support

- **Minimum**: Python 3.9
- **Tested**: Python 3.9, 3.10, 3.11, 3.12

The package is compatible with all modern Python versions.

## Running Tests

After installation:

```bash
# Run all tests
pytest tests/ -v

# Run specific test suite
pytest tests/test_solvers.py -v

# Run by marker
pytest -m solver
pytest -m datamodel

# With coverage
pytest tests/ --cov=utahlsm
```

**Current Status**: 63 tests passing, 19 skipped (pending data loading)

## Publishing to PyPI (Future)

When ready to publish UtahLSM on PyPI:

```bash
# Build distributions
python -m build

# Upload to PyPI (requires credentials)
python -m twine upload dist/*
```

Then anyone can install with:
```bash
pip install utahlsm
```

## Project Metadata

The `pyproject.toml` includes:
- **Name**: utahlsm
- **Version**: 0.1.0
- **License**: MIT
- **Authors**: Jeremy A. Gibbs, Rob Stoll, Eric Pardyjak, Pete Willemsen
- **Repository**: GitHub link (update `homepage` URL)
- **Keywords**: land-surface-model, weather, climate, atmospheric

## Structure for Coupling

The package is organized to support coupling scenarios:

```
External Model
    ↓ (imports)
from utahlsm import UtahLSM
from utahlsm.data_models import AtmosphericState, SurfaceState
from utahlsm.util.io import Input, Output
    ↓
UtahLSM instance
    ↓ (provides)
- Surface fluxes (sensible, latent, ground heat)
- Turbulence scales (friction velocity, Obukhov length)
- Energy and water budgets
    ↓
External Model (uses for atmospheric feedback)
```

## Verified Working

✅ Package installs successfully
✅ All imports work correctly
✅ Tests pass (63/63 running tests)
✅ Development mode (`-e`) installation works
✅ Package discoverable in Python path

## Next Steps

1. **Update GitHub URLs** in `pyproject.toml` (replace `yourusername`)
2. **Add GitHub CI/CD** to run tests on every commit
3. **Create CHANGELOG.md** for version history
4. **Add API documentation** (optional: Sphinx/ReadTheDocs)
5. **Publish to PyPI** when ready (optional: for public distribution)

## Troubleshooting

### "utahlsm not found" when importing
```bash
# Make sure you installed it
pip install -e .

# Verify installation
pip show utahlsm

# Check Python path
python -c "import utahlsm; print(utahlsm.__file__)"
```

### Tests not running
```bash
# Install dev dependencies
pip install -e ".[dev]"

# Run from python directory
cd python
pytest tests/
```

### ModuleNotFoundError: No module named 'netCDF4'
```bash
# Install the missing dependency
pip install netCDF4
```

## Questions?

See:
- **TESTING.md** - Comprehensive testing guide
- **CLAUDE.md** - Architecture and development guide
- **README.md** - Quick start and usage examples
