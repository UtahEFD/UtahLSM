# Testing

UtahLSM uses `pytest` for the Python test suite. The commands below assume you are running from the `python/` directory.

## Install Test Tools

```bash
pip install -e ".[dev]"
```

## Run the Full Test Suite

```bash
pytest
```

## Useful Targeted Runs

```bash
pytest tests/test_solvers.py -v
pytest tests/test_data_models.py -v
pytest tests/test_soil_models.py -v
```

Marker-based selection is also available:

```bash
pytest -m solver
pytest -m soil
pytest -m canopy
pytest -m integration
```

## Documentation Smoke Test

After installing the docs extra, verify that the site builds cleanly:

```bash
pip install -e ".[docs]"
mkdocs build
```

This is the fastest way to catch broken docstring imports, invalid MkDocs configuration, or bad internal links.

## Offline Driver Regression Check

A useful integration sanity check is to run the bundled case:

```bash
python utahlsm_offline.py -c gabls3 -o lsm_gabls3.nc
```

That exercises namelist validation, NetCDF I/O, the timestep loop, and output writing in one path.
