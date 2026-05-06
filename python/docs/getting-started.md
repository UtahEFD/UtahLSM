# Getting Started

This guide covers local installation, running the bundled offline case, and building the documentation site.

## Prerequisites

UtahLSM currently targets:

- Python 3.9 or newer
- NumPy and netCDF4 runtime dependencies
- A checkout of the full `UtahLSM` repository, since the sample case files live in `../cases/`

## Installation

From the `python/` directory:

```bash
pip install -e .
```

For development tools:

```bash
pip install -e ".[dev]"
```

For documentation tools:

```bash
pip install -e ".[docs]"
```

## Run the Bundled GABLS3 Case

The repository includes an offline reference case in `../cases/gabls3/`. The offline driver expects a case directory under `../cases/<case>/` containing:

- `lsm_namelist.json`
- `lsm_init.nc`
- `lsm_offline.nc`

`lsm_offline.nc` may include an optional `seb_storage` forcing term
[W m-2]. When present, the surface energy budget is solved as
`Rn - H - LE - G - seb_storage = 0`; when omitted, this term defaults to zero.

```bash
python utahlsm_offline.py -c gabls3 -o lsm_gabls3.nc
```

This driver will:

1. Read `../cases/gabls3/lsm_namelist.json`
2. Load `lsm_init.nc` and `lsm_offline.nc`
3. Construct `Input`, `Output`, and `UtahLSM`
4. March through the forcing time series and write NetCDF output

If `-o` is omitted, the default output name is `lsm_gabls3_py.nc`.

## Use UtahLSM as a Library

The core package API is centered around `Input`, `Output`, and `UtahLSM`:

```python
from utahlsm import Input, Output, UtahLSM

input_lsm = Input(
    "../cases/gabls3/lsm_namelist.json",
    "../cases/gabls3/lsm_init.nc",
    "../cases/gabls3/lsm_offline.nc",
)
output_lsm = Output("lsm_gabls3.nc", enabled=input_lsm.output.save)
lsm = UtahLSM(input_lsm, output_lsm)
```

For an offline run, iterate over `input_lsm.forcing.atmos`, call `lsm.update(...)`, `lsm.run()`, and `lsm.save(...)`.

## Build the Documentation Site

The documentation configuration lives in `python/mkdocs.yml`, mirroring the structure already used in `pyburgers`.

```bash
mkdocs serve
```

That starts a local preview server. For a build artifact suitable for deployment:

```bash
mkdocs build
```

## Where to Look Next

- [Architecture](architecture.md) for the model layout
- [Namelist Reference](namelist.md) for configuration details
- [API Reference](reference.md) for docstring-generated class and function docs
