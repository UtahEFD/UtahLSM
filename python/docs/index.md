# UtahLSM Documentation

UtahLSM is a fast-response land-surface model for simulating the exchange of heat, moisture, and momentum between land and atmosphere. This site is built with MkDocs and uses `mkdocstrings` to render API documentation directly from the package docstrings in `utahlsm/`.

## Core Capabilities

- Offline or coupled workflows, with the offline driver living in `utahlsm_offline.py`
- Modular surface, soil, radiation, and canopy parameterizations
- Column-oriented state handling that scales from a single site to `nx × ny` column layouts
- NetCDF input and output plus schema-validated JSON namelists

## Documentation Map

- [Getting Started](getting-started.md): install the package, run the bundled GABLS3 case, and build the docs locally
- [Architecture](architecture.md): high-level diagrams of the model setup and timestep flow
- [Namelist Reference](namelist.md): section-by-section summary of the JSON configuration format
- [API Reference](reference.md): docstring-generated reference for the public model classes and utilities
- [Testing](testing.md): pytest commands and a docs build smoke test

## Docstring-Driven API Docs

The API reference is generated from the source tree rather than maintained by hand. Updating a class or function docstring in `utahlsm/` changes the rendered reference page on the next MkDocs build.

## Local Docs Workflow

Run the docs commands from the `python/` directory:

```bash
pip install -e ".[docs]"
mkdocs serve
```

For a non-interactive build:

```bash
mkdocs build
```

The expected published location for this project site is `https://utahefd.github.io/UtahLSM/`.
