# Architecture

UtahLSM separates input loading, timestep orchestration, physics parameterizations, and NetCDF output. The main `UtahLSM` class coordinates those pieces while concrete physics modules are selected from the namelist.

## High-Level Structure

```mermaid
flowchart TB
    subgraph entry[Entry Points]
        Driver[utahlsm_offline.py]
        Host[External host model]
    end

    subgraph input[Configuration and Input]
        NL[lsm_namelist.json]
        INIT[lsm_init.nc]
        FORCE[lsm_offline.nc]
        SOIL[data/soil/*.json]
        Input[utahlsm.util.io.input.Input]
    end

    subgraph state[Typed State and Config]
        Config[data_models config dataclasses]
        Forcing[ForcingData]
        Initial[Initial SoilState]
    end

    subgraph core[Core Orchestrator]
        Model[utahlsm.core.UtahLSM]
        Coupling[Surface coupling loop]
        Diffusion[Soil heat and moisture diffusion]
    end

    subgraph physics[Physics Modules]
        Surface[SurfaceMOST]
        Soil[BrooksCorey / Campbell / VanGenuchten]
        Radiation[RadBasic or disabled]
        Canopy[Jarvis canopy or bare soil]
        Thermo[thermo.py helpers]
    end

    subgraph output[Output]
        Output[utahlsm.util.io.output.Output]
        NC[NetCDF results]
    end

    Driver --> Input
    Host --> Model
    NL --> Input
    INIT --> Input
    FORCE --> Input
    SOIL --> Input
    Input --> Config
    Input --> Forcing
    Input --> Initial
    Config --> Model
    Forcing --> Model
    Initial --> Model
    Model --> Coupling
    Coupling --> Surface
    Coupling --> Soil
    Coupling --> Radiation
    Coupling --> Canopy
    Surface --> Thermo
    Soil --> Thermo
    Radiation --> Thermo
    Coupling --> Diffusion
    Diffusion --> Output
    Model --> Output
    Output --> NC
```

## Timestep Flow

```mermaid
flowchart TD
    Start(["Start"]) --> Load["Load namelist, initial state, and forcing"]
    Load --> Build["Instantiate Input, Output, and UtahLSM"]
    Build --> Loop{"More forcing records?"}
    Loop -->|Yes| Update["Update model with forcing"]
    Update --> SurfaceState["Refresh surface state from top soil layer"]
    SurfaceState --> Radiation["Optional radiation update"]
    Radiation --> Run["Run core physics"]
    Run --> Coupled["Iterative SEB and SMB surface coupling"]
    Coupled --> Heat["Soil heat diffusion solve"]
    Heat --> Moisture["Soil moisture diffusion solve"]
    Moisture --> Save["Write timestep output"]
    Save --> Loop
    Loop -->|No| Close["Close NetCDF output"]
    Close --> End(["End"])
```

## Model Selection from the Namelist

| Component | Selector | Current implementations |
| --- | --- | --- |
| Surface layer | `surface.model` | `1 = SurfaceMOST` |
| Soil hydrology | `soil.model` | `1 = BrooksCorey`, `2 = Campbell`, `3 = VanGenuchten` |
| Radiation | `radiation.model` | `0 = disabled`, `1 = RadBasic` |
| Canopy | `canopy.model` | `"none"`, `"jarvis"` |

## Design Notes

- `Input` is responsible for schema validation, NetCDF loading, and assembling typed dataclasses before the core model runs.
- `UtahLSM` owns the mutable runtime state and invokes the physics modules in the correct order.
- Soil, surface, radiation, and canopy implementations live behind abstract base classes so new parameterizations can be added without changing the timestep driver.
- The API reference page is generated directly from these docstrings, so the code remains the source of truth for method behavior.
