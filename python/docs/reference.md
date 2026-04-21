# API Reference

The sections below are rendered from the `utahlsm` docstrings by `mkdocstrings`.

## Package Overview

::: utahlsm
    options:
      members: false
      show_root_heading: false

## Core Model

::: utahlsm.core.UtahLSM

## Data Models

::: utahlsm.data_models

## Input and Output

::: utahlsm.util.io.input.Input

::: utahlsm.util.io.output.Output

::: utahlsm.util.io.soil_properties_loader.SoilPropertiesLoader

::: utahlsm.util.io.logging_helper

## Physics Interfaces

::: utahlsm.physics.surface.sfc.Surface

::: utahlsm.physics.soil.soil.SoilProperties

::: utahlsm.physics.soil.soil.Soil

::: utahlsm.physics.radiation.radiation.Radiation

::: utahlsm.physics.canopy.canopy.Canopy

## Physics Implementations

::: utahlsm.physics.surface.sfc_most.SurfaceMOST

::: utahlsm.physics.soil.soil_brookscorey.BrooksCorey

::: utahlsm.physics.soil.soil_campbell.Campbell

::: utahlsm.physics.soil.soil_vangenuchten.VanGenuchten

::: utahlsm.physics.radiation.rad_basic.RadBasic

::: utahlsm.physics.canopy.canopy_jarvis.CanopyJarvis

## Factory Modules

::: utahlsm.physics.surface.factory

::: utahlsm.physics.soil.factory

::: utahlsm.physics.radiation.factory

::: utahlsm.physics.canopy.factory

## Thermodynamics and Numerics

::: utahlsm.physics.thermo

::: utahlsm.util.solvers

## Constants and Exceptions

::: utahlsm.util.constants
    options:
      show_if_no_docstring: true
      members: true

::: utahlsm.exceptions
