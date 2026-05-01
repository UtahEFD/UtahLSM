#
# UtahLSM
#
# Copyright (c) 2017–2025 Jeremy A. Gibbs
# Copyright (c) 2017–2025 Rob Stoll
# Copyright (c) 2017–2025 Eric Pardyjak
# Copyright (c) 2017–2025 Pete Willemsen
#
# This file is part of UtahLSM.
#
# This software is free and is distributed under the MIT License.
# See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
#
"""Pytest fixtures and configuration for UtahLSM testing.

This module provides reusable test fixtures and utilities for the UtahLSM
test suite. Fixtures defined here can be used across all test files.

Key fixtures:
- Atmospheric states (neutral, stable, unstable conditions)
- Soil states (various moisture/temperature profiles)
- Configuration dataclasses
"""

import numpy as np
import pytest

from utahlsm.data_models import (
    AtmosphericState,
    GridConfig,
    SoilState,
    SurfaceConfig,
    SurfaceFluxes,
    SurfaceState,
    TurbulenceScales,
)

# ============================================================================
# Atmospheric State Fixtures
# ============================================================================

@pytest.fixture
def atm_neutral():
    """Neutral atmospheric conditions (typical for testing).

    Returns:
        AtmosphericState with moderate wind, typical temperature,
        and neutral radiation balance.
    """
    return AtmosphericState(
        wind_speed=5.0,              # m/s
        temperature=283.15,          # K (~10°C)
        specific_humidity=0.005,     # kg/kg
        pressure=101325.0,           # Pa (standard sea level)
        sw_in=500.0,                 # W/m^2
        radiation_net=500.0          # W/m^2
    )


@pytest.fixture
def atm_stable():
    """Stable atmospheric conditions (cold, calm, strong inversion).

    Returns:
        AtmosphericState with low wind and temperature inversion.
    """
    return AtmosphericState(
        wind_speed=2.0,              # m/s (calm)
        temperature=275.15,          # K (cold, ~2°C)
        specific_humidity=0.003,     # kg/kg (dry)
        pressure=101325.0,           # Pa
        sw_in=100.0,                 # W/m^2
        radiation_net=100.0          # W/m^2 (weak radiation)
    )


@pytest.fixture
def atm_unstable():
    """Unstable atmospheric conditions (warm, calm, strong convection).

    Returns:
        AtmosphericState with high temperature and weak wind.
    """
    return AtmosphericState(
        wind_speed=1.0,              # m/s (weak wind)
        temperature=298.15,          # K (warm, ~25°C)
        specific_humidity=0.015,     # kg/kg (moist)
        pressure=101325.0,           # Pa
        sw_in=800.0,                 # W/m^2
        radiation_net=800.0          # W/m^2 (strong radiation)
    )


# ============================================================================
# Soil State Fixtures
# ============================================================================

@pytest.fixture
def soil_state_simple():
    """Simple soil state with uniform properties.

    Returns:
        SoilState with 5 layers, uniform temperature and moisture.
    """
    nz = 5
    return SoilState(
        temperature=np.full(nz, 293.15),  # 20°C everywhere
        moisture=np.full(nz, 0.3),         # 30% saturation
        type=np.array(['clay'] * nz)       # Same soil type
    )


@pytest.fixture
def soil_state_profile():
    """Soil state with realistic vertical profile.

    Returns:
        SoilState with decreasing temperature with depth,
        and varying moisture content.
    """
    nz = 10
    # Temperature decreases with depth (realistic gradient)
    temp = np.linspace(293.15, 288.15, nz)  # 20°C at surface, 15°C at depth
    # Moisture is higher near surface (drying with depth)
    mois = np.linspace(0.35, 0.25, nz)

    return SoilState(
        temperature=temp,
        moisture=mois,
        type=np.array(['clay'] * nz)
    )


@pytest.fixture
def soil_state_dry():
    """Dry soil state (low moisture).

    Returns:
        SoilState with low moisture content throughout.
    """
    nz = 5
    return SoilState(
        temperature=np.full(nz, 293.15),  # 20°C
        moisture=np.full(nz, 0.15),       # 15% (dry)
        type=np.array(['sand'] * nz)
    )


@pytest.fixture
def soil_state_saturated():
    """Saturated soil state (high moisture).

    Returns:
        SoilState with high moisture content throughout.
    """
    nz = 5
    return SoilState(
        temperature=np.full(nz, 293.15),  # 20°C
        moisture=np.full(nz, 0.45),       # 45% (wet)
        type=np.array(['clay'] * nz)
    )


# ============================================================================
# Surface State Fixtures
# ============================================================================

@pytest.fixture
def surface_state_base():
    """Basic surface state for testing.

    Returns:
        SurfaceState with default values and single-layer fluxes.
    """
    return SurfaceState(
        temperature=293.15,           # 20°C
        fluxes=SurfaceFluxes(
            kinematic_heat=np.array([0.01]),
            kinematic_moisture=np.array([1e-5]),
            sensible_heat=np.array([100.0]),
            latent_heat=np.array([200.0]),
            ground_heat=np.array([-50.0])
        ),
        turbulence=TurbulenceScales(
            friction_velocity=np.array([0.3]),
            obukhov_length=np.array([10.0])
        )
    )


# ============================================================================
# Configuration Fixtures
# ============================================================================

@pytest.fixture
def grid_config_5layers():
    """Grid configuration with 5 soil layers.

    Returns:
        GridConfig with typical soil depth and layer spacing.
    """
    # Cumulative depths for 5 layers
    z = np.array([0.05, 0.20, 0.55, 1.30, 2.80])  # Cumulative layer depths (m)

    return GridConfig(
        nx=1,
        ny=1,
        nz=5,
        z=z
    )


@pytest.fixture
def grid_config_10layers():
    """Grid configuration with 10 soil layers (finer resolution).

    Returns:
        GridConfig with finer vertical resolution.
    """
    nz = 10
    # Geometric spacing: shallower layers near surface
    z = np.zeros(nz)
    z[0] = 0.02
    for i in range(1, nz):
        z[i] = z[i-1] * 1.5  # 1.5x increase per layer

    return GridConfig(
        nx=1,
        ny=1,
        nz=nz,
        z=z
    )


@pytest.fixture
def surface_config():
    """Surface configuration for MOST calculations.

    Returns:
        SurfaceConfig with typical roughness lengths.
    """
    return SurfaceConfig(
        z_o=0.1,              # Aerodynamic roughness length [m]
        z_t=0.001,            # Thermal roughness length [m]
        z_m=10.0,             # Measurement height for wind [m]
        z_s=2.0,              # Measurement height for temp/humidity [m]
        albedo=0.23,          # Surface albedo
        emissivity=0.98,      # Surface emissivity
        model=1               # Surface layer model ID
    )


# ============================================================================
# Helper Functions for Test Data Generation
# ============================================================================

@pytest.fixture
def create_tridiagonal_system():
    """Factory fixture to create test tridiagonal systems.

    Returns:
        A function that creates tridiagonal systems with known solutions.
    """
    def _create_system(n: int, condition='well-conditioned') -> tuple:
        """Create a tridiagonal system Ax = b.

        Args:
            n: Size of the matrix
            condition: 'well-conditioned', 'diagonal-dominant', or 'ill-conditioned'

        Returns:
            Tuple of (a, b, c, r, x_expected) where the system Ax=r
            has solution x=x_expected.
        """
        if condition == 'well-conditioned':
            # Strong diagonal dominance
            a = np.full(n, -1.0)
            b = np.full(n, 4.0)
            c = np.full(n, -1.0)
        elif condition == 'diagonal-dominant':
            # Just barely diagonal dominant
            a = np.full(n, -0.4)
            b = np.full(n, 1.0)
            c = np.full(n, -0.4)
        else:  # ill-conditioned
            a = np.full(n, -0.6)
            b = np.full(n, 1.0)
            c = np.full(n, -0.6)

        # Create solution x, then compute b such that Ax = b
        x_expected = np.linspace(1.0, 2.0, n)

        # Compute r = Ax
        r = np.zeros(n)
        r[0] = b[0] * x_expected[0] + c[0] * x_expected[1]
        for i in range(1, n-1):
            r[i] = (a[i] * x_expected[i-1] +
                   b[i] * x_expected[i] +
                   c[i] * x_expected[i+1])
        r[n-1] = a[n-1] * x_expected[n-2] + b[n-1] * x_expected[n-1]

        return a, b, c, r, x_expected

    return _create_system


@pytest.fixture
def create_root_function():
    """Factory fixture to create test functions for root-finding.

    Returns:
        A function that creates test functions with known roots.
    """
    def _create_function(root: float, type_: str = 'quadratic'):
        """Create a test function with a known root.

        Args:
            root: The location of the known root
            type_: 'quadratic', 'cubic', 'sine', 'rational'

        Returns:
            Tuple of (f, bracket) where f is the function
            and bracket=[a, b] surrounds the root.
        """
        if type_ == 'quadratic':
            # f(x) = (x - root)^2 - 1
            f = lambda x: (x - root)**2 - 1
            bracket = [root - 2, root + 2]
        elif type_ == 'cubic':
            # f(x) = (x - root)^3 - 1
            f = lambda x: (x - root)**3 - 1
            bracket = [root - 2, root + 2]
        elif type_ == 'sine':
            # f(x) = sin(x - root)
            f = lambda x: np.sin(x - root)
            bracket = [root - np.pi/2, root + np.pi/2]
        else:  # rational
            # f(x) = 1/(x - root) - 1
            f = lambda x: 1 / (x - root + 1e-6) - 1
            bracket = [root - 1, root + 1]

        return f, bracket

    return _create_function


# ============================================================================
# Assertion Helpers
# ============================================================================

@pytest.fixture
def assert_physically_reasonable():
    """Fixture providing assertion helpers for physical constraints.

    Returns:
        A function that checks physical constraints.
    """
    def _check(value: float, name: str, min_: float = None, max_: float = None):
        """Check that a value is within physically reasonable bounds.

        Args:
            value: The value to check
            name: Name of the quantity (for error messages)
            min_: Minimum reasonable value
            max_: Maximum reasonable value
        """
        if min_ is not None and value < min_:
            raise AssertionError(
                f"{name} = {value} is below minimum {min_}")
        if max_ is not None and value > max_:
            raise AssertionError(
                f"{name} = {value} is above maximum {max_}")

    return _check
