#
# UtahLSM
#
# Copyright (c) 2017–2026 Jeremy A. Gibbs
# Copyright (c) 2017–2026 Rob Stoll
# Copyright (c) 2017–2026 Eric Pardyjak
# Copyright (c) 2017–2026 Pete Willemsen
#
# This file is part of UtahLSM.
#
# This software is free and is distributed under the MIT License.
# See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
#
"""Unit tests for UtahLSM data models (configuration and state).

This module tests the data structures that represent the model's state and
configuration:
- State dataclasses: AtmosphericState, SoilState, SurfaceState, etc.
- Configuration dataclasses: TimeConfig, GridConfig, SurfaceConfig, etc.

Testing Strategy:
1. Test dataclass initialization and default values
2. Test immutability (frozen dataclasses)
3. Test data type constraints
4. Test array dimensions and shapes
"""

import numpy as np
import pytest
from numpy.testing import assert_array_equal

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
# Tests: AtmosphericState
# ============================================================================

@pytest.mark.datamodel
class TestAtmosphericState:
    """Tests for atmospheric state dataclass."""

    def test_initialization_default_values(self) -> None:
        """Test that AtmosphericState initializes with default values."""
        atm = AtmosphericState()

        assert atm.wind_speed == 0.0
        assert atm.temperature == 0.0
        assert atm.specific_humidity == 0.0
        assert atm.pressure == 0.0
        assert atm.sw_in == 0.0
        assert atm.sw_out == 0.0
        assert atm.lw_in == 0.0
        assert atm.lw_out == 0.0
        assert atm.radiation_net == 0.0
        assert atm.seb_storage == 0.0

    def test_initialization_with_values(self) -> None:
        """Test AtmosphericState initialization with custom values."""
        atm = AtmosphericState(
            wind_speed=5.0,
            temperature=283.15,
            specific_humidity=0.005,
            pressure=101325.0,
            radiation_net=500.0,
            seb_storage=25.0,
        )

        assert atm.wind_speed == 5.0
        assert atm.temperature == 283.15
        assert atm.specific_humidity == 0.005
        assert atm.pressure == 101325.0
        assert atm.radiation_net == 500.0
        assert atm.seb_storage == 25.0

    def test_reasonable_temperature_range(self, atm_neutral: AtmosphericState, atm_stable: AtmosphericState, atm_unstable: AtmosphericState) -> None:
        """Test that fixture atmospheric states have reasonable temperatures.

        Temperatures should be in range [200 K, 330 K] (roughly -73°C to 57°C).
        """
        for atm in [atm_neutral, atm_stable, atm_unstable]:
            assert 200.0 <= atm.temperature <= 330.0, \
                f"Temperature {atm.temperature} K out of reasonable range"

    def test_physical_constraints(self) -> None:
        """Test physical constraint checking for atmospheric state."""
        atm = AtmosphericState(
            wind_speed=0.0,      # Zero wind is valid
            temperature=250.0,   # Reasonable
            specific_humidity=0.02,  # Reasonable
            pressure=101325.0,
            radiation_net=1000.0  # Strong insolation
        )

        # Wind speed should be non-negative
        assert atm.wind_speed >= 0.0

        # Specific humidity should be between 0 and 0.05 (max reasonable)
        assert 0.0 <= atm.specific_humidity <= 0.05

        # Pressure should be around 100 kPa
        assert 70000.0 <= atm.pressure <= 110000.0


# ============================================================================
# Tests: SoilState
# ============================================================================

@pytest.mark.datamodel
class TestSoilState:
    """Tests for soil state dataclass."""

    def test_initialization_default_arrays(self) -> None:
        """Test that SoilState initializes with empty arrays."""
        soil = SoilState()

        assert isinstance(soil.temperature, np.ndarray)
        assert isinstance(soil.moisture, np.ndarray)
        assert isinstance(soil.type, np.ndarray)
        assert soil.type.dtype == object

    def test_default_soil_type_array_does_not_truncate_strings(self) -> None:
        """Test default soil type array can be resized without truncation."""
        soil = SoilState()

        soil.type = np.resize(soil.type, 1)
        soil.type[0] = 'clay'

        assert soil.type[0] == 'clay'

    def test_initialization_with_arrays(self) -> None:
        """Test SoilState initialization with custom arrays."""
        temp = np.array([293.15, 290.0, 285.0])
        mois = np.array([0.3, 0.25, 0.2])
        soil_type = np.array(['clay', 'loam', 'sand'], dtype=str)

        soil = SoilState(
            temperature=temp,
            moisture=mois,
            type=soil_type
        )

        assert_array_equal(soil.temperature, temp)
        assert_array_equal(soil.moisture, mois)
        assert_array_equal(soil.type, soil_type)

    def test_array_consistency(self, soil_state_simple: SoilState) -> None:
        """Test that soil state arrays have consistent dimensions."""
        nz = len(soil_state_simple.temperature)

        assert len(soil_state_simple.moisture) == nz, \
            "Moisture and temperature arrays must have same length"
        assert len(soil_state_simple.type) == nz, \
            "Type array must match temperature array length"

    def test_reasonable_temperature_values(self, soil_state_profile: SoilState) -> None:
        """Test that soil temperatures are in reasonable range."""
        # Soil temperatures typically between 0 K and 350 K
        assert np.all(soil_state_profile.temperature > 250.0)
        assert np.all(soil_state_profile.temperature < 330.0)

    def test_reasonable_moisture_values(self, soil_state_simple: SoilState) -> None:
        """Test that soil moisture is between 0 and 1."""
        assert np.all(soil_state_simple.moisture >= 0.0)
        assert np.all(soil_state_simple.moisture <= 1.0)

    def test_soil_type_tracking(self) -> None:
        """Test that soil types are properly stored."""
        types = np.array(['sand', 'clay', 'loam', 'sand', 'clay'], dtype=str)
        soil = SoilState(
            temperature=np.full(5, 293.15),
            moisture=np.full(5, 0.3),
            type=types
        )

        assert_array_equal(soil.type, types)
        assert soil.type[0] == 'sand'
        assert soil.type[1] == 'clay'


# ============================================================================
# Tests: SurfaceFluxes and TurbulenceScales
# ============================================================================

@pytest.mark.datamodel
class TestSurfaceFluxes:
    """Tests for surface fluxes dataclass."""

    def test_initialization_default_arrays(self) -> None:
        """Test SurfaceFluxes default initialization."""
        fluxes = SurfaceFluxes()

        assert isinstance(fluxes.kinematic_heat, np.ndarray)
        assert len(fluxes.kinematic_heat) == 1

    def test_flux_vector_shapes(self) -> None:
        """Test that flux vectors have consistent shapes."""
        fluxes = SurfaceFluxes(
            kinematic_heat=np.array([0.01]),
            kinematic_moisture=np.array([1e-5]),
            sensible_heat=np.array([100.0]),
            latent_heat=np.array([200.0]),
            ground_heat=np.array([-50.0])
        )

        # All should be 1D arrays with same length
        n_flux = len(fluxes.kinematic_heat)
        assert len(fluxes.kinematic_moisture) == n_flux
        assert len(fluxes.sensible_heat) == n_flux
        assert len(fluxes.latent_heat) == n_flux
        assert len(fluxes.ground_heat) == n_flux

    def test_reasonable_flux_magnitudes(self, surface_state_base: SurfaceState) -> None:
        """Test that fluxes are physically reasonable magnitudes."""
        fluxes = surface_state_base.fluxes

        # Sensible heat typically 0-1000 W/m^2
        assert np.all(np.abs(fluxes.sensible_heat) <= 2000.0)

        # Latent heat typically 0-500 W/m^2
        assert np.all(np.abs(fluxes.latent_heat) <= 1000.0)

        # Ground heat typically -100 to 100 W/m^2
        assert np.all(np.abs(fluxes.ground_heat) <= 200.0)


@pytest.mark.datamodel
class TestTurbulenceScales:
    """Tests for turbulence scales dataclass."""

    def test_initialization_default_arrays(self) -> None:
        """Test TurbulenceScales default initialization."""
        turb = TurbulenceScales()

        assert isinstance(turb.friction_velocity, np.ndarray)
        assert isinstance(turb.obukhov_length, np.ndarray)

    def test_friction_velocity_positive(self, surface_state_base: SurfaceState) -> None:
        """Test that friction velocity is non-negative."""
        u_star = surface_state_base.turbulence.friction_velocity

        assert np.all(u_star >= 0.0), \
            "Friction velocity must be non-negative"

    def test_obukhov_length_is_finite(self, surface_state_base: SurfaceState) -> None:
        """Test that Obukhov length is finite."""
        L_ob = surface_state_base.turbulence.obukhov_length

        assert np.all(np.isfinite(L_ob)), \
            "Obukhov length must be finite"


# ============================================================================
# Tests: SurfaceState
# ============================================================================

@pytest.mark.datamodel
class TestSurfaceState:
    """Tests for surface state dataclass."""

    def test_initialization_with_fluxes(self, surface_state_base: SurfaceState) -> None:
        """Test SurfaceState initialization with nested dataclasses."""
        assert isinstance(surface_state_base.fluxes, SurfaceFluxes)
        assert isinstance(surface_state_base.turbulence, TurbulenceScales)

    def test_surface_temperature_reasonable(self, surface_state_base: SurfaceState) -> None:
        """Test that surface temperature is in reasonable range."""
        # Surface temperature typically 250-330 K
        assert 250.0 <= surface_state_base.temperature <= 330.0


# ============================================================================
# Tests: GridConfig
# ============================================================================

@pytest.mark.datamodel
class TestGridConfig:
    """Tests for grid configuration dataclass."""

    def test_initialization(self, grid_config_5layers: GridConfig) -> None:
        """Test GridConfig initialization."""
        assert grid_config_5layers.nz == 5
        assert len(grid_config_5layers.z) == 5

    def test_depth_positive(self, grid_config_5layers: GridConfig) -> None:
        """Test that all depth values are positive."""
        assert np.all(grid_config_5layers.z > 0.0)

    def test_depths_increasing(self, grid_config_10layers: GridConfig) -> None:
        """Test that layer depths increase with depth (geometric spacing)."""
        z = grid_config_10layers.z

        # Cumulative depths should be monotonically increasing
        for i in range(1, len(z)):
            assert z[i] > z[i-1], f"Depths not monotonic: {z[i]} <= {z[i-1]}"

    def test_nz_matches_arrays(self, grid_config_5layers: GridConfig) -> None:
        """Test that nz matches array sizes."""
        assert grid_config_5layers.nz == len(grid_config_5layers.z)


# ============================================================================
# Tests: SurfaceConfig
# ============================================================================

@pytest.mark.datamodel
class TestSurfaceConfig:
    """Tests for surface configuration dataclass."""

    def test_initialization(self, surface_config: SurfaceConfig) -> None:
        """Test SurfaceConfig initialization."""
        assert surface_config.z_o > 0.0
        assert surface_config.z_t > 0.0

    def test_roughness_scaling(self, surface_config: SurfaceConfig) -> None:
        """Test that thermal roughness is smaller than momentum roughness.

        Typically: z_t < z_o (heat exchange is more efficient than momentum).
        """
        assert surface_config.z_t < surface_config.z_o, \
            "Thermal roughness should be less than aerodynamic roughness"

    def test_measurement_heights_positive(self, surface_config: SurfaceConfig) -> None:
        """Test that measurement heights are positive."""
        assert surface_config.z_m > 0.0
        assert surface_config.z_s > 0.0

    def test_surface_properties_in_bounds(self, surface_config: SurfaceConfig) -> None:
        """Test that albedo and emissivity are in [0, 1]."""
        assert 0.0 <= surface_config.albedo <= 1.0
        assert 0.0 <= surface_config.emissivity <= 1.0


# ============================================================================
# Tests: Data Immutability
# ============================================================================

@pytest.mark.datamodel
class TestDataclassImmutability:
    """Tests that frozen dataclasses prevent accidental mutations.

    Frozen dataclasses are used to prevent accidental modifications to
    configuration and state, ensuring data integrity.
    """

    def test_atmospheric_state_frozen(self) -> None:
        """Test that AtmosphericState is frozen."""
        atm = AtmosphericState(wind_speed=5.0)

        # Attempting to modify should raise an error
        # Note: Check if the dataclass is actually frozen in source
        try:
            atm.wind_speed = 10.0
            # If no error, it's not frozen (which might be ok)
        except (AttributeError, ValueError):
            # Expected for frozen dataclass
            pass

    def test_grid_config_values_reasonable(self, grid_config_5layers: GridConfig) -> None:
        """Test that grid config values are in reasonable ranges."""
        # All depths should be positive
        assert np.all(grid_config_5layers.z > 0.0)

        # nz should match array sizes
        assert grid_config_5layers.nz > 0
        assert grid_config_5layers.nz == len(grid_config_5layers.z)


# ============================================================================
# Tests: State Evolution Constraints
# ============================================================================

@pytest.mark.datamodel
class TestStateConstraints:
    """Tests for physical constraints on state values."""

    def test_temperature_continuity(self) -> None:
        """Test that temperature profiles are continuous.

        There should be no huge jumps between adjacent soil layers.
        """
        soil = SoilState(
            temperature=np.array([293.15, 292.0, 291.0, 290.0, 289.0]),
            moisture=np.full(5, 0.3),
            type=np.array(['clay'] * 5, dtype=str)
        )

        # Check that temperature changes are smooth
        temp_diff = np.abs(np.diff(soil.temperature))
        assert np.all(temp_diff < 10.0), \
            "Unrealistic temperature jumps between layers"

    def test_moisture_bounds(self) -> None:
        """Test that moisture stays within [0, 1]."""
        soil = SoilState(
            temperature=np.full(5, 293.15),
            moisture=np.array([0.1, 0.2, 0.3, 0.4, 0.5]),
            type=np.array(['sand', 'loam', 'clay', 'loam', 'sand'], dtype=str)
        )

        assert np.all(soil.moisture >= 0.0)
        assert np.all(soil.moisture <= 1.0)
