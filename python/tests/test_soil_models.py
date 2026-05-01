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
"""Unit tests for soil physics models in UtahLSM.

This module tests soil-specific physics calculations:
- Water potential relationships (psi-theta curves)
- Hydraulic conductivity functions
- Thermal conductivity and diffusivity
- Soil property constraints

Testing Strategy:
1. Verify physics constraints (monotonicity, sign, bounds)
2. Test with known analytical solutions where available
3. Test with realistic soil parameter ranges
4. Verify numerical stability
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

# Import soil models and dependencies
from utahlsm.physics.soil.soil_brookscorey import BrooksCorey
from utahlsm.physics.soil.soil_campbell import Campbell
from utahlsm.physics.soil.soil_vangenuchten import VanGenuchten
from utahlsm.util.io.soil_properties_loader import SoilPropertiesLoader

# ============================================================================
# Fixtures: Load Real Soil Properties
# ============================================================================

@pytest.fixture(
    params=['clapp-hornberger', 'cosby', 'rawls-brakensiek', 'cabauw-heinen']
)
def soil_props(request):
    """Load bundled soil properties dataset.

    Parametrized fixture that yields a tuple of (properties_dict, dataset_name)
    for all available bundled datasets:
    - clapp-hornberger: Clapp & Hornberger model parameters
    - cosby: Cosby et al. pedotransfer functions
    - rawls-brakensiek: Rawls-Brakensiek parameterization
    - cabauw-heinen: Cabauw experimental site data

    Each test using this fixture will run once per dataset.
    """
    properties_dict = SoilPropertiesLoader.load(request.param)
    return (properties_dict, request.param)


@pytest.fixture
def soil_types_to_test(soil_props):
    """Soil type names appropriate for the current soil_props dataset.

    Returns common soil types for most datasets, or site-specific types
    for cabauw-heinen. This fixture is coupled to soil_props and will
    automatically select the right types.
    """
    properties_dict, dataset_name = soil_props

    # cabauw-heinen uses site-specific soil names
    if dataset_name == 'cabauw-heinen':
        # Return some available types from the cabauw dataset
        available = list(properties_dict.keys())
        # Return up to 3 types
        return available[:3]

    # Standard USDA soil type names
    return ['sand', 'loam', 'clay']


# ============================================================================
# Tests: Water Potential Functions (psi-theta curves)
# ============================================================================

@pytest.mark.soil
class TestWaterPotential:
    """Tests for water potential calculations across soil models."""

    def test_brookscorey_water_potential_monotonicity(self, soil_props, soil_types_to_test):
        """Test that water potential is roughly monotonic with moisture.

        Physical principle: wetter soil should have less negative (higher) water potential.
        In the Brooks-Corey model near saturation, the relationship may be weakly monotonic.
        """
        properties_dict, dataset_name = soil_props
        model = BrooksCorey(properties_dict, soil_types_to_test, dataset_name)
        nz = len(model.properties.porosity)
        layer = 0  # Test first layer

        # Test moisture levels from mid-range to saturation (avoid both extremes)
        theta_min = model.properties.residual[layer] + 0.1 * (model.properties.porosity[layer] - model.properties.residual[layer])
        theta = np.linspace(
            theta_min,
            model.properties.porosity[layer],
            10
        )

        psi = np.array([model.water_potential(t, level=layer) for t in theta])

        # Water potential should be monotonically INCREASING (becoming less negative)
        # as moisture increases from mid-range to saturation
        # Allow for small tolerance due to numerical precision
        for i in range(len(psi) - 1):
            assert psi[i] <= psi[i+1] + 1e-6, \
                f"Water potential not monotonic at indices {i},{i+1}: {psi[i]}, {psi[i+1]}"

    def test_campbell_water_potential_saturation(self, soil_props, soil_types_to_test):
        """Test that water potential at saturation is near zero (or specified value).

        At saturation (theta = porosity), water potential should be at maximum
        (least negative).
        """
        properties_dict, dataset_name = soil_props
        model = Campbell(properties_dict, soil_types_to_test, dataset_name)
        nz = len(model.properties.porosity)
        layer = 0

        # At saturation
        psi_sat = model.water_potential(model.properties.porosity[layer], level=layer)

        # Should be close to zero or slightly negative
        assert psi_sat < 0.1, f"Water potential at saturation should be near 0: {psi_sat}"

    def test_vangenuchten_water_potential_residual(self, soil_props, soil_types_to_test):
        """Test that water potential at residual moisture is very negative.

        At residual moisture, water is tightly bound to soil particles.
        """
        properties_dict, dataset_name = soil_props
        model = VanGenuchten(properties_dict, soil_types_to_test, dataset_name)
        nz = len(model.properties.porosity)
        layer = 0

        # At residual moisture
        psi_res = model.water_potential(model.properties.residual[layer], level=layer)

        # Should be very negative (strong suction)
        assert psi_res < -100.0, \
            f"Water potential at residual should be very negative: {psi_res}"

    def test_water_potential_bounds(self, soil_props, soil_types_to_test):
        """Test that water potential stays within physical bounds.

        For any moisture level, psi should remain between saturation and
        residual water potential extremes.
        """
        properties_dict, dataset_name = soil_props
        for ModelClass in [BrooksCorey, Campbell, VanGenuchten]:
            model = ModelClass(properties_dict, soil_types_to_test, dataset_name)
            layer = 0

            # Range of moisture values (avoid exactly at residual to prevent singularities)
            theta_test = np.linspace(
                model.properties.residual[layer] + 0.001,
                model.properties.porosity[layer],
                20
            )

            for theta in theta_test:
                psi = model.water_potential(theta, level=layer)

                # Water potential should be reasonable (not inf, not nan)
                assert np.isfinite(psi), \
                    f"{ModelClass.__name__}: Non-finite psi={psi} at theta={theta}"

                # Should be negative or very close to zero (at saturation may be near 0)
                assert psi <= 1e-6, \
                    f"{ModelClass.__name__}: Significantly positive psi={psi} at theta={theta}"

    def test_water_content_inverts_water_potential(self, soil_props, soil_types_to_test):
        """Test that θ(ψ(θ)) recovers the original moisture value."""
        properties_dict, dataset_name = soil_props
        for ModelClass in [BrooksCorey, Campbell, VanGenuchten]:
            model = ModelClass(properties_dict, soil_types_to_test, dataset_name)
            layer = 0
            theta_test = np.linspace(
                model.properties.residual[layer] + 0.01,
                model.properties.porosity[layer] - 0.01,
                10,
            )

            psi = np.array(
                [model.water_potential(theta, level=layer) for theta in theta_test]
            )
            theta_back = np.array(
                [model.water_content(psi_i, level=layer) for psi_i in psi]
            )

            assert_allclose(theta_back, theta_test, rtol=1e-5, atol=1e-6)

    def test_moisture_capacity_positive(self, soil_props, soil_types_to_test):
        """Test that dθ/dψ is finite and non-negative in the unsaturated range."""
        properties_dict, dataset_name = soil_props
        for ModelClass in [BrooksCorey, Campbell, VanGenuchten]:
            model = ModelClass(properties_dict, soil_types_to_test, dataset_name)
            layer = 0
            theta_test = np.linspace(
                model.properties.residual[layer] + 0.01,
                model.properties.porosity[layer] - 0.01,
                10,
            )
            psi = np.array(
                [model.water_potential(theta, level=layer) for theta in theta_test]
            )
            capacity = np.array(
                [model.moisture_capacity(psi_i, level=layer) for psi_i in psi]
            )

            assert np.all(np.isfinite(capacity)), \
                f"{ModelClass.__name__}: Non-finite moisture capacity"
            assert np.all(capacity >= 0.0), \
                f"{ModelClass.__name__}: Negative moisture capacity"


# ============================================================================
# Tests: Hydraulic Conductivity
# ============================================================================

@pytest.mark.soil
class TestHydraulicConductivity:
    """Tests for hydraulic conductivity calculations."""

    def test_conductivity_moisture_saturation(self, soil_props, soil_types_to_test):
        """Test that conductivity at saturation equals K_sat.

        At saturation (theta = porosity), moisture conductivity should equal
        the saturated value.
        """
        properties_dict, dataset_name = soil_props
        for ModelClass in [BrooksCorey, Campbell, VanGenuchten]:
            model = ModelClass(properties_dict, soil_types_to_test, dataset_name)
            layer = 0

            # Conductivity at saturation
            K_at_sat = model.conductivity_moisture(
                model.properties.porosity[layer], level=layer
            )
            K_sat = model.properties.K_sat[layer]

            # Should match (with small tolerance for numerical issues)
            assert_allclose(K_at_sat, K_sat, rtol=0.01)

    def test_conductivity_moisture_decreases_with_dryness(self, soil_props, soil_types_to_test):
        """Test that conductivity increases as soil wets.

        Physical constraint: K(theta) is monotonically increasing with theta
        """
        properties_dict, dataset_name = soil_props
        for ModelClass in [BrooksCorey, Campbell, VanGenuchten]:
            model = ModelClass(properties_dict, soil_types_to_test, dataset_name)
            layer = 0

            # Test moisture levels
            theta = np.linspace(
                model.properties.residual[layer] + 0.001,  # Avoid residual singularity
                model.properties.porosity[layer],
                15
            )

            K = np.array([model.conductivity_moisture(t, level=layer) for t in theta])

            # Conductivity should be monotonically increasing with moisture
            for i in range(len(K) - 1):
                assert K[i] <= K[i+1] * 1.01, \
                    f"{ModelClass.__name__}: K not monotonic at {i},{i+1}: {K[i]}, {K[i+1]}"

    def test_conductivity_moisture_positive(self, soil_props, soil_types_to_test):
        """Test that moisture conductivity is non-negative."""
        properties_dict, dataset_name = soil_props
        for ModelClass in [BrooksCorey, Campbell, VanGenuchten]:
            model = ModelClass(properties_dict, soil_types_to_test, dataset_name)
            layer = 0

            theta_test = np.linspace(
                model.properties.residual[layer] + 0.01,  # Avoid singularity at residual
                model.properties.porosity[layer],
                20
            )

            for theta in theta_test:
                K = model.conductivity_moisture(theta, level=layer)

                # Should be non-negative (very small values near residual can be ~0)
                assert K >= -1e-10, \
                    f"{ModelClass.__name__}: Significantly negative K={K} at theta={theta}"
                assert np.isfinite(K), \
                    f"{ModelClass.__name__}: Non-finite K={K} at theta={theta}"


# ============================================================================
# Tests: Thermal Properties
# ============================================================================

@pytest.mark.soil
class TestThermalProperties:
    """Tests for thermal conductivity and diffusivity."""

    def test_thermal_conductivity_positive(self, soil_props, soil_types_to_test):
        """Test that thermal conductivity is always positive."""
        properties_dict, dataset_name = soil_props
        for ModelClass in [BrooksCorey, Campbell, VanGenuchten]:
            model = ModelClass(properties_dict, soil_types_to_test, dataset_name)

            # Create full moisture profile
            nz = len(model.properties.porosity)
            theta_profile = np.linspace(0.05, 0.4, nz)

            K_th = model.conductivity_thermal(theta_profile)

            assert np.all(K_th > 0.0), \
                f"{ModelClass.__name__}: Non-positive K_th found"
            assert np.all(np.isfinite(K_th)), \
                f"{ModelClass.__name__}: Non-finite K_th found"

    def test_thermal_diffusivity_positive(self, soil_props, soil_types_to_test):
        """Test that thermal diffusivity is always positive.

        Thermal diffusivity = K_th / (density * specific_heat)
        """
        properties_dict, dataset_name = soil_props
        for ModelClass in [BrooksCorey, Campbell, VanGenuchten]:
            model = ModelClass(properties_dict, soil_types_to_test, dataset_name)

            # Create full moisture profile
            nz = len(model.properties.porosity)
            theta_profile = np.linspace(0.05, 0.4, nz)

            D_th = model.diffusivity_thermal(theta_profile)

            assert np.all(D_th > 0.0), \
                f"{ModelClass.__name__}: Non-positive D_th found"
            assert np.all(np.isfinite(D_th)), \
                f"{ModelClass.__name__}: Non-finite D_th found"

    def test_moisture_diffusivity_valid(self, soil_props, soil_types_to_test):
        """Test that moisture diffusivity returns valid values."""
        properties_dict, dataset_name = soil_props
        for ModelClass in [BrooksCorey, Campbell, VanGenuchten]:
            model = ModelClass(properties_dict, soil_types_to_test, dataset_name)

            # Create full moisture profile (well above residual to avoid singularities)
            nz = len(model.properties.porosity)
            theta_profile = np.linspace(0.15, 0.4, nz)

            D_m = model.diffusivity_moisture(theta_profile)

            # Should be finite (diffusivity can be very small but should be reasonable)
            assert np.all(np.isfinite(D_m)), \
                f"{ModelClass.__name__}: Non-finite D_m found"


# ============================================================================
# Tests: Soil Parameter Constraints
# ============================================================================

@pytest.mark.soil
class TestSoilParameterConstraints:
    """Tests for soil property bounds and physical consistency."""

    def test_residual_less_than_porosity(self, soil_props, soil_types_to_test):
        """Test that residual moisture < porosity for all layers."""
        properties_dict, dataset_name = soil_props
        for ModelClass in [BrooksCorey, Campbell, VanGenuchten]:
            model = ModelClass(properties_dict, soil_types_to_test, dataset_name)

            assert np.all(model.properties.residual < model.properties.porosity), \
                f"{ModelClass.__name__}: Residual >= Porosity in some layers"

    def test_porosity_in_physical_range(self, soil_props, soil_types_to_test):
        """Test that porosity is between 0 and 1."""
        properties_dict, dataset_name = soil_props
        for ModelClass in [BrooksCorey, Campbell, VanGenuchten]:
            model = ModelClass(properties_dict, soil_types_to_test, dataset_name)

            assert np.all(model.properties.porosity > 0.0)
            assert np.all(model.properties.porosity < 1.0)

    def test_saturated_conductivity_positive(self, soil_props, soil_types_to_test):
        """Test that K_sat is positive for all layers."""
        properties_dict, dataset_name = soil_props
        for ModelClass in [BrooksCorey, Campbell, VanGenuchten]:
            model = ModelClass(properties_dict, soil_types_to_test, dataset_name)

            assert np.all(model.properties.K_sat > 0.0), \
                f"{ModelClass.__name__}: Non-positive K_sat"

    def test_heat_capacity_positive(self, soil_props, soil_types_to_test):
        """Test that volumetric heat capacity is positive."""
        properties_dict, dataset_name = soil_props
        for ModelClass in [BrooksCorey, Campbell, VanGenuchten]:
            model = ModelClass(properties_dict, soil_types_to_test, dataset_name)

            assert np.all(model.properties.ci > 0.0), \
                f"{ModelClass.__name__}: Non-positive heat capacity"


# ============================================================================
# Tests: Model Consistency
# ============================================================================

@pytest.mark.soil
class TestModelConsistency:
    """Tests for consistency across different soil models."""

    def test_all_models_compute_same_layer_count(self, soil_props, soil_types_to_test):
        """Test that all models handle the same number of layers."""
        properties_dict, dataset_name = soil_props
        models = [
            BrooksCorey(properties_dict, soil_types_to_test, dataset_name),
            Campbell(properties_dict, soil_types_to_test, dataset_name),
            VanGenuchten(properties_dict, soil_types_to_test, dataset_name),
        ]

        nz_expected = len(models[0].properties.porosity)

        for model in models:
            assert len(model.properties.porosity) == nz_expected
            assert len(model.properties.residual) == nz_expected
            assert len(model.properties.K_sat) == nz_expected

    def test_moisture_in_physical_bounds(self, soil_props, soil_types_to_test):
        """Test that test moisture values stay in [residual, porosity] bounds."""
        properties_dict, dataset_name = soil_props
        model = BrooksCorey(properties_dict, soil_types_to_test, dataset_name)
        layer = 0

        # Moisture values should be bounded
        theta_test = np.linspace(
            model.properties.residual[layer],
            model.properties.porosity[layer],
            5
        )

        for theta in theta_test:
            assert theta >= model.properties.residual[layer] - 1e-10
            assert theta <= model.properties.porosity[layer] + 1e-10


# ============================================================================
# Tests: Numerical Stability
# ============================================================================

@pytest.mark.soil
class TestNumericalStability:
    """Tests for numerical stability near boundaries."""

    def test_no_nan_near_saturation(self, soil_props, soil_types_to_test):
        """Test that calculations don't produce NaN near saturation."""
        properties_dict, dataset_name = soil_props
        model = BrooksCorey(properties_dict, soil_types_to_test, dataset_name)
        layer = 0

        # Values very close to saturation
        porosity = model.properties.porosity[layer]
        theta_near_sat = np.linspace(porosity * 0.95, porosity * 0.9999, 5)

        for theta in theta_near_sat:
            psi = model.water_potential(theta, level=layer)
            K = model.conductivity_moisture(theta, level=layer)

            assert np.isfinite(psi), "NaN in water potential near saturation"
            assert np.isfinite(K), "NaN in conductivity near saturation"

    def test_no_inf_near_residual(self, soil_props, soil_types_to_test):
        """Test that calculations stay finite near residual moisture.

        Note: Some models may have singularities exactly at residual.
        """
        properties_dict, dataset_name = soil_props
        model = VanGenuchten(properties_dict, soil_types_to_test, dataset_name)
        layer = 0

        # Values slightly above residual (to avoid singularities)
        residual = model.properties.residual[layer]
        theta_near_res = np.linspace(residual * 1.1, residual * 1.5, 5)

        for theta in theta_near_res:
            if theta > residual:  # Only test above residual
                psi = model.water_potential(theta, level=layer)
                assert np.isfinite(psi), "Inf/NaN in water potential near residual"
