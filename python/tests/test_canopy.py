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
"""Unit tests for the canopy subsystem.

Covers:
- The ``Canopy`` base class' root-distribution construction and
  root-weighted reduction.
- The ``CanopyJarvis`` multiplicative stress-factor response (bounds,
  monotonicity, wilt/FC anchors).
- The ``get_canopy_model`` factory (bare-soil pass-through, scalar-to-
  per-column broadcasting, invalid-model error path).
- End-to-end partitioning inside ``_solve_most`` /
  ``_finalize_canopy_partition``: total LH = LH_soil + LH_veg, and the
  partition reduces to pure bare-soil when f_veg == 0 or r_s is huge.
"""

import logging
from types import SimpleNamespace
from typing import Any, cast

import numpy as np
import pytest
from numpy.typing import NDArray

from utahlsm.core import UtahLSM
from utahlsm.data_models import (
    AtmosphericState,
    CanopyConfig,
    CanopyState,
    SoilState,
    SurfaceState,
    TurbulenceScales,
)
from utahlsm.exceptions import NamelistError
from utahlsm.physics import thermo
from utahlsm.physics.canopy.canopy import Canopy
from utahlsm.physics.canopy.canopy_jarvis import CanopyJarvis
from utahlsm.physics.canopy.factory import get_canopy_model

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def z_layers() -> NDArray[np.float64]:
    """Soil node depths used throughout these tests (negative-downward)."""
    # 11 layers from 0 to -3 m, geometric-ish spacing.
    return -np.array([0.0, 0.02, 0.06, 0.13, 0.25, 0.45, 0.75, 1.20,
                      1.80, 2.40, 3.00])


@pytest.fixture
def canopy_params_single(z_layers: NDArray[np.float64]) -> dict[str, Any]:
    """Per-column parameters for a single-column grassland canopy."""
    ncol = 1
    return {
        "lai": np.full(ncol, 3.0),
        "veg_fraction": np.full(ncol, 0.95),
        "rooting_depth": np.full(ncol, 0.4),
        "beta": np.full(ncol, 0.965),
        "rs_min": np.full(ncol, 40.0),
        "rs_max": np.full(ncol, 5000.0),
        "rg_half": np.full(ncol, 30.0),
        "vpd_coef": np.full(ncol, 1.0e-4),
        "t_opt": np.full(ncol, 298.0),
        "t_coef": np.full(ncol, 1.6e-3),
        "r_ground": np.full(ncol, 0.0),
        "water_capacity_lai": np.full(ncol, 0.2),
        "wet_cooling_max": np.full(ncol, 3.0),
        "z": z_layers,
    }


@pytest.fixture
def jarvis_single(canopy_params_single: dict[str, Any]) -> CanopyJarvis:
    return CanopyJarvis(**canopy_params_single)


@pytest.fixture
def canopy_params_3col(z_layers: NDArray[np.float64]) -> dict[str, Any]:
    ncol = 3
    return {
        "lai": np.array([1.0, 3.0, 5.0]),
        "veg_fraction": np.array([0.2, 0.6, 1.0]),
        "rooting_depth": np.array([0.3, 0.5, 1.0]),
        "beta": np.full(ncol, 0.965),
        "rs_min": np.full(ncol, 40.0),
        "rs_max": np.full(ncol, 5000.0),
        "rg_half": np.full(ncol, 30.0),
        "vpd_coef": np.full(ncol, 1.0e-4),
        "t_opt": np.full(ncol, 298.0),
        "t_coef": np.full(ncol, 1.6e-3),
        "r_ground": np.full(ncol, 0.0),
        "water_capacity_lai": np.full(ncol, 0.2),
        "wet_cooling_max": np.full(ncol, 3.0),
        "z": z_layers,
    }


@pytest.fixture
def jarvis_3col(canopy_params_3col: dict[str, Any]) -> CanopyJarvis:
    return CanopyJarvis(**canopy_params_3col)


# ---------------------------------------------------------------------------
# Canopy base class: root distribution and reductions
# ---------------------------------------------------------------------------

@pytest.mark.canopy
class TestRootDistribution:
    """Covers the Jackson-1996 truncated root-fraction construction."""

    def test_columns_sum_to_one(self, jarvis_single: CanopyJarvis) -> None:
        rf = jarvis_single.root_fraction
        assert rf.shape == (11, 1)
        assert np.isclose(rf.sum(axis=0), 1.0).all()

    def test_positive_everywhere(self, jarvis_single: CanopyJarvis) -> None:
        assert np.all(jarvis_single.root_fraction >= 0.0)

    def test_roots_respect_rooting_depth(self, jarvis_single: CanopyJarvis, z_layers: NDArray[np.float64]) -> None:
        rf = jarvis_single.root_fraction[:, 0]
        # rooting_depth = 0.4m, so layers whose top is beyond 0.4m
        # should have zero fraction. Use node depths (positive) as a proxy.
        deep_layers = -z_layers > 0.5
        assert np.all(rf[deep_layers] == 0.0)

    def test_per_column_rooting_respected(self, jarvis_3col: CanopyJarvis, z_layers: NDArray[np.float64]) -> None:
        rf = jarvis_3col.root_fraction
        # All columns must still sum to 1.
        assert np.allclose(rf.sum(axis=0), 1.0)
        # Column 0 (0.3m) should have a more concentrated top-heavy
        # profile than column 2 (1.0m).
        assert rf[1, 0] > rf[1, 2]

    def test_degenerate_rooting_depth_zero(self, canopy_params_single: dict[str, Any]) -> None:
        params = dict(canopy_params_single)
        params['rooting_depth'] = np.array([0.0])
        canopy = CanopyJarvis(**params)
        # Should not raise; all-zero normalizer falls back to 1.0 and
        # produces an all-zero profile.
        assert np.all(canopy.root_fraction == 0.0)


@pytest.mark.canopy
class TestRootZoneMean:
    """Covers the weighted-mean utility used by f4(θ_root)."""

    def test_constant_field_preserved(self, jarvis_single: CanopyJarvis) -> None:
        nz = 11
        field = np.full(nz, 0.3)
        mean = jarvis_single.root_zone_mean(field)
        assert mean.shape == (1,)
        assert np.isclose(mean[0], 0.3)

    def test_reduces_across_columns(self, jarvis_3col: CanopyJarvis) -> None:
        field = np.full((11, 3), 0.25)
        mean = jarvis_3col.root_zone_mean(field)
        assert mean.shape == (3,)
        assert np.allclose(mean, 0.25)

    def test_weights_shallow_for_shallow_roots(
        self, canopy_params_single: dict[str, Any], z_layers: NDArray[np.float64]
    ) -> None:
        # Moisture profile with a dry top and wet bottom: shallow-rooted
        # column should see lower θ_root than a deep-rooted column.
        nz = len(z_layers)
        moisture = np.linspace(0.1, 0.4, nz)

        shallow = dict(canopy_params_single)
        shallow['rooting_depth'] = np.array([0.1])
        deep = dict(canopy_params_single)
        deep['rooting_depth'] = np.array([2.0])

        c_shallow = CanopyJarvis(**shallow)
        c_deep = CanopyJarvis(**deep)

        m_shallow = c_shallow.root_zone_mean(moisture)[0]
        m_deep = c_deep.root_zone_mean(moisture)[0]
        assert m_shallow < m_deep


# ---------------------------------------------------------------------------
# Jarvis stress functions
# ---------------------------------------------------------------------------

@pytest.mark.canopy
class TestJarvisStressFunctions:
    """Individual f1..f4 factors and their physical limits."""

    def test_f_radiation_night_zero(self, jarvis_single: CanopyJarvis) -> None:
        f = jarvis_single._f_radiation(np.array([-50.0]))
        assert np.all(f == 0.0)

    def test_f_radiation_saturates(self, jarvis_single: CanopyJarvis) -> None:
        f_low = jarvis_single._f_radiation(np.array([10.0]))
        f_high = jarvis_single._f_radiation(np.array([2000.0]))
        assert f_low[0] < f_high[0]
        assert f_high[0] < 1.0  # saturating form never reaches 1
        assert f_high[0] > 0.9

    def test_f_temperature_optimum(self, jarvis_single: CanopyJarvis) -> None:
        leaf_T = np.array([298.0])  # exactly t_opt
        f = jarvis_single._f_temperature(leaf_T)
        assert np.isclose(f, 1.0)

    def test_f_temperature_cold_clip(self, jarvis_single: CanopyJarvis) -> None:
        # 50K below optimum drives the parabola negative → clipped to 0.
        leaf_T = np.array([248.0])
        f = jarvis_single._f_temperature(leaf_T)
        assert f[0] == 0.0

    def test_f_vpd_saturated_air(self, jarvis_single: CanopyJarvis) -> None:
        from utahlsm.physics import thermo
        T: NDArray[np.float64] = np.array([293.15])
        p: NDArray[np.float64] = np.array([101325.0])
        q_sat = thermo.saturation_specific_humidity(T, p)
        atm = AtmosphericState(
            temperature=T, pressure=p, specific_humidity=q_sat,
            wind_speed=np.array([3.0]),
            sw_in=np.array([400.0]),
        )
        f = jarvis_single._f_vpd(atm)
        assert np.isclose(f, 1.0)

    def test_f_vpd_dry_air_reduces(self, jarvis_single: CanopyJarvis) -> None:
        from utahlsm.physics import thermo
        T: NDArray[np.float64] = np.array([293.15])
        p: NDArray[np.float64] = np.array([101325.0])
        q_sat = thermo.saturation_specific_humidity(T, p)
        atm = AtmosphericState(
            temperature=T, pressure=p,
            specific_humidity=0.1 * q_sat,
            wind_speed=np.array([3.0]),
            sw_in=np.array([400.0]),
        )
        f = jarvis_single._f_vpd(atm)
        assert 0.0 < f[0] < 1.0

    def test_f_moisture_wilt_collapses(self, jarvis_single: CanopyJarvis) -> None:
        nz = 11
        theta = np.full(nz, 0.1)
        wilt = np.full(nz, 0.1)
        fc = np.full(nz, 0.3)
        f = jarvis_single._f_moisture(theta, wilt, fc)
        assert np.isclose(f, 0.0)

    def test_f_moisture_above_fc_ones(self, jarvis_single: CanopyJarvis) -> None:
        nz = 11
        theta = np.full(nz, 0.45)
        wilt = np.full(nz, 0.1)
        fc = np.full(nz, 0.3)
        f = jarvis_single._f_moisture(theta, wilt, fc)
        assert np.isclose(f, 1.0)


# ---------------------------------------------------------------------------
# Compute_resistance end-to-end
# ---------------------------------------------------------------------------

@pytest.mark.canopy
class TestComputeResistance:
    """Tests of the combined r_s = rs_min / (LAI · Π f_i) response."""

    @staticmethod
    def _states(nz: int) -> tuple[AtmosphericState, SurfaceState, SoilState]:
        atm = AtmosphericState(
            wind_speed=np.array([3.0]),
            temperature=np.array([298.0]),
            specific_humidity=np.array([0.01]),
            pressure=np.array([101325.0]),
            sw_in=np.array([500.0]),
            radiation_net=np.array([400.0]),
        )
        sfc = SurfaceState(temperature=np.array([298.0]))
        soil = SoilState(
            temperature=np.full(nz, 293.15),
            moisture=np.full(nz, 0.25),
            type=np.array(['clay'] * nz),
        )
        return atm, sfc, soil

    def test_r_s_bounded_by_rs_max(self, jarvis_single: CanopyJarvis) -> None:
        atm, sfc, soil = self._states(11)
        # Shut the radiation stress → total stress collapses.
        atm.sw_in = np.array([0.0])
        r_s = jarvis_single.compute_resistance(
            atm, sfc, soil,
            theta_wilt=np.full(11, 0.1),
            theta_fc=np.full(11, 0.3),
        )
        assert r_s.shape == (1,)
        assert np.isclose(r_s[0], 5000.0)

    def test_r_s_bounded_below_by_rs_min_over_lai(self, jarvis_single: CanopyJarvis) -> None:
        """At optimal conditions, r_s ≈ rs_min / LAI."""
        atm, sfc, soil = self._states(11)
        # High radiation, optimal temp, saturated air, wet soil.
        from utahlsm.physics import thermo
        atm.sw_in = np.array([5_000.0])  # saturates f1
        sfc.temperature = np.array([298.0])      # f3 = 1 (leaf ≡ sfc)
        q_sat = thermo.saturation_specific_humidity(
            atm.temperature, atm.pressure)
        atm.specific_humidity = q_sat  # f2 = 1
        soil.moisture = np.full(11, 0.45)  # above FC → f4 = 1
        r_s = jarvis_single.compute_resistance(
            atm, sfc, soil,
            theta_wilt=np.full(11, 0.1),
            theta_fc=np.full(11, 0.3),
        )
        # rs_min / LAI = 40 / 3 ≈ 13.3 (f1 < 1 slightly). Loose upper
        # check against rs_min / LAI * (1 / f1).
        assert r_s[0] >= 40.0 / 3.0
        assert r_s[0] < 40.0  # well below rs_min single-leaf

    def test_r_s_moisture_monotonic(self, jarvis_single: CanopyJarvis) -> None:
        atm, sfc, soil = self._states(11)
        wilt = np.full(11, 0.15)
        fc = np.full(11, 0.30)

        soil.moisture = np.full(11, 0.16)  # near wilt
        r_wet_near = jarvis_single.compute_resistance(
            atm, sfc, soil, wilt, fc)

        soil.moisture = np.full(11, 0.28)  # near FC
        r_wet_fc = jarvis_single.compute_resistance(
            atm, sfc, soil, wilt, fc)

        # Wetter soil → lower resistance.
        assert r_wet_fc[0] < r_wet_near[0]


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------

@pytest.mark.canopy
class TestFactory:
    """Tests ``get_canopy_model`` namelist dispatch and broadcasting."""

    def test_none_returns_none(self, z_layers: NDArray[np.float64]) -> None:
        cfg = CanopyConfig(model='none')
        canopy = get_canopy_model(cfg, z_layers, ncol=1)
        assert canopy is None

    def test_jarvis_scalar_broadcasts_to_ncol(self, z_layers: NDArray[np.float64]) -> None:
        cfg = CanopyConfig(
            model='jarvis', lai=3.0, veg_fraction=0.9,
            rooting_depth=0.4,
        )
        canopy = get_canopy_model(cfg, z_layers, ncol=4)
        assert isinstance(canopy, CanopyJarvis)
        assert canopy.lai.shape == (4,)
        assert np.allclose(canopy.lai, 3.0)
        assert canopy.veg_fraction.shape == (4,)

    def test_jarvis_sequence_respected(self, z_layers: NDArray[np.float64]) -> None:
        cfg = CanopyConfig(
            model='jarvis', lai=[1.0, 2.0, 3.0], veg_fraction=0.5,
            rooting_depth=0.4,
        )
        canopy = get_canopy_model(cfg, z_layers, ncol=3)
        assert isinstance(canopy, CanopyJarvis)
        assert np.allclose(canopy.lai, [1.0, 2.0, 3.0])

    def test_invalid_model_raises(self, z_layers: NDArray[np.float64]) -> None:
        cfg = CanopyConfig(model='penman')
        with pytest.raises(NamelistError):
            get_canopy_model(cfg, z_layers, ncol=1)

    def test_size_mismatch_raises(self, z_layers: NDArray[np.float64]) -> None:
        cfg = CanopyConfig(
            model='jarvis', lai=[1.0, 2.0], veg_fraction=0.5,
            rooting_depth=0.4,
        )
        with pytest.raises(NamelistError):
            get_canopy_model(cfg, z_layers, ncol=3)

    @pytest.mark.parametrize(
        ('overrides', 'match'),
        [
            ({'lai': -1.0}, 'lai'),
            ({'veg_fraction': 1.1}, 'veg_fraction'),
            ({'rooting_depth': -0.1}, 'rooting_depth'),
            ({'beta': 1.0}, 'beta'),
            ({'rs_min': 0.0}, 'rs_min'),
            ({'rs_min': 50.0, 'rs_max': 40.0}, 'rs_max'),
            ({'rg_half': 0.0}, 'rg_half'),
            ({'vpd_coef': -1.0e-4}, 'vpd_coef'),
            ({'t_coef': -1.0e-3}, 't_coef'),
        ],
    )
    def test_invalid_parameters_raise(
        self, z_layers: NDArray[np.float64], overrides: dict[str, Any], match: str
    ) -> None:
        cfg_kwargs: dict[str, Any] = {
            'model': 'jarvis',
            'lai': 3.0,
            'veg_fraction': 0.9,
            'rooting_depth': 0.4,
        }
        cfg_kwargs.update(overrides)
        cfg = CanopyConfig(**cfg_kwargs)
        with pytest.raises(NamelistError, match=match):
            get_canopy_model(cfg, z_layers, ncol=1)


# ---------------------------------------------------------------------------
# Canopy base-class is abstract
# ---------------------------------------------------------------------------

@pytest.mark.canopy
class TestCanopyABC:
    def test_cannot_instantiate_directly(self, canopy_params_single: dict[str, Any]) -> None:
        canopy_cls = cast(Any, Canopy)
        with pytest.raises(TypeError):
            canopy_cls(**canopy_params_single)


@pytest.mark.canopy
class TestCanopyParameterShapes:
    """Direct canopy construction still normalizes scalar column parameters."""

    def test_direct_constructor_broadcasts_scalars_to_inferred_ncol(
        self, canopy_params_3col: dict[str, Any]
    ) -> None:
        params = dict(canopy_params_3col)
        params['veg_fraction'] = 0.5
        params['r_ground'] = 3.0
        params['rg_half'] = 40.0

        canopy = CanopyJarvis(**params)

        assert canopy.ncol == 3
        assert canopy.veg_fraction.shape == (3,)
        assert np.allclose(canopy.veg_fraction, 0.5)
        assert canopy.r_ground.shape == (3,)
        assert np.allclose(canopy.r_ground, 3.0)
        assert canopy.rg_half.shape == (3,)
        assert np.allclose(canopy.rg_half, 40.0)

    def test_direct_constructor_rejects_mismatched_column_sizes(
        self, canopy_params_3col: dict[str, Any]
    ) -> None:
        params = dict(canopy_params_3col)
        params['veg_fraction'] = np.array([0.2, 0.8])

        with pytest.raises(ValueError, match='sizes disagree'):
            CanopyJarvis(**params)


@pytest.mark.canopy
def test_supersaturated_air_does_not_create_negative_root_uptake(
    jarvis_single: CanopyJarvis,
    z_layers: NDArray[np.float64],
) -> None:
    """Vegetation dew condensation must not be routed backward into roots."""
    model: Any = UtahLSM.__new__(UtahLSM)
    model.logger = logging.getLogger("test")
    model.ncol = 1
    model.canopy = jarvis_single
    model.tstep = 1800.0

    sfc_T = np.array([280.0])
    atm_p = np.array([101325.0])
    q_sat_sfc = thermo.saturation_specific_humidity(sfc_T, atm_p)
    atm_q = 1.2 * q_sat_sfc

    model.atm_state = AtmosphericState(
        wind_speed=np.array([3.0]),
        temperature=np.array([280.0]),
        specific_humidity=atm_q,
        pressure=atm_p,
        sw_in=np.array([0.0]),
        lw_in=np.array([280.0]),
        radiation_net=np.array([100.0]),
    )
    model.sfc_state = SurfaceState(
        temperature=sfc_T,
        moisture=np.array([0.25]),
        air_density=atm_p / (287.04 * np.array([280.0]) * (1.0 + 0.608 * atm_q)),
        turbulence=TurbulenceScales(
            friction_velocity=np.array([0.3]),
            obukhov_length=np.array([100.0]),
        ),
    )
    model.soil_state = SoilState(
        temperature=np.full((z_layers.size, 1), 280.0),
        moisture=np.full((z_layers.size, 1), 0.25),
        type=np.array(["clay"] * z_layers.size, dtype=object)[:, None],
    )
    model.canopy_state = CanopyState(
        resistance=np.array([100.0]),
        theta_root=np.array([0.25]),
        transpiration=np.zeros(1),
        wet_evaporation=np.zeros(1),
        evap_soil=np.zeros(1),
        water_storage=np.zeros(1),
        water_capacity=(
            jarvis_single.veg_fraction
            * jarvis_single.lai
            * jarvis_single.water_capacity_lai
        ),
        latent_veg=np.zeros(1),
        latent_wet=np.zeros(1),
        latent_soil=np.zeros(1),
        root_uptake=np.zeros((z_layers.size, 1)),
    )
    model.input = SimpleNamespace(
        surface=SimpleNamespace(z_s=2.0, z_t=0.01),
        grid=SimpleNamespace(nz=z_layers.size, z=z_layers),
    )
    model.sfc = SimpleNamespace(
        fh=lambda _z_s, _z_t, L: np.full_like(L, 0.1),
    )
    model.soil = SimpleNamespace(
        surface_specific_humidity=lambda _T, _q, _p: atm_q,
    )

    flux = model._partition_flux_wq(
        sfc_T, atm_q, atm_q, atm_p, np.array([0.3]), np.array([0.1])
    )
    model._finalize_canopy_partition()

    assert flux[0] < 0.0
    assert np.allclose(model.canopy_state.transpiration, 0.0)
    assert np.allclose(model.canopy_state.latent_veg, 0.0)
    assert model.canopy_state.latent_wet[0] < 0.0
    assert np.allclose(model.canopy_state.latent_soil, 0.0)
    assert model.canopy_state.water_storage[0] > 0.0
    assert np.all(model.canopy_state.root_uptake >= 0.0)
    assert np.allclose(model.canopy_state.root_uptake, 0.0)
