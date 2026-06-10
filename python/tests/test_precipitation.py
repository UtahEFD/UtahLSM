"""Tests for liquid precipitation forcing and saturation-excess runoff.

Exercises:
- The SMB residual with a P_infil term (sign convention, runoff branch).
- The NetCDF forcing loader: optional precip variable, validation,
  cold-temperature warning.
"""

from __future__ import annotations

import logging
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import netCDF4 as nc
import numpy as np
import pytest

from utahlsm.core import UtahLSM
from utahlsm.data_models import (
    AtmosphericState,
    CanopyState,
    IterationsConfig,
    NumericsConfig,
    SoilState,
    SurfaceFluxes,
    SurfaceState,
    TolerancesConfig,
    TurbulenceScales,
)
from utahlsm.physics.soil.soil_brookscorey import BrooksCorey
from utahlsm.util.io.input import Input
from utahlsm.util.io.soil_properties_loader import SoilPropertiesLoader


def _make_smb_model(
    *,
    theta_profile: float = 0.30,
    precipitation: float = 0.0,
    soil_type: str = "loam",
    tstep: float = 60.0,
    z: np.ndarray | None = None,
) -> Any:
    """Build a minimal UtahLSM with the dependencies needed by _solve_smb."""
    model: Any = UtahLSM.__new__(UtahLSM)
    model.logger = logging.getLogger("test.precip")
    model.ncol = 1
    model.tstep = tstep
    if z is None:
        z = np.array([0.0, -0.15, -0.30])

    model.input = SimpleNamespace(
        grid=SimpleNamespace(nx=1, ny=1, nz=3, z=z),
        surface=SimpleNamespace(z_s=2.0, z_t=0.01),
        numerics=NumericsConfig(
            heat_diffusion_back_weight=1.0,
            iterations=IterationsConfig(
                sfc_flux=10,
                seb_bracket=10,
                seb_root=10,
                smb_flux=100,
                moisture_picard=50,
                coupling=10,
            ),
            tolerances=TolerancesConfig(
                sfc_flux=1e-3,
                seb_root=1e-6,
                smb_flux=1e-10,
                moisture_picard=1e-8,
                moisture_bounds=1e-6,
                coupling_temp=1e-3,
                coupling_mois=1e-8,
            ),
        ),
    )

    props = SoilPropertiesLoader.load("cosby")
    model.soil = BrooksCorey(props, [soil_type] * 3, "cosby")

    model.soil_state = SoilState(
        temperature=np.array([[293.15], [293.15], [293.15]]),
        moisture=np.array([[theta_profile], [theta_profile], [theta_profile]]),
        type=np.array([soil_type] * 3, dtype=object),
    )

    fluxes = SurfaceFluxes(
        kinematic_heat=np.zeros(1),
        kinematic_moisture=np.zeros(1),
        sensible_heat=np.zeros(1),
        latent_heat=np.zeros(1),
        ground_heat=np.zeros(1),
        runoff=np.zeros(1),
    )
    turb = TurbulenceScales(
        friction_velocity=np.array([0.3]),
        obukhov_length=np.array([1e6]),
    )
    model.sfc_state = SurfaceState(
        temperature=np.array([293.15]),
        moisture=np.array([theta_profile]),
        air_density=np.array([1.205]),
        fluxes=fluxes,
        turbulence=turb,
    )

    model.atm_state = AtmosphericState(
        wind_speed=np.array([3.0]),
        temperature=np.array([290.0]),
        specific_humidity=np.array([0.006]),
        pressure=np.array([101325.0]),
        precipitation=np.array([precipitation]),
    )

    # Neutral stability stub: fh constant. f_veg=0 (no canopy).
    model.sfc = SimpleNamespace(
        fh=lambda _z_s, _z_t, L: np.full_like(L, 0.1),
    )
    model.canopy = None

    return model


def _infiltration_capacity(model: Any) -> float:
    """Green-Ampt infiltration capacity [m/s], mirroring ``_solve_smb``.

    The gravity-limited saturated conductivity at the wet bracket
    endpoint, enhanced by the wetting-front suction term
    ``1 + delta_theta * psi_f / F`` with F the cumulative-infiltration
    event state (``model._ga_cum_infil``, floored). Scaled by rho_w it is
    the maximum mass infiltration rate before infiltration-excess runoff
    is shed.
    """
    residual_q = float(model.soil.properties.residual[0])
    porosity = float(model.soil.properties.porosity[0])
    span = porosity - residual_q
    theta_wet = np.full_like(model.sfc_state.temperature, porosity - 1e-6 * span)
    K0_sat = np.asarray(
        model.soil.conductivity_moisture(theta_wet, level=0), dtype=float
    )
    theta_ante = np.asarray(model.soil_state.moisture, dtype=float)[0]
    F = np.maximum(
        np.asarray(getattr(model, '_ga_cum_infil', 0.0), dtype=float),
        UtahLSM._GA_MIN_CUM_INFIL,
    )
    psi_f = np.asarray(model._wetting_front_suction(theta_ante), dtype=float)
    delta_theta = np.maximum(porosity - theta_ante, 0.0)
    capacity = K0_sat * (1.0 + delta_theta * psi_f / F)
    return float(capacity[0])


def _add_simple_canopy(
    model: Any,
    *,
    veg_fraction: float = 0.5,
    storage: float = 0.0,
    capacity: float = 0.1,
) -> None:
    """Attach the minimal canopy attributes used by precipitation routing."""
    model.canopy = SimpleNamespace(
        veg_fraction=np.array([veg_fraction]),
        wet_cooling_max=np.array([0.0]),
    )
    model.canopy_state = CanopyState(
        water_storage=np.array([storage], dtype=float),
        water_capacity=np.array([capacity], dtype=float),
    )
    model._soil_precipitation = np.zeros(1)
    model._precipitation_prepared = False


# ---------------------------------------------------------------------------
# SMB physics
# ---------------------------------------------------------------------------


@pytest.mark.precip
@pytest.mark.smb
def test_zero_precip_baseline_runoff_is_zero() -> None:
    """With P=0, SMB behaves as before: runoff stays zero, theta_sfc in bounds."""
    model = _make_smb_model(theta_profile=0.30, precipitation=0.0)
    porosity = float(model.soil.properties.porosity[0])
    residual_q = float(model.soil.properties.residual[0])

    model._solve_smb()

    assert np.all(model.sfc_state.fluxes.runoff == 0.0)
    assert np.all(model.sfc_state.moisture >= residual_q)
    assert np.all(model.sfc_state.moisture <= porosity)


@pytest.mark.precip
@pytest.mark.smb
def test_light_rain_no_runoff() -> None:
    """A light rain rate produces no runoff (soil absorbs it all)."""
    model = _make_smb_model(theta_profile=0.30, precipitation=1e-5)
    porosity = float(model.soil.properties.porosity[0])

    model._solve_smb()

    assert np.all(model.sfc_state.fluxes.runoff == 0.0)
    assert float(model.sfc_state.moisture[0]) <= porosity


@pytest.mark.precip
@pytest.mark.canopy
def test_canopy_intercepts_vegetated_rain_before_smb() -> None:
    """Rain falling over vegetation fills canopy storage before throughfall."""
    model = _make_smb_model(precipitation=1.0e-3, tstep=60.0)
    _add_simple_canopy(model, veg_fraction=0.5, storage=0.0, capacity=0.1)

    model._prepare_precipitation_for_timestep()

    assert model._precipitation_prepared
    assert model.canopy_state.water_storage[0] == pytest.approx(0.03)
    assert model._soil_precipitation[0] == pytest.approx(5.0e-4)


@pytest.mark.precip
@pytest.mark.canopy
def test_canopy_interception_limited_by_remaining_storage() -> None:
    """Once the canopy bucket fills, excess rain becomes throughfall."""
    model = _make_smb_model(precipitation=1.0e-3, tstep=60.0)
    _add_simple_canopy(model, veg_fraction=0.95, storage=0.0, capacity=0.01)

    model._prepare_precipitation_for_timestep()

    intercepted_rate = 0.01 / 60.0
    assert model.canopy_state.water_storage[0] == pytest.approx(0.01)
    assert model._soil_precipitation[0] == pytest.approx(1.0e-3 - intercepted_rate)


@pytest.mark.precip
@pytest.mark.canopy
def test_full_canopy_storage_passes_all_rain_to_soil() -> None:
    """A saturated canopy does not remove additional rain from throughfall."""
    model = _make_smb_model(precipitation=1.0e-3, tstep=60.0)
    _add_simple_canopy(model, veg_fraction=0.95, storage=0.1, capacity=0.1)

    model._prepare_precipitation_for_timestep()

    assert model.canopy_state.water_storage[0] == pytest.approx(0.1)
    assert model._soil_precipitation[0] == pytest.approx(1.0e-3)


@pytest.mark.precip
@pytest.mark.canopy
def test_intercepted_canopy_water_suppresses_dry_transpiration() -> None:
    """Rain-filled canopy storage reduces the dry transpiring fraction."""
    dry = _make_smb_model()
    wet = _make_smb_model(precipitation=1.0e-3, tstep=60.0)
    _add_simple_canopy(dry, veg_fraction=1.0, storage=0.0, capacity=0.06)
    _add_simple_canopy(wet, veg_fraction=1.0, storage=0.0, capacity=0.06)
    dry.canopy_state.resistance = np.array([0.0])
    wet.canopy_state.resistance = np.array([0.0])
    wet._prepare_precipitation_for_timestep()

    sfc_T = np.array([300.0])
    atm_p = np.array([101325.0])
    atm_q = np.array([0.005])
    gnd_q = np.array([0.010])
    ust = np.array([0.3])
    fh = np.array([0.1])

    _, dry_t, _ = dry._partition_flux_wq_components(
        sfc_T, gnd_q, atm_q, atm_p, ust, fh
    )
    _, wet_t, _ = wet._partition_flux_wq_components(
        sfc_T, gnd_q, atm_q, atm_p, ust, fh
    )

    assert wet.canopy_state.water_storage[0] == pytest.approx(0.06)
    assert wet_t[0] == pytest.approx(0.0)
    assert dry_t[0] > wet_t[0]


@pytest.mark.precip
@pytest.mark.smb
def test_heavy_rain_produces_runoff() -> None:
    """Rain far above the late-event infiltration capacity runs off.

    With a large cumulative-infiltration state the Green-Ampt suction
    enhancement has decayed, so only roughly the gravity-limited capacity
    rho_w*K_sat enters the soil; the remainder is infiltration-excess
    (Hortonian) runoff.
    """
    RHO_W = 1000.0
    model = _make_smb_model(theta_profile=0.30)
    model._ga_cum_infil = np.array([1.0])  # deep into the event: sealed
    infil_cap = RHO_W * _infiltration_capacity(model)

    P_heavy = infil_cap * 20.0
    model.atm_state.precipitation = np.array([P_heavy])
    model._solve_smb()

    runoff = float(model.sfc_state.fluxes.runoff[0])
    infiltrated = P_heavy - runoff
    assert infiltrated == pytest.approx(infil_cap, rel=1e-6)
    # Sanity: runoff cannot exceed total precipitation.
    assert 0.0 < runoff <= P_heavy


@pytest.mark.precip
@pytest.mark.smb
def test_rain_above_late_event_capacity_runs_off() -> None:
    """Rain above the sealed-surface capacity runs off late in an event.

    Once sustained infiltration has exhausted the wetting-front suction
    enhancement, a rate exceeding the (near-K_sat) capacity is shed as
    infiltration-excess runoff rather than buffered as storage.
    """
    RHO_W = 1000.0
    model = _make_smb_model(theta_profile=0.10)
    model._ga_cum_infil = np.array([1.0])
    infil_cap = RHO_W * _infiltration_capacity(model)
    P_intense = infil_cap * 10.0

    model.atm_state.precipitation = np.array([P_intense])
    model._solve_smb()

    runoff = float(model.sfc_state.fluxes.runoff[0])
    assert runoff > 0.0
    assert (P_intense - runoff) == pytest.approx(infil_cap, rel=1e-6)


@pytest.mark.precip
@pytest.mark.smb
def test_dry_soil_fresh_event_absorbs_intense_rain() -> None:
    """Early-storm rain on dry soil infiltrates without Hortonian runoff.

    At the start of an event (cumulative infiltration at its floor) the
    Green-Ampt suction term dominates on dry soil, so a burst well above
    rho_w*K_sat is absorbed - matching observed behaviour where dry
    fields swallow intense convective rain.
    """
    RHO_W = 1000.0
    model = _make_smb_model(theta_profile=0.10)
    model._ga_cum_infil = np.zeros(1)
    K_sat_top = float(model.soil.properties.K_sat[0])

    model.atm_state.precipitation = np.array([10.0 * RHO_W * K_sat_top])
    model._solve_smb()

    assert float(model.sfc_state.fluxes.runoff[0]) == 0.0


@pytest.mark.precip
@pytest.mark.smb
def test_capacity_decays_with_cumulative_infiltration() -> None:
    """The Green-Ampt capacity decreases monotonically toward K_sat."""
    model = _make_smb_model(theta_profile=0.10)
    K_sat_top = float(model.soil.properties.K_sat[0])

    model._ga_cum_infil = np.array([1.0e-4])
    cap_early = _infiltration_capacity(model)
    model._ga_cum_infil = np.array([0.1])
    cap_late = _infiltration_capacity(model)

    assert cap_early > cap_late > K_sat_top


@pytest.mark.precip
@pytest.mark.smb
def test_infiltration_history_accumulates_and_decays() -> None:
    """The event state grows while raining and relaxes after rain stops."""
    RHO_W = 1000.0
    model = _make_smb_model(theta_profile=0.30, tstep=60.0)
    model._ga_cum_infil = np.array([0.01])
    model._infiltration_flux = np.array([2.0e-3])

    model.atm_state.precipitation = np.array([2.0e-3])
    model._update_infiltration_history()
    expected = 0.01 + 2.0e-3 * 60.0 / RHO_W
    assert model._ga_cum_infil[0] == pytest.approx(expected)

    model.atm_state.precipitation = np.array([0.0])
    model._update_infiltration_history()
    decayed = expected * np.exp(-60.0 / UtahLSM._GA_REDIST_TAU)
    assert model._ga_cum_infil[0] == pytest.approx(decayed)


# ---------------------------------------------------------------------------
# Macropore bypass flow
# ---------------------------------------------------------------------------


def _enable_macropores(
    model: Any,
    *,
    fraction: float = 0.5,
    z_top: float = 0.10,
    z_bottom: float = 0.20,
    e_folding: float = 0.2,
) -> None:
    """Attach a macropore soil config to a partial SMB test rig."""
    model.input.soil = SimpleNamespace(
        macropore_fraction=fraction,
        macropore_z_top=z_top,
        macropore_z_bottom=z_bottom,
        macropore_e_folding=e_folding,
    )


@pytest.mark.precip
@pytest.mark.smb
def test_bypass_capture_diverts_throughfall_fraction() -> None:
    """The macropore fraction of rain skips the matrix surface balance."""
    P_rain = 1.0e-3
    model = _make_smb_model(theta_profile=0.30, precipitation=P_rain)
    _enable_macropores(model, fraction=0.4)

    model._solve_smb()

    assert float(model._bypass_flux[0]) == pytest.approx(0.4 * P_rain)
    # The matrix only ever sees the non-bypass remainder.
    assert float(model._infiltration_flux[0]) <= 0.6 * P_rain + 1e-12


@pytest.mark.precip
@pytest.mark.smb
def test_bypass_disabled_by_default() -> None:
    """Without macropore config, no rain is diverted."""
    P_rain = 1.0e-3
    model = _make_smb_model(theta_profile=0.30, precipitation=P_rain)

    model._solve_smb()

    assert float(np.asarray(model._bypass_flux)[0]) == 0.0


@pytest.mark.precip
@pytest.mark.soil
def test_bypass_distribution_deposits_in_zone() -> None:
    """Captured water is deposited only in the configured crack zone."""
    model = _make_smb_model(theta_profile=0.10)  # below field capacity
    _enable_macropores(model, fraction=0.5, z_top=0.10, z_bottom=0.20)
    M = 1.0e-4
    model._bypass_flux = np.array([M])

    model._distribute_bypass()

    deposit = model._bypass_deposit
    assert deposit is not None
    assert deposit.shape == (3, 1)
    assert deposit[0, 0] == 0.0          # surface BC layer excluded
    assert deposit[1, 0] == pytest.approx(M)   # node at 0.15 m, in zone
    assert deposit[2, 0] == 0.0          # node at 0.30 m, below zone
    assert float(model.sfc_state.fluxes.runoff[0]) == 0.0


@pytest.mark.precip
@pytest.mark.soil
def test_bypass_overflow_goes_to_runoff() -> None:
    """Bypass water the zone cannot absorb ponds and runs off."""
    model = _make_smb_model(theta_profile=0.10)  # below field capacity
    _enable_macropores(model, fraction=0.5, z_top=0.10, z_bottom=0.20)

    # Zone headroom: fill node 1 from theta to field capacity this step.
    RHO_W = 1000.0
    dz = 0.15
    theta_fc1 = float(np.asarray(model.soil.theta_fc)[1])
    headroom = (theta_fc1 - 0.10) * dz * RHO_W / model.tstep
    M = 3.0 * headroom
    model._bypass_flux = np.array([M])

    model._distribute_bypass()

    deposit = model._bypass_deposit
    assert deposit is not None
    assert deposit[1, 0] == pytest.approx(headroom, rel=1e-6)
    assert float(model.sfc_state.fluxes.runoff[0]) == pytest.approx(
        M - headroom, rel=1e-6
    )
    # The reported bypass flux is what was actually deposited.
    assert float(model._bypass_flux[0]) == pytest.approx(headroom, rel=1e-6)


@pytest.mark.precip
@pytest.mark.soil
def test_no_bypass_flux_no_deposit() -> None:
    """With zero captured flux the deposition profile is omitted."""
    model = _make_smb_model(theta_profile=0.20)
    _enable_macropores(model)
    model._bypass_flux = np.zeros(1)

    model._distribute_bypass()

    assert model._bypass_deposit is None


@pytest.mark.precip
@pytest.mark.soil
def test_bypass_deposit_cools_deep_layer() -> None:
    """Cold bypass water deposited at depth cools that layer directly.

    With no matrix infiltration and a uniform warm column, the only
    heating term is the bypass deposit at rain temperature in layer 2 -
    a reach the surface-input term alone cannot produce.
    """
    model = _make_smb_model(theta_profile=0.30)
    model.soil_state.temperature[:] = 303.0
    model.atm_state.temperature = np.array([290.0])
    model._infiltration_flux = np.zeros(1)
    deposit = np.zeros((3, 1))
    deposit[2, 0] = 1.0e-4
    model._bypass_deposit = deposit

    source = model._rain_advection_source()

    assert source is not None
    assert source[0, 0] == 0.0
    assert source[1, 0] == pytest.approx(0.0)
    assert source[2, 0] < 0.0

    c_w = 4.184e6 / 1000.0
    dz = 0.15
    c_vol2 = float(
        np.asarray(model.soil.heat_capacity(model.soil_state.moisture))[2, 0]
    )
    expected = c_w * 1.0e-4 * (290.0 - 303.0) / (c_vol2 * dz)
    assert source[2, 0] == pytest.approx(expected, rel=1e-6)


@pytest.mark.precip
@pytest.mark.smb
def test_smb_storage_tendency_limits_short_step_wetting() -> None:
    """Top-layer moisture changes scale with dt through the storage term."""
    z = np.array([0.0, -0.01, -0.02])
    # Keep the rate below the saturated infiltration capacity so no
    # infiltration-excess runoff is shed; this isolates the dt-dependent
    # storage tendency. cosby loam K_sat ~ 5.1e-6 m/s -> cap ~ 5.1e-3.
    P_rain = 2.0e-3
    theta_initial = 0.10

    short = _make_smb_model(
        theta_profile=theta_initial,
        precipitation=P_rain,
        tstep=1.0,
        z=z,
    )
    long = _make_smb_model(
        theta_profile=theta_initial,
        precipitation=P_rain,
        tstep=3600.0,
        z=z,
    )

    short._solve_smb()
    long._solve_smb()

    theta_short = float(short.sfc_state.moisture[0])
    theta_long = float(long.sfc_state.moisture[0])
    assert float(short.sfc_state.fluxes.runoff[0]) == 0.0
    assert float(long.sfc_state.fluxes.runoff[0]) == 0.0
    assert theta_initial < theta_short < theta_long


@pytest.mark.precip
@pytest.mark.smb
def test_rain_below_capacity_no_runoff() -> None:
    """Rain below rho_w*K_sat still produces zero runoff."""
    probe = _make_smb_model()
    K_sat_top = float(probe.soil.properties.K_sat[0])
    RHO_W = 1000.0
    P_modest = 0.5 * RHO_W * K_sat_top  # half of capacity

    model = _make_smb_model(theta_profile=0.30, precipitation=P_modest)
    model._solve_smb()

    assert float(model.sfc_state.fluxes.runoff[0]) == 0.0


@pytest.mark.precip
@pytest.mark.smb
def test_precip_field_propagates_via_load_atm_state() -> None:
    """update() -> _load_atm_state copies precipitation onto self.atm_state."""
    model = _make_smb_model(precipitation=0.0)
    new_atm = AtmosphericState(
        wind_speed=np.array([3.0]),
        temperature=np.array([290.0]),
        specific_humidity=np.array([0.006]),
        pressure=np.array([101325.0]),
        sw_in=np.array([0.0]),
        lw_in=np.array([300.0]),
        seb_storage=np.array([0.0]),
        precipitation=np.array([2.5e-4]),
    )

    model._load_atm_state(new_atm)

    assert np.allclose(model.atm_state.precipitation, 2.5e-4)


# ---------------------------------------------------------------------------
# Rain heat advection
# ---------------------------------------------------------------------------


@pytest.mark.precip
@pytest.mark.soil
def test_rain_advection_cools_warm_soil() -> None:
    """Cold infiltrating rain deposits a cooling tendency in the top layer.

    The source goes only into layer 1 (layer 0 is the Dirichlet skin BC),
    and its magnitude is c_w*P*(T_rain-T)/(C_vol*dz).
    """
    model = _make_smb_model(theta_profile=0.30)
    model.soil_state.temperature[:] = 303.0
    model.atm_state.temperature = np.array([290.0])
    model._infiltration_flux = np.array([1.0e-3])

    source = model._rain_advection_source()

    assert source is not None
    assert source.shape == (3, 1)
    assert source[0, 0] == 0.0          # skin BC layer untouched
    assert source[2, 0] == 0.0          # below the deposit layer untouched
    assert source[1, 0] < 0.0           # cold rain cools the top layer

    c_w = 4.184e6 / 1000.0
    dz = 0.15
    c_vol = float(
        np.asarray(model.soil.heat_capacity(model.soil_state.moisture))[1, 0]
    )
    expected = c_w * 1.0e-3 * (290.0 - 303.0) / (c_vol * dz)
    assert source[1, 0] == pytest.approx(expected, rel=1e-6)


@pytest.mark.precip
@pytest.mark.soil
def test_rain_advection_sign_tracks_temperature_difference() -> None:
    """Rain warmer than the soil warms the top layer; equal temps do nothing."""
    model = _make_smb_model(theta_profile=0.30)
    model.soil_state.temperature[:] = 290.0
    model._infiltration_flux = np.array([1.0e-3])

    model.atm_state.temperature = np.array([295.0])
    warm = model._rain_advection_source()
    assert warm is not None and warm[1, 0] > 0.0

    model.atm_state.temperature = np.array([290.0])
    iso = model._rain_advection_source()
    assert iso is not None and iso[1, 0] == pytest.approx(0.0)


@pytest.mark.precip
@pytest.mark.soil
def test_no_infiltration_no_rain_source() -> None:
    """With no infiltration the rain heat source is omitted entirely."""
    model = _make_smb_model(theta_profile=0.30)
    model._infiltration_flux = np.zeros(1)
    assert model._rain_advection_source() is None


@pytest.mark.precip
@pytest.mark.soil
def test_interior_advection_cools_layer_below_top() -> None:
    """Term 2: the downward Darcy flux advects cold into a deeper layer.

    With a cold top layer over warm subsoil, the gravity-driven downward
    water flux carries the top layer's deficit into the layer below - a
    reach the old surface-input-only source (which wrote layer 1 only)
    could not produce.
    """
    model = _make_smb_model(theta_profile=0.30)
    # Top prognostic layer cold, layer below warm. T_rain == T_1 isolates
    # the interior term by zeroing layer 1's surface-input contribution.
    model.soil_state.temperature[:] = np.array([[300.0], [290.0], [303.0]])
    model.atm_state.temperature = np.array([290.0])
    model._infiltration_flux = np.array([1.0e-3])

    source = model._rain_advection_source()

    assert source is not None
    assert source[0, 0] == 0.0          # skin BC layer untouched
    assert source[1, 0] == pytest.approx(0.0)   # T_rain == T_1: no input
    assert source[2, 0] < 0.0           # cold water advected downward

    # Exact value: c_w * m_face1 * (T_1 - T_2) / (C_vol_2 * dz), with the
    # interior face mass flux from the start-of-step Darcy state.
    c_w = 4.184e6 / 1000.0
    rho_w = 1000.0
    dz = 0.15
    theta = model.soil_state.moisture
    psi = np.asarray(model.soil.water_potential(theta))
    k_node = np.asarray(model.soil.conductivity_moisture(theta))
    k_face = model._moisture_face_conductivity(
        k_node[:-1], k_node[1:], psi[:-1], psi[1:], dz
    )
    m_face1 = rho_w * float(k_face[1, 0]) * (1.0 - (psi[2, 0] - psi[1, 0]) / dz)
    c_vol2 = float(np.asarray(model.soil.heat_capacity(theta))[2, 0])
    expected = c_w * m_face1 * (290.0 - 303.0) / (c_vol2 * dz)
    assert source[2, 0] == pytest.approx(expected, rel=1e-6)


# ---------------------------------------------------------------------------
# NetCDF forcing loader
# ---------------------------------------------------------------------------


def _make_input_loader() -> Input:
    """Build a bare Input instance with just the attrs the loader needs."""
    obj = Input.__new__(Input)
    obj.logger = logging.getLogger("tests.precip.input")
    obj.grid = SimpleNamespace(nx=1, ny=1)
    return obj


def _write_forcing_nc(
    path: Path,
    *,
    atm_T: np.ndarray,
    precip: np.ndarray | None,
    ntime: int = 3,
) -> None:
    """Write a minimal offline forcing NetCDF for loader tests."""
    with nc.Dataset(path, "w") as ds:
        ds.createDimension("t", ntime)
        ds.createDimension("scalar", 1)
        ds.createVariable("tstep", "f8", ("scalar",))[:] = [60.0]
        for name, vals in (
            ("atm_U", np.full(ntime, 3.0)),
            ("atm_T", atm_T),
            ("atm_q", np.full(ntime, 0.006)),
            ("atm_p", np.full(ntime, 101325.0)),
            ("sw_in", np.full(ntime, 300.0)),
            ("lw_in", np.full(ntime, 350.0)),
        ):
            ds.createVariable(name, "f8", ("t",))[:] = vals
        if precip is not None:
            ds.createVariable("precip", "f8", ("t",))[:] = precip


@pytest.mark.precip
@pytest.mark.unit
def test_precip_optional_in_netcdf(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """A forcing file without `precip` loads cleanly with zero defaults."""
    path = tmp_path / "lsm_offline.nc"
    _write_forcing_nc(path, atm_T=np.full(3, 290.0), precip=None)

    loader = _make_input_loader()
    with caplog.at_level(logging.INFO, logger=loader.logger.name):
        loader._load_offline_data(str(path))

    assert loader.forcing is not None
    assert loader.forcing.ntime == 3
    for atm in loader.forcing.atmos:
        assert np.allclose(np.asarray(atm.precipitation), 0.0)
    assert any("precip not in forcing" in rec.getMessage() for rec in caplog.records)


@pytest.mark.precip
@pytest.mark.unit
def test_precip_validation_rejects_negative(tmp_path: Path) -> None:
    """Negative precipitation in the NetCDF must raise."""
    path = tmp_path / "lsm_offline.nc"
    precip = np.array([0.0, -1.0, 0.0])
    _write_forcing_nc(path, atm_T=np.full(3, 290.0), precip=precip)

    loader = _make_input_loader()
    with pytest.raises(ValueError, match="precip"):
        loader._load_offline_data(str(path))


@pytest.mark.precip
@pytest.mark.unit
def test_cold_temp_precip_warning(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """T < 273.15 K with P > 0 fires a single load-time WARNING."""
    path = tmp_path / "lsm_offline.nc"
    atm_T = np.array([270.0, 268.0, 285.0])
    precip = np.array([1e-5, 2e-5, 1e-5])
    _write_forcing_nc(path, atm_T=atm_T, precip=precip)

    loader = _make_input_loader()
    with caplog.at_level(logging.WARNING, logger=loader.logger.name):
        loader._load_offline_data(str(path))

    warnings = [r for r in caplog.records
                if r.levelno == logging.WARNING
                and "Cold-temperature precipitation" in r.getMessage()]
    assert len(warnings) == 1
    assert "2 sample" in warnings[0].getMessage()
