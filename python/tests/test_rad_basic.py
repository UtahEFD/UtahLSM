"""Tests for the basic clear-sky radiation parameterization.

``RadBasic`` is not exercised by any bundled case (both set
``radiation.model = "forcing"``), so its parameterizations need direct
coverage. These tests pin the Brutsaert (1975) unit convention, which is
easy to get wrong: the 1.24 coefficient is calibrated for vapor pressure
in hPa, and feeding it Pa silently produces an emissivity above unity.
"""

import numpy as np
import pytest

from utahlsm.data_models import AtmosphericState
from utahlsm.physics.radiation.rad_basic import RadBasic
from utahlsm.util import constants as c

SB = c.radiation.STEFAN_BOLTZMANN


@pytest.fixture
def rad() -> RadBasic:
    """A RadBasic instance with representative surface optics."""
    return RadBasic(latitude=36.6049, longitude=-97.4856,
                    albedo=0.19, emissivity=0.96)


def _state(T: float, q: float, p: float = 100000.0) -> AtmosphericState:
    """Builds a single-column atmospheric state."""
    return AtmosphericState(
        temperature=np.array([T]),
        specific_humidity=np.array([q]),
        pressure=np.array([p]),
    )


class TestBrutsaertLongwave:
    """Clear-sky downwelling longwave via Brutsaert (1975)."""

    def test_matches_closed_form_with_vapor_pressure_in_hpa(
        self, rad: RadBasic
    ) -> None:
        """lw_in reproduces 1.24*(e_hPa/T)^(1/7)*sigma*T^4."""
        T, q, p = 290.0, 0.010, 100000.0
        e_pa = p * q / (c.thermodynamic.EPSILON + q)
        expected_eps = 1.24 * (e_pa / 100.0 / T) ** (1.0 / 7.0)

        lw_in = float(np.atleast_1d(rad._longwave_in(_state(T, q, p)))[0])

        assert lw_in == pytest.approx(expected_eps * SB * T**4, rel=1e-12)

    def test_effective_emissivity_stays_below_unity(
        self, rad: RadBasic
    ) -> None:
        """A physical emissivity never exceeds 1 across the forcing envelope."""
        for T in (250.0, 273.15, 290.0, 315.0):
            for q in (1e-4, 1e-3, 5e-3, 0.010, 0.025):
                lw_in = float(
                    np.atleast_1d(rad._longwave_in(_state(T, q)))[0]
                )
                eps_eff = lw_in / (SB * T**4)
                assert 0.0 < eps_eff < 1.0, (
                    f"eps_eff={eps_eff:.3f} outside (0, 1) at T={T}, q={q}"
                )

    def test_never_exceeds_blackbody_emission(self, rad: RadBasic) -> None:
        """Downwelling longwave stays below sigma*T_air^4."""
        for T in (250.0, 290.0, 315.0):
            for q in (1e-3, 0.010, 0.025):
                lw_in = float(
                    np.atleast_1d(rad._longwave_in(_state(T, q)))[0]
                )
                assert lw_in < SB * T**4

    def test_increases_with_moisture_and_temperature(
        self, rad: RadBasic
    ) -> None:
        """lw_in rises with both vapor pressure and air temperature."""
        dry = float(np.atleast_1d(rad._longwave_in(_state(290.0, 2e-3)))[0])
        wet = float(np.atleast_1d(rad._longwave_in(_state(290.0, 0.015)))[0])
        assert wet > dry

        cool = float(np.atleast_1d(rad._longwave_in(_state(275.0, 0.008)))[0])
        warm = float(np.atleast_1d(rad._longwave_in(_state(305.0, 0.008)))[0])
        assert warm > cool

    def test_magnitude_is_physically_reasonable(self, rad: RadBasic) -> None:
        """A humid 290 K atmosphere emits roughly 300-360 W/m^2 downward."""
        lw_in = float(np.atleast_1d(rad._longwave_in(_state(290.0, 0.010)))[0])
        assert 300.0 < lw_in < 360.0

    def test_vectorizes_over_columns(self, rad: RadBasic) -> None:
        """Multi-column input returns per-column values, not a scalar."""
        state = AtmosphericState(
            temperature=np.array([280.0, 290.0, 300.0]),
            specific_humidity=np.array([0.004, 0.010, 0.018]),
            pressure=np.full(3, 100000.0),
        )
        lw_in = np.asarray(rad._longwave_in(state))

        assert lw_in.shape == (3,)
        assert np.all(np.diff(lw_in) > 0.0)


def _peak_sw_hour(rad: RadBasic, julian_day: int) -> float:
    """Returns the UTC hour at which modelled SW_in peaks, to 30 s."""
    hours = np.arange(0.0, 24.0, 1.0 / 120.0)
    sw = np.array(
        [
            float(np.atleast_1d(rad._shortwave_in(julian_day, h * 3600.0))[0])
            for h in hours
        ]
    )
    return float(hours[int(sw.argmax())])


class TestSolarGeometry:
    """Clear-sky shortwave, focused on the east-positive longitude convention."""

    # (name, latitude, longitude east-positive, julian day)
    SITES = [
        ("Cabauw", 51.9711, 4.9267, 183),      # gabls3, eastern hemisphere
        ("ARM SGP", 36.6049, -97.4856, 168),   # arm, western hemisphere
        ("Greenwich", 51.4800, 0.0000, 172),   # prime meridian control
        ("Sydney", -33.8688, 151.2093, 15),    # southern + far east
    ]

    @pytest.mark.parametrize("name,lat,lon,jday", SITES)
    def test_peak_matches_true_solar_noon(
        self, name: str, lat: float, lon: float, jday: int
    ) -> None:
        """Peak SW_in occurs at 12 - lon/15 UTC, wrapped to [0, 24).

        This pins the sign of the longitude term. Getting it backwards
        displaces the diurnal cycle by 2*(lon/15) hours, which is silent
        at Greenwich and 13 h wrong at ARM.
        """
        rad = RadBasic(latitude=lat, longitude=lon,
                       albedo=0.2, emissivity=0.97)
        expected = (12.0 - lon / 15.0) % 24.0

        peak = _peak_sw_hour(rad, jday)

        # Circular difference, tolerant of the 30 s sampling grid.
        error = (peak - expected + 12.0) % 24.0 - 12.0
        assert abs(error) < 0.05, (
            f"{name}: peak at {peak:.2f} h UTC, expected {expected:.2f} h "
            f"(error {error:+.2f} h)"
        )

    def test_sun_is_up_at_local_noon_and_down_at_local_midnight(self) -> None:
        """ARM sees full sun at 18:30 UTC and none at 05:30 UTC in June."""
        rad = RadBasic(latitude=36.6049, longitude=-97.4856,
                       albedo=0.19, emissivity=0.96)

        noon = float(np.atleast_1d(rad._shortwave_in(168, 18.5 * 3600))[0])
        midnight = float(np.atleast_1d(rad._shortwave_in(168, 5.5 * 3600))[0])

        assert noon > 800.0
        assert midnight == 0.0

    def test_never_negative_and_zero_when_sun_is_down(self) -> None:
        """SW_in is clipped to zero below the horizon, never negative."""
        rad = RadBasic(latitude=36.6049, longitude=-97.4856,
                       albedo=0.19, emissivity=0.96)
        sw = np.array(
            [
                float(np.atleast_1d(rad._shortwave_in(168, h * 3600.0))[0])
                for h in np.arange(0.0, 24.0, 0.05)
            ]
        )

        assert np.all(sw >= 0.0)
        assert np.any(sw == 0.0)
        assert sw.max() < c.radiation.SOLAR_CONSTANT

    def test_summer_solstice_peaks_higher_than_winter(self) -> None:
        """Northern-hemisphere peak insolation is larger near day 173."""
        rad = RadBasic(latitude=51.9711, longitude=4.9267,
                       albedo=0.2, emissivity=0.97)

        summer = max(
            float(np.atleast_1d(rad._shortwave_in(173, h * 3600.0))[0])
            for h in np.arange(0.0, 24.0, 0.05)
        )
        winter = max(
            float(np.atleast_1d(rad._shortwave_in(356, h * 3600.0))[0])
            for h in np.arange(0.0, 24.0, 0.05)
        )

        assert summer > winter
