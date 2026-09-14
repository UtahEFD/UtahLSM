"""Tests for the soil-heat theta-scheme monotonicity floor.

The theta scheme's amplification factor for a Fourier mode is
``G = (1 - theta_f*a*s) / (1 + theta_b*a*s)`` with ``s`` reaching 4 at the
grid-scale sawtooth. Any ``theta_b >= 0.5`` keeps ``|G| <= 1``, but ``G``
turns negative once ``theta_f*4*alpha > 1`` and the mode then decays by
flipping sign every step. Crank-Nicolson is non-oscillatory only up to
``alpha = 0.5``; a 600 s step on a 1 cm soil grid runs at ``alpha ~ 3.8``.
"""

import logging
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest

from utahlsm.core import UtahLSM
from utahlsm.data_models import (
    IterationsConfig,
    NumericsConfig,
    TolerancesConfig,
)

CASES_ROOT = (Path(__file__).resolve().parents[2] / "cases")


def _load_case(name: str) -> Any:
    """Loads a bundled case, skipping when its files are not present."""
    import utahlsm

    case = CASES_ROOT / name
    if not case.exists():
        pytest.skip(f"{name} case files not available.")
    return utahlsm.Input(
        str(case / "lsm_namelist.json"),
        str(case / "lsm_init.nc"),
        str(case / "lsm_offline.nc"),
    )


def _numerics(back_weight: float, monotone: bool) -> NumericsConfig:
    """Builds a NumericsConfig with the heat-diffusion knobs under test."""
    return NumericsConfig(
        heat_diffusion_back_weight=back_weight,
        heat_diffusion_monotone=monotone,
        iterations=IterationsConfig(
            sfc_flux=1, seb_bracket=1, seb_root=1,
            smb_flux=1, moisture_picard=1, coupling=1,
        ),
        tolerances=TolerancesConfig(
            sfc_flux=1e-6, seb_root=1e-12, smb_flux=1e-12,
            moisture_picard=1e-8, moisture_bounds=1e-12,
            coupling_temp=1e-6, coupling_mois=1e-9,
        ),
    )


def _stub(back_weight: float = 0.5, monotone: bool = True) -> Any:
    """A bare UtahLSM carrying only what _monotone_back_weight reads."""
    model: Any = UtahLSM.__new__(UtahLSM)
    model.logger = logging.getLogger("test")
    model.input = SimpleNamespace(numerics=_numerics(back_weight, monotone))
    model._monotone_back_weight_applied = None
    model._monotone_warning_issued = False
    return model


class TestMonotoneFloorFormula:
    """The floor is 1 - 1/(2 * max_row(alpha_up + alpha_dn))."""

    @pytest.mark.parametrize("alpha", [0.6, 1.0, 3.817, 25.0])
    def test_raises_weight_to_the_analytic_floor(self, alpha: float) -> None:
        """Above threshold the returned weight matches the closed form."""
        model = _stub(back_weight=0.5)
        row_sum = np.full((4, 1), 2.0 * alpha)

        theta_b = model._monotone_back_weight(0.5, row_sum)

        assert theta_b == pytest.approx(1.0 - 1.0 / (2.0 * 2.0 * alpha))

    @pytest.mark.parametrize("alpha", [0.05, 0.25, 0.4999])
    def test_is_a_noop_below_threshold(self, alpha: float) -> None:
        """Crank-Nicolson survives untouched wherever it is safe."""
        model = _stub(back_weight=0.5)
        row_sum = np.full((4, 1), 2.0 * alpha)

        assert model._monotone_back_weight(0.5, row_sum) == 0.5
        assert model._monotone_back_weight_applied is None

    def test_crank_nicolson_threshold_sits_at_alpha_one_half(self) -> None:
        """alpha = 0.5 is the exact CN boundary."""
        model = _stub(back_weight=0.5)

        below = model._monotone_back_weight(0.5, np.full((3, 1), 2.0 * 0.499))
        above = model._monotone_back_weight(0.5, np.full((3, 1), 2.0 * 0.501))

        assert below == 0.5
        assert above > 0.5

    def test_never_exceeds_fully_implicit(self) -> None:
        """The floor saturates at backward Euler."""
        model = _stub(back_weight=0.5)
        theta_b = model._monotone_back_weight(0.5, np.full((3, 1), 1e12))

        assert theta_b <= 1.0

    def test_backward_euler_is_already_monotone(self) -> None:
        """theta_b = 1 is never raised, at any diffusion number."""
        model = _stub(back_weight=1.0)
        assert model._monotone_back_weight(1.0, np.full((3, 1), 1e6)) == 1.0

    def test_uses_the_worst_row_not_the_mean(self) -> None:
        """A single stiff layer governs the whole column."""
        model = _stub(back_weight=0.5)
        row_sum = np.array([[0.1], [0.1], [20.0], [0.1]])

        theta_b = model._monotone_back_weight(0.5, row_sum)

        assert theta_b == pytest.approx(1.0 - 1.0 / (2.0 * 20.0))

    def test_disabled_flag_holds_the_configured_weight(self) -> None:
        """heat_diffusion_monotone=False leaves the scheme exactly as set."""
        model = _stub(back_weight=0.5, monotone=False)

        assert model._monotone_back_weight(0.5, np.full((3, 1), 100.0)) == 0.5
        assert model._monotone_back_weight_applied is None

    def test_degenerate_row_sums_are_ignored(self) -> None:
        """Zero or empty conductivity must not produce a divide-by-zero."""
        model = _stub(back_weight=0.5)

        assert model._monotone_back_weight(0.5, np.zeros((3, 1))) == 0.5
        assert model._monotone_back_weight(0.5, np.zeros((0, 1))) == 0.5

    def test_warns_once_not_every_step(self, caplog: Any) -> None:
        """A 1440-step run must not emit 1440 identical warnings."""
        model = _stub(back_weight=0.5)
        row_sum = np.full((3, 1), 8.0)

        with caplog.at_level(logging.WARNING):
            for _ in range(5):
                model._monotone_back_weight(0.5, row_sum)

        assert len(caplog.records) == 1
        assert "heat_diffusion_back_weight" in caplog.records[0].message


class TestSawtoothDecay:
    """End-to-end: the grid-scale mode must not alternate in sign."""

    @staticmethod
    def _sawtooth_history(monotone: bool, steps: int = 8) -> np.ndarray:
        """Seeds a +/-1 K grid-scale mode and tracks its projected amplitude."""
        import utahlsm

        inp = _load_case("gabls3")
        inp.numerics = replace(
            inp.numerics, heat_diffusion_monotone=monotone
        )
        assert inp.forcing is not None
        lsm = UtahLSM(inp, utahlsm.Output("/dev/null", enabled=False))
        lsm.tstep = inp.forcing.tstep

        field = lsm._ensure_writable_soil_field("temperature")
        base = np.array(field, copy=True)
        saw = ((-1.0) ** np.arange(field.shape[0]))[:, None]
        field[:] = base + saw
        lsm.sfc_state.soil_top_temperature = np.array(field[0], copy=True)

        amplitude = []
        for _ in range(steps):
            lsm._solve_soil_heat()
            delta = np.asarray(lsm.soil_state.temperature) - base
            amplitude.append(
                float(
                    np.sum(delta[1:-1] * saw[1:-1])
                    / np.sum(saw[1:-1] ** 2)
                )
            )
            lsm._ensure_writable_soil_field("temperature")
        return np.array(amplitude)

    @staticmethod
    def _sign_flips(series: np.ndarray) -> int:
        """Counts sign reversals between consecutive samples."""
        return int(np.sum(np.sign(series[1:]) * np.sign(series[:-1]) < 0))

    def test_pure_crank_nicolson_rings_at_gabls3_timestep(self) -> None:
        """Baseline: without the floor the mode flips sign every step."""
        amplitude = self._sawtooth_history(monotone=False)

        assert self._sign_flips(amplitude) >= len(amplitude) - 2

    def test_floor_removes_the_oscillation(self) -> None:
        """With the floor the mode decays without alternating."""
        amplitude = self._sawtooth_history(monotone=True)

        assert self._sign_flips(amplitude) <= 1
        assert abs(amplitude[-1]) < abs(amplitude[0])


class TestCaseConfigurations:
    """Which bundled cases actually trip the floor."""

    def test_gabls3_engages_the_floor(self) -> None:
        """600 s on a 1 cm grid needs far more implicitness than CN."""
        import utahlsm

        inp = _load_case("gabls3")
        assert inp.forcing is not None
        lsm = UtahLSM(inp, utahlsm.Output("/dev/null", enabled=False))
        lsm.tstep = inp.forcing.tstep
        lsm._solve_soil_heat()

        applied = lsm._monotone_back_weight_applied
        assert applied is not None
        assert applied > inp.numerics.heat_diffusion_back_weight
        assert applied == pytest.approx(0.93, abs=0.02)

    def test_arm_keeps_crank_nicolson(self) -> None:
        """60 s on the same grid stays comfortably below alpha = 0.5."""
        import utahlsm

        inp = _load_case("arm")
        assert inp.forcing is not None
        lsm = UtahLSM(inp, utahlsm.Output("/dev/null", enabled=False))
        lsm.tstep = inp.forcing.tstep
        lsm._solve_soil_heat()

        assert lsm._monotone_back_weight_applied is None
