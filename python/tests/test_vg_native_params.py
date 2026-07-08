"""Tests for native van Genuchten (alpha/n) soil property parameters.

Datasets may define retention either as Campbell-style ('b', 'psi_sat') or
as native van Genuchten ('alpha' [1/m], 'n'). The loader derives the exact
internal equivalents b = 1/(n-1) and psi_sat = -1/alpha (the model's
Mualem-constrained parameterization makes this mapping exact) and tags the
entry so non-VG soil models reject it.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from utahlsm.exceptions import NamelistError
from utahlsm.physics.soil.factory import get_soil_model
from utahlsm.util.io.soil_properties_loader import SoilPropertiesLoader


def _native_entry(**overrides: float) -> dict[str, float]:
    entry = {
        "alpha": 2.0,
        "n": 1.41,
        "porosity": 0.45,
        "residual": 0.067,
        "K_sat": 1.25e-6,
        "ci": 1.27e6,
    }
    entry.update(overrides)
    return entry


def _write_dataset(path: Path, soil_types: dict[str, object]) -> str:
    path.write_text(
        json.dumps(
            {
                "metadata": {"dataset_name": "test-native-vg"},
                "soil_types": soil_types,
            }
        ),
        encoding="utf-8",
    )
    return str(path)


@pytest.mark.soil
class TestNativeRetentionLoading:
    """Loader derives exact Campbell-style equivalents for alpha/n entries."""

    def test_native_entry_derives_b_and_psi_sat(self, tmp_path: Path) -> None:
        spec = _write_dataset(tmp_path / "vg.json", {"silty_loam": _native_entry()})
        types = SoilPropertiesLoader.load(spec)
        props = types["silty_loam"]
        assert props["b"] == pytest.approx(1.0 / (1.41 - 1.0))
        assert props["psi_sat"] == pytest.approx(-1.0 / 2.0)
        assert props["parameterization"] == "van-genuchten"

    def test_campbell_entry_is_untagged(self, tmp_path: Path) -> None:
        entry = {
            "b": 5.3,
            "psi_sat": -0.786,
            "porosity": 0.485,
            "residual": 0.0,
            "K_sat": 7.2e-6,
            "ci": 1.27e6,
        }
        spec = _write_dataset(tmp_path / "cb.json", {"silty_loam": entry})
        types = SoilPropertiesLoader.load(spec)
        assert "parameterization" not in types["silty_loam"]

    def test_mixed_retention_pairs_rejected(self, tmp_path: Path) -> None:
        spec = _write_dataset(
            tmp_path / "mixed.json",
            {"silty_loam": _native_entry(b=5.3, psi_sat=-0.786)},
        )
        with pytest.raises(NamelistError, match="validation failed"):
            SoilPropertiesLoader.load(spec)

    def test_incomplete_native_pair_rejected(self, tmp_path: Path) -> None:
        entry = _native_entry()
        del entry["n"]
        spec = _write_dataset(tmp_path / "incomplete.json", {"silty_loam": entry})
        with pytest.raises(NamelistError, match="validation failed"):
            SoilPropertiesLoader.load(spec)

    def test_carsel_parrish_bundled_dataset(self) -> None:
        types = SoilPropertiesLoader.load("carsel-parrish")
        usda_classes = (
            "sand", "loamy_sand", "sandy_loam", "loam", "silt", "silty_loam",
            "sandy_clay_loam", "clay_loam", "silty_clay_loam", "sandy_clay",
            "silty_clay", "clay",
        )
        for name in usda_classes:
            props = types[name]
            assert props["parameterization"] == "van-genuchten"
            assert props["b"] == pytest.approx(1.0 / (props["n"] - 1.0))
            assert props["psi_sat"] == pytest.approx(-1.0 / props["alpha"])
            assert "quartz_fraction" in props
        # Spot-check silt loam against Carsel & Parrish (1988) Table 3.
        assert types["silty_loam"]["alpha"] == pytest.approx(2.0)
        assert types["silty_loam"]["n"] == pytest.approx(1.41)
        assert types["silty_loam"]["K_sat"] == pytest.approx(0.45e-2 / 3600.0,
                                                             rel=1e-3)


@pytest.mark.soil
class TestNativeRetentionModels:
    """The VG model reproduces the analytic van Genuchten curve; others gate."""

    def test_vg_matches_analytic_curve(self, tmp_path: Path) -> None:
        spec = _write_dataset(tmp_path / "vg.json", {"silty_loam": _native_entry()})
        types = SoilPropertiesLoader.load(spec)
        model = get_soil_model("van-genuchten", types, ["silty_loam"])

        alpha, n = 2.0, 1.41
        m = 1.0 - 1.0 / n
        porosity, residual = 0.45, 0.067
        theta = np.linspace(0.1, 0.43, 8)
        Se = (theta - residual) / (porosity - residual)
        psi_expected = -(1.0 / alpha) * (Se ** (-1.0 / m) - 1.0) ** (1.0 / n)
        psi_model = np.array(
            [model.water_potential(t, level=0) for t in theta]
        )
        np.testing.assert_allclose(psi_model, psi_expected, rtol=1e-12)

    def test_campbell_rejects_native_types(self, tmp_path: Path) -> None:
        spec = _write_dataset(tmp_path / "vg.json", {"silty_loam": _native_entry()})
        types = SoilPropertiesLoader.load(spec)
        with pytest.raises(NamelistError, match="native van Genuchten"):
            get_soil_model("campbell", types, ["silty_loam"])

    def test_brooks_corey_rejects_native_types(self, tmp_path: Path) -> None:
        spec = _write_dataset(tmp_path / "vg.json", {"silty_loam": _native_entry()})
        types = SoilPropertiesLoader.load(spec)
        with pytest.raises(NamelistError, match="native van Genuchten"):
            get_soil_model("brooks-corey", types, ["silty_loam"])

    def test_campbell_allows_campbell_types(self) -> None:
        types = SoilPropertiesLoader.load("clapp-hornberger")
        model = get_soil_model("campbell", types, ["silty_loam"])
        assert model is not None
