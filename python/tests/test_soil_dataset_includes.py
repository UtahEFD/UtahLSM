"""Tests for bundled soil supplements and custom dataset includes."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from utahlsm.exceptions import NamelistError
from utahlsm.util.io.soil_properties_loader import SoilPropertiesLoader


PUBLIC_DATASETS = ["clapp-hornberger", "cosby", "rawls-brakensiek"]


def _base_props(**overrides: float) -> dict[str, float]:
    props = {
        "b": 1.0,
        "psi_sat": -0.1,
        "porosity": 0.4,
        "residual": 0.0,
        "K_sat": 1e-5,
        "ci": 1.0e6,
    }
    props.update(overrides)
    return props


def _write_dataset(path: Path, data: dict[str, object]) -> None:
    path.write_text(json.dumps(data), encoding="utf-8")


def test_public_bundled_datasets_are_base_hydraulic_tables() -> None:
    """Users choose only from the three public base datasets."""
    assert SoilPropertiesLoader.get_bundled_datasets() == PUBLIC_DATASETS


@pytest.mark.parametrize("dataset", PUBLIC_DATASETS)
def test_bundled_base_datasets_resolve_quartz_and_peat(
    dataset: str,
) -> None:
    """Bundled bases carry internal quartz and Letts peat supplements."""
    types = SoilPropertiesLoader.load(dataset)

    for mineral in ("sand", "loam", "clay", "clay_loam"):
        assert mineral in types
    for peat in ("peat_fibric", "peat_hemic", "peat_sapric"):
        assert peat in types

    assert types["sand"]["quartz_fraction"] == pytest.approx(0.92)
    assert types["loam"]["quartz_fraction"] == pytest.approx(0.40)
    assert types["clay_loam"]["quartz_fraction"] == pytest.approx(0.35)
    assert types["clay"]["quartz_fraction"] == pytest.approx(0.25)
    assert types["peat_hemic"]["ci"] == pytest.approx(2.5e6)
    assert "quartz_fraction" not in types["peat_hemic"]


@pytest.mark.parametrize(
    "removed_alias",
    ["clapp-hornberger_letts", "cosby_letts", "rawls-brakensiek_letts", "letts"],
)
def test_old_composite_and_internal_names_are_not_public(
    removed_alias: str,
) -> None:
    """The namelist-facing bundled namespace is intentionally small."""
    with pytest.raises(NamelistError, match="Available bundled datasets"):
        SoilPropertiesLoader.load(removed_alias)


def test_conflict_default_raises(tmp_path: Path) -> None:
    """Default on_conflict='error' raises when sources overlap."""
    a = tmp_path / "a.json"
    b = tmp_path / "b.json"
    composite = tmp_path / "c.json"

    base = _base_props()
    _write_dataset(
        a,
        {
            "metadata": {"dataset_name": "a"},
            "soil_types": {"loam": base},
        },
    )
    _write_dataset(
        b,
        {
            "metadata": {"dataset_name": "b"},
            "soil_types": {"loam": {**base, "b": 2.0}},
        },
    )
    _write_dataset(
        composite,
        {
            "metadata": {"dataset_name": "c"},
            "includes": [str(a), str(b)],
        },
    )

    with pytest.raises(NamelistError, match="appears in more than one source"):
        SoilPropertiesLoader.load(str(composite))


def test_conflict_prefer_last_overrides(tmp_path: Path) -> None:
    """on_conflict='prefer_last' lets later sources override earlier ones."""
    a = tmp_path / "a.json"
    b = tmp_path / "b.json"
    composite = tmp_path / "c.json"

    base = _base_props()
    _write_dataset(
        a,
        {
            "metadata": {"dataset_name": "a"},
            "soil_types": {"loam": base},
        },
    )
    _write_dataset(
        b,
        {
            "metadata": {"dataset_name": "b"},
            "soil_types": {"loam": {**base, "b": 99.0}},
        },
    )
    _write_dataset(
        composite,
        {
            "metadata": {"dataset_name": "c"},
            "includes": [str(a), str(b)],
            "on_conflict": "prefer_last",
        },
    )

    types = SoilPropertiesLoader.load(str(composite))
    assert types["loam"]["b"] == pytest.approx(99.0)


def test_circular_include_detected(tmp_path: Path) -> None:
    """Cycles in includes must raise rather than recurse forever."""
    a = tmp_path / "a.json"
    b = tmp_path / "b.json"

    _write_dataset(
        a,
        {
            "metadata": {"dataset_name": "a"},
            "includes": [str(b)],
        },
    )
    _write_dataset(
        b,
        {
            "metadata": {"dataset_name": "b"},
            "includes": [str(a)],
        },
    )

    with pytest.raises(NamelistError, match="Circular include"):
        SoilPropertiesLoader.load(str(a))


def test_dataset_with_only_includes_has_no_own_types(tmp_path: Path) -> None:
    """A pure manifest (includes only, no soil_types) is valid."""
    base = tmp_path / "base.json"
    composite = tmp_path / "composite.json"

    _write_dataset(
        base,
        {
            "metadata": {"dataset_name": "base"},
            "soil_types": {"loam": _base_props(b=3.0)},
        },
    )
    _write_dataset(
        composite,
        {
            "metadata": {"dataset_name": "composite"},
            "includes": [str(base)],
        },
    )

    types = SoilPropertiesLoader.load(str(composite))
    assert len(types) > 0
    assert types["loam"]["b"] == pytest.approx(3.0)
