"""Tests for the composite-include mechanism in SoilPropertiesLoader."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from utahlsm.exceptions import NamelistError
from utahlsm.util.io.soil_properties_loader import SoilPropertiesLoader


def test_cosby_letts_resolves_mineral_and_peat() -> None:
    """cosby-letts should expose Cosby mineral classes plus Letts peat tiers."""
    types = SoilPropertiesLoader.load('cosby-letts')

    for mineral in ('sand', 'loam', 'clay'):
        assert mineral in types
    for peat in ('peat_fibric', 'peat_hemic', 'peat_sapric'):
        assert peat in types

    assert types['clay']['porosity'] == pytest.approx(0.468)
    assert types['peat_hemic']['ci'] == pytest.approx(2.5e6)


def test_clapp_hornberger_and_rawls_composites_resolve() -> None:
    """The other shipped composites must also resolve cleanly."""
    for composite in ('clapp-hornberger-letts', 'rawls-brakensiek-letts'):
        types = SoilPropertiesLoader.load(composite)
        assert 'clay' in types
        assert 'peat_sapric' in types


def test_conflict_default_raises(tmp_path: Path) -> None:
    """Default on_conflict='error' raises when sources overlap."""
    a = tmp_path / 'a.json'
    b = tmp_path / 'b.json'
    composite = tmp_path / 'c.json'

    base = {'b': 1.0, 'psi_sat': -0.1, 'porosity': 0.4,
            'residual': 0.0, 'K_sat': 1e-5, 'ci': 1.0e6}
    a.write_text(json.dumps({
        'metadata': {'dataset_name': 'a'},
        'soil_types': {'loam': base},
    }))
    b.write_text(json.dumps({
        'metadata': {'dataset_name': 'b'},
        'soil_types': {'loam': {**base, 'b': 2.0}},
    }))
    composite.write_text(json.dumps({
        'metadata': {'dataset_name': 'c'},
        'includes': [str(a), str(b)],
    }))

    with pytest.raises(NamelistError, match='appears in more than one source'):
        SoilPropertiesLoader.load(str(composite))


def test_conflict_prefer_last_overrides(tmp_path: Path) -> None:
    """on_conflict='prefer_last' lets later sources override earlier ones."""
    a = tmp_path / 'a.json'
    b = tmp_path / 'b.json'
    composite = tmp_path / 'c.json'

    base = {'b': 1.0, 'psi_sat': -0.1, 'porosity': 0.4,
            'residual': 0.0, 'K_sat': 1e-5, 'ci': 1.0e6}
    a.write_text(json.dumps({
        'metadata': {'dataset_name': 'a'},
        'soil_types': {'loam': base},
    }))
    b.write_text(json.dumps({
        'metadata': {'dataset_name': 'b'},
        'soil_types': {'loam': {**base, 'b': 99.0}},
    }))
    composite.write_text(json.dumps({
        'metadata': {'dataset_name': 'c'},
        'includes': [str(a), str(b)],
        'on_conflict': 'prefer_last',
    }))

    types = SoilPropertiesLoader.load(str(composite))
    assert types['loam']['b'] == pytest.approx(99.0)


def test_circular_include_detected(tmp_path: Path) -> None:
    """Cycles in includes must raise rather than recurse forever."""
    a = tmp_path / 'a.json'
    b = tmp_path / 'b.json'

    a.write_text(json.dumps({
        'metadata': {'dataset_name': 'a'},
        'includes': [str(b)],
    }))
    b.write_text(json.dumps({
        'metadata': {'dataset_name': 'b'},
        'includes': [str(a)],
    }))

    with pytest.raises(NamelistError, match='Circular include'):
        SoilPropertiesLoader.load(str(a))


def test_dataset_with_only_includes_has_no_own_types() -> None:
    """A pure manifest (includes only, no soil_types) is valid."""
    types = SoilPropertiesLoader.load('cosby-letts')
    assert len(types) > 0
