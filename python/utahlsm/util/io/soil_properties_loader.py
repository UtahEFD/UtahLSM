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
"""Loader for externalized soil property datasets.

This module provides functionality to load soil property datasets from JSON files.
Properties can be loaded from public bundled datasets (referenced by name) or
from custom files specified by file path.

A dataset may also reference other datasets via an ``includes`` field; the
loader resolves these recursively and merges their ``soil_types`` tables. The
including dataset's own entries take precedence (i.e. it is merged last). When
the same soil-type name appears in multiple sources the ``on_conflict`` policy
applies: ``error`` (default), ``prefer_first``, or ``prefer_last``.
"""

import json
from importlib import resources
from pathlib import Path
from typing import Any, Protocol, cast

import jsonschema

from ...exceptions import NamelistError


class ResourcePath(Protocol):
    """Minimal packaged-resource interface used by this loader."""

    @property
    def name(self) -> str:
        """Return the resource base name."""
        ...

    def is_file(self) -> bool:
        """Return whether this resource points to a file."""
        ...

    def read_text(self, encoding: str | None = None) -> str:
        """Read the resource as text."""
        ...


class SoilPropertiesLoader:
    """Loads and validates soil property datasets from JSON files.

    Supports loading from:
    - Public bundled hydraulic datasets referenced by name (e.g., 'cosby')
    - Custom files specified by absolute or relative path
    - Composite custom datasets that ``include`` other datasets
    """

    # User-facing bundled dataset names (must match filenames in
    # utahlsm/data/soil/). Internal supplements such as Letts peat properties
    # are intentionally omitted from this public list.
    PUBLIC_BUNDLED_DATASETS = (
        "clapp-hornberger",
        "cosby",
        "rawls-brakensiek",
        "carsel-parrish",
    )
    BUNDLED_DATASETS: list[str] = list(PUBLIC_BUNDLED_DATASETS)
    _INTERNAL_BUNDLED_DATASETS = (
        *PUBLIC_BUNDLED_DATASETS,
        "letts",
    )
    MINERAL_BASE_DATASETS = set(PUBLIC_BUNDLED_DATASETS)
    # Peters-Lidard et al. (1998) texture-class quartz fractions used by
    # the Noah/Peters-Lidard Johansen thermal-conductivity implementation.
    PETERS_LIDARD_QUARTZ_FRACTIONS = {
        "sand": 0.92,
        "loamy_sand": 0.82,
        "sandy_loam": 0.60,
        "silty_loam": 0.25,
        "silt": 0.10,
        "loam": 0.40,
        "sandy_clay_loam": 0.60,
        "silty_clay_loam": 0.10,
        "clay_loam": 0.35,
        "sandy_clay": 0.52,
        "silty_clay": 0.10,
        "clay": 0.25,
    }

    @staticmethod
    def load(properties_spec: str) -> dict[str, dict[str, Any]]:
        """Load soil properties from a bundled dataset, file, or composite.

        Args:
            properties_spec: Either a bundled dataset name (e.g., 'cosby')
                or a file path to a custom JSON property file.

        Returns:
            Dictionary mapping soil type names (lowercase) to property dicts.
            Each property dict contains keys: 'b', 'psi_sat', 'porosity',
            'residual', 'K_sat', 'ci'.

        Raises:
            NamelistError: If dataset/file not found, file is invalid JSON,
                properties fail schema validation, includes form a cycle, or
                a soil-type conflict cannot be resolved by the on_conflict
                policy.
        """
        merged, _ = SoilPropertiesLoader._resolve(
            properties_spec, _seen=set(), _allow_internal=False
        )
        return merged

    @staticmethod
    def get_bundled_datasets() -> list[str]:
        """Return list of available bundled dataset names."""
        return SoilPropertiesLoader.BUNDLED_DATASETS.copy()

    @staticmethod
    def is_file_spec(spec: str) -> bool:
        """Return True when ``spec`` denotes a file path, not a bundled name."""
        return _is_file_path(spec)

    @staticmethod
    def _resolve(
        spec: str,
        _seen: set[str],
        _allow_internal: bool = False,
    ) -> tuple[dict[str, dict[str, Any]], str]:
        """Recursively resolve a dataset spec to merged soil_types.

        Returns the merged ``soil_types`` table and a canonical key used for
        cycle detection.
        """
        canonical = SoilPropertiesLoader._canonical_key(spec)
        if canonical in _seen:
            chain = " -> ".join([*_seen, canonical])
            raise NamelistError(
                f"Circular include detected in soil property datasets: {chain}"
            )
        _seen = _seen | {canonical}

        data, source = SoilPropertiesLoader._load_raw(spec, _allow_internal)
        SoilPropertiesLoader._validate(data, source)

        on_conflict = data.get("on_conflict", "error")
        merged: dict[str, dict[str, Any]] = {}

        for include_spec in data.get("includes", []):
            child_types, _ = SoilPropertiesLoader._resolve(
                include_spec, _seen, _allow_internal=_allow_internal
            )
            SoilPropertiesLoader._merge(
                merged, child_types, on_conflict, source, include_spec
            )

        own_types = cast(dict[str, dict[str, Any]], data.get("soil_types", {}))
        if own_types:
            SoilPropertiesLoader._normalize_retention(own_types)
            SoilPropertiesLoader._merge(merged, own_types, on_conflict, source, source)

        if not _is_file_path(spec):
            SoilPropertiesLoader._apply_bundled_supplements(spec, merged, _seen)

        if not merged:
            raise NamelistError(
                f"Soil property dataset {source} resolved to no soil types"
            )

        return merged, canonical

    @staticmethod
    def _normalize_retention(soil_types: dict[str, dict[str, Any]]) -> None:
        """Derive Campbell-style fields for native van Genuchten entries.

        Dataset entries may specify retention either as ('b', 'psi_sat') or
        as native van Genuchten ('alpha' [1/m], 'n'); the schema enforces
        exactly one pair. For native entries the internal equivalents are
        exact under the model's Mualem-constrained parameterization
        (m = 1/(1+b) = 1-1/n and air entry psi_sat = -1/alpha):

            b = 1/(n-1),    psi_sat = -1/alpha

        The entry is tagged with ``parameterization: 'van-genuchten'`` so the
        soil-model factory can reject these types for Campbell/Brooks-Corey,
        where the derived 'b' is not a measured pore-size index.
        """
        for props in soil_types.values():
            if "alpha" in props:
                alpha = float(props["alpha"])
                n = float(props["n"])
                props["b"] = 1.0 / (n - 1.0)
                props["psi_sat"] = -1.0 / alpha
                props["parameterization"] = "van-genuchten"

    @staticmethod
    def _merge(
        target: dict[str, dict[str, Any]],
        incoming: dict[str, dict[str, Any]],
        policy: str,
        owner: str,
        source: str,
    ) -> None:
        """Merge ``incoming`` soil types into ``target`` using ``policy``."""
        for soil_type, props in incoming.items():
            if soil_type not in target:
                target[soil_type] = props
                continue
            if policy == "prefer_first":
                continue
            if policy == "prefer_last":
                target[soil_type] = props
                continue
            raise NamelistError(
                f"Soil type '{soil_type}' appears in more than one source "
                f"while resolving dataset {owner} (conflict from {source}). "
                f"Set 'on_conflict' to 'prefer_first' or 'prefer_last' to "
                f"resolve."
            )

    @staticmethod
    def _apply_bundled_supplements(
        dataset_name: str,
        soil_types: dict[str, dict[str, Any]],
        seen: set[str],
    ) -> None:
        """Attach bundled thermal/organic data to mineral base datasets.

        The user-facing bundled choices stay as the three mineral hydraulic
        datasets. Internally they carry the extra data required by supported
        thermal and organic soil options:

        * Peters-Lidard texture-class quartz fractions for Johansen thermal
          conductivity.
        * Letts peat tiers for organic layers.
        """
        if dataset_name not in SoilPropertiesLoader.MINERAL_BASE_DATASETS:
            return

        for (
            soil_type,
            quartz_fraction,
        ) in SoilPropertiesLoader.PETERS_LIDARD_QUARTZ_FRACTIONS.items():
            if soil_type in soil_types:
                soil_types[soil_type].setdefault("quartz_fraction", quartz_fraction)

        letts_types, _ = SoilPropertiesLoader._resolve(
            "letts", seen, _allow_internal=True
        )
        for soil_type, props in letts_types.items():
            soil_types.setdefault(soil_type, props)

    @staticmethod
    def _load_raw(
        spec: str, _allow_internal: bool = False
    ) -> tuple[dict[str, Any], str]:
        """Load and JSON-parse a dataset spec without validation or merging."""
        if _is_file_path(spec):
            return SoilPropertiesLoader._load_file_raw(spec)
        return SoilPropertiesLoader._load_bundled_raw(spec, _allow_internal)

    @staticmethod
    def _load_bundled_raw(
        dataset_name: str, _allow_internal: bool = False
    ) -> tuple[dict[str, Any], str]:
        """Load a bundled dataset by name (no validation)."""
        allowed_datasets = (
            SoilPropertiesLoader._INTERNAL_BUNDLED_DATASETS
            if _allow_internal
            else SoilPropertiesLoader.PUBLIC_BUNDLED_DATASETS
        )
        if dataset_name not in allowed_datasets:
            available = ", ".join(SoilPropertiesLoader.PUBLIC_BUNDLED_DATASETS)
            raise NamelistError(
                f"Soil property dataset {dataset_name} not found. "
                f"Available bundled datasets: {available}"
            )

        bundled_path = SoilPropertiesLoader._get_bundled_path(dataset_name)
        source = f"utahlsm/data/soil/{dataset_name}.json"

        if not bundled_path.is_file():
            raise NamelistError(
                f"Bundled soil property file not found: {bundled_path}\n"
                f"Expected location: {source}"
            )

        try:
            data = json.loads(bundled_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as e:
            raise NamelistError(
                f"Invalid JSON in soil property resource {source}: {e}"
            ) from e
        except OSError as e:
            raise NamelistError(
                f"Error reading soil property resource {source}: {e}"
            ) from e
        return data, source

    @staticmethod
    def _load_file_raw(file_path: str) -> tuple[dict[str, Any], str]:
        """Load a JSON file from disk (no validation)."""
        path = Path(file_path).expanduser().resolve()

        if not path.exists():
            raise NamelistError(
                f"Soil property file not found: {file_path}\nResolved to: {path}"
            )
        if not path.is_file():
            raise NamelistError(f"Path is not a file: {path}")

        try:
            with open(path, encoding="utf-8") as f:
                data = json.load(f)
        except json.JSONDecodeError as e:
            raise NamelistError(
                f"Invalid JSON in soil property file {path}: {e}"
            ) from e
        except OSError as e:
            raise NamelistError(f"Error reading soil property file {path}: {e}") from e
        return data, str(path)

    @staticmethod
    def _validate(data: dict[str, Any], source: str = "properties") -> None:
        """Validate loaded properties against schema."""
        try:
            schema_resource = resources.files("utahlsm.util.io").joinpath(
                "schema_soil_properties.json"
            )
            schema = json.loads(schema_resource.read_text(encoding="utf-8"))
        except Exception as e:
            raise NamelistError(f"Failed to load soil properties schema: {e}") from e

        try:
            jsonschema.validate(instance=data, schema=schema)
        except jsonschema.ValidationError as e:
            raise NamelistError(
                f"Soil properties validation failed for {source}:\n"
                f"  Path: {list(e.path)}\n"
                f"  Message: {e.message}"
            ) from e
        except jsonschema.SchemaError as e:
            raise NamelistError(
                f"Internal error: soil properties schema is invalid: {e}"
            ) from e

    @staticmethod
    def _get_bundled_path(dataset_name: str) -> ResourcePath:
        """Get the packaged resource for a bundled dataset file."""
        return (
            resources.files("utahlsm")
            .joinpath("data")
            .joinpath("soil")
            .joinpath(f"{dataset_name}.json")
        )

    @staticmethod
    def _canonical_key(spec: str) -> str:
        """Return a canonical key for cycle detection."""
        if _is_file_path(spec):
            return f"file:{Path(spec).expanduser().resolve()}"
        return f"bundled:{spec}"


def _is_file_path(spec: str) -> bool:
    """Determine if a spec string is a file path or dataset name.

    A spec is considered a file path if it contains path separators or
    if it ends with '.json'.
    """
    return "/" in spec or "\\" in spec or spec.endswith(".json")
