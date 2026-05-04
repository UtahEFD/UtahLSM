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
Properties can be loaded from bundled default datasets (referenced by name) or
from custom files specified by file path.
"""

import json
from importlib import resources
from importlib.abc import Traversable
from pathlib import Path
from typing import Any

import jsonschema

from ...exceptions import NamelistError


class SoilPropertiesLoader:
    """Loads and validates soil property datasets from JSON files.

    Supports loading from:
    - Bundled datasets referenced by name (e.g., 'cosby-1984')
    - Custom files specified by absolute or relative path
    """

    # Bundled dataset names (must match filenames in utahlsm/data/soil/)
    BUNDLED_DATASETS = [
        'clapp-hornberger',
        'cosby',
        'rawls-brakensiek',
        'cabauw-heinen'
    ]

    @staticmethod
    def load(properties_spec: str) -> dict[str, dict[str, float]]:
        """Load soil properties from bundled dataset or custom file.

        Args:
            properties_spec: Either a bundled dataset name (e.g., 'cosby-1984')
                or a file path to a custom JSON property file.

        Returns:
            Dictionary mapping soil type names (lowercase) to property dicts.
            Each property dict contains keys: 'b', 'psi_sat', 'porosity',
            'residual', 'K_sat', 'ci'.

        Raises:
            NamelistError: If dataset/file not found, file is invalid JSON,
                or properties fail schema validation.
        """
        # Determine if spec is a file path or bundled dataset name
        if _is_file_path(properties_spec):
            return SoilPropertiesLoader._load_from_file(properties_spec)

        return SoilPropertiesLoader._load_bundled(properties_spec)

    @staticmethod
    def get_bundled_datasets() -> list[str]:
        """Return list of available bundled dataset names."""
        return SoilPropertiesLoader.BUNDLED_DATASETS.copy()

    @staticmethod
    def _load_bundled(dataset_name: str) -> dict[str, dict[str, float]]:
        """Load a bundled dataset by name.

        Args:
            dataset_name: Name of bundled dataset (e.g., 'cosby-1984')

        Returns:
            Dictionary mapping soil type names to property dicts.

        Raises:
            NamelistError: If dataset not found or fails validation.
        """
        if dataset_name not in SoilPropertiesLoader.BUNDLED_DATASETS:
            available = ', '.join(SoilPropertiesLoader.BUNDLED_DATASETS)
            raise NamelistError(
                f'Soil property dataset {dataset_name} not found. '
                f'Available bundled datasets: {available}'
            )

        bundled_path: Traversable = SoilPropertiesLoader._get_bundled_path(dataset_name)

        if not bundled_path.is_file():
            raise NamelistError(
                f'Bundled soil property file not found: {bundled_path}\n'
                f'Expected location: utahlsm/data/soil/{dataset_name}.json'
            )

        return SoilPropertiesLoader._load_from_resource(
            bundled_path, f'utahlsm/data/soil/{dataset_name}.json')

    @staticmethod
    def _load_from_file(file_path: str) -> dict[str, dict[str, float]]:
        """Load properties from a JSON file (bundled or custom).

        Args:
            file_path: Path to JSON file (absolute or relative).

        Returns:
            Dictionary mapping soil type names to property dicts.

        Raises:
            NamelistError: If file not found, is invalid JSON, or fails validation.
        """
        path = Path(file_path).expanduser().resolve()

        if not path.exists():
            raise NamelistError(
                f'Soil property file not found: {file_path}\n'
                f'Resolved to: {path}'
            )

        if not path.is_file():
            raise NamelistError(f'Path is not a file: {path}')

        try:
            with open(path, encoding='utf-8') as f:
                data = json.load(f)
        except json.JSONDecodeError as e:
            raise NamelistError(
                f'Invalid JSON in soil property file {path}: {e}'
            ) from e
        except OSError as e:
            raise NamelistError(
                f'Error reading soil property file {path}: {e}'
            ) from e

        SoilPropertiesLoader._validate(data, str(path))
        return data['soil_types']

    @staticmethod
    def _load_from_resource(resource: Traversable, source: str) -> dict[str, dict[str, float]]:
        """Load packaged JSON resources from the installed utahlsm package.

        Args:
            resource: Packaged resource object returned by
                ``importlib.resources.files(...).joinpath(...)``.
            source: Resource description used in error messages.

        Returns:
            Dictionary mapping soil type names to property dicts.

        Raises:
            NamelistError: If the resource cannot be read, is invalid JSON, or
                fails schema validation.
        """
        try:
            with resource.open('r', encoding='utf-8') as f:
                data: dict[str, Any] = json.load(f)
        except json.JSONDecodeError as e:
            raise NamelistError(
                f'Invalid JSON in soil property resource {source}: {e}'
            ) from e
        except OSError as e:
            raise NamelistError(
                f'Error reading soil property resource {source}: {e}'
            ) from e

        SoilPropertiesLoader._validate(data, source)
        return data['soil_types']

    @staticmethod
    def _validate(data: dict[str, Any], source: str = 'properties') -> None:
        """Validate loaded properties against schema.

        Args:
            data: Loaded JSON data (should contain 'metadata' and 'soil_types')
            source: Description of source for error messages

        Raises:
            NamelistError: If data fails schema validation.
        """
        try:
            schema_resource = resources.files(__package__).joinpath(
                'schema_soil_properties.json')
            with schema_resource.open('r', encoding='utf-8') as f:
                schema = json.load(f)
        except Exception as e:
            raise NamelistError(
                f'Failed to load soil properties schema: {e}'
            ) from e

        try:
            jsonschema.validate(instance=data, schema=schema)
        except jsonschema.ValidationError as e:
            raise NamelistError(
                f'Soil properties validation failed for {source}:\n'
                f'  Path: {list(e.path)}\n'
                f'  Message: {e.message}'
            ) from e
        except jsonschema.SchemaError as e:
            raise NamelistError(
                f'Internal error: soil properties schema is invalid: {e}'
            ) from e

    @staticmethod
    def _get_bundled_path(dataset_name: str) -> Traversable:
        """Get the packaged resource for a bundled dataset file.

        Args:
            dataset_name: Name of bundled dataset (e.g., 'cosby')

        Returns:
            Traversable resource for the JSON file in utahlsm/data/soil/
        """
        return resources.files('utahlsm').joinpath(
            'data', 'soil', f'{dataset_name}.json')


def _is_file_path(spec: str) -> bool:
    """Determine if a spec string is a file path or dataset name.

    A spec is considered a file path if it contains path separators or
    if it ends with '.json'.

    Args:
        spec: Specification string

    Returns:
        True if spec appears to be a file path, False otherwise.
    """
    return '/' in spec or '\\' in spec or spec.endswith('.json')
