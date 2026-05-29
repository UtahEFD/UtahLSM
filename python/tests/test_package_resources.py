"""Packaging and packaged-resource regression tests."""

from __future__ import annotations

from importlib import resources
from pathlib import Path

from utahlsm.util.io.soil_properties_loader import SoilPropertiesLoader


def test_packaged_resources_exist() -> None:
    """Checks that packaged schemas, soil datasets, and typing marker exist."""
    package_root = resources.files("utahlsm")
    io_root = package_root.joinpath("util").joinpath("io")
    soil_root = package_root.joinpath("data").joinpath("soil")

    assert package_root.joinpath("py.typed").is_file()
    assert io_root.joinpath("schema_namelist.json").is_file()
    assert io_root.joinpath("schema_soil_properties.json").is_file()

    for dataset in SoilPropertiesLoader.get_bundled_datasets():
        assert soil_root.joinpath(f"{dataset}.json").is_file()
    assert soil_root.joinpath("letts.json").is_file()


def test_bundled_soil_resource_is_packaged() -> None:
    """Ensures bundled soil datasets resolve through packaged resources."""
    bundled = SoilPropertiesLoader._get_bundled_path("cosby")

    assert bundled.is_file()
    assert bundled.name == "cosby.json"


def test_pyproject_discovers_all_utahlsm_packages() -> None:
    """Ensures setuptools package discovery includes all utahlsm subpackages."""
    pyproject = Path(__file__).resolve().parents[1] / "pyproject.toml"
    contents = pyproject.read_text(encoding="utf-8")

    assert "[tool.setuptools.packages.find]" in contents
    assert 'include = ["utahlsm*"]' in contents
    assert '"data/soil/*.json"' in contents
    assert '"util/io/*.json"' in contents
    assert '"py.typed"' in contents
