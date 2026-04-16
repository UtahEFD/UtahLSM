"""Integration test for GABLS3 2x2 regression."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from tests import compare_gabls3_2x2 as cmp


@pytest.mark.integration
@pytest.mark.slow
def test_gabls3_2x2_matches_single(tmp_path: Path) -> None:
    """Ensures each 2x2 column matches the single-column GABLS3 output."""
    repo_root = Path(__file__).resolve().parents[1]
    cases_root = (repo_root / ".." / "cases").resolve()
    case_single = cases_root / "gabls3"
    case_multi = cases_root / "gabls3_2x2"

    if not case_single.exists() or not case_multi.exists():
        pytest.skip("GABLS3 case files not available.")

    single_out = tmp_path / "lsm_gabls3_py.nc"
    multi_out = tmp_path / "lsm_gabls3_2x2_py.nc"

    subprocess.run(
        [sys.executable, "utahlsm_offline.py", "-c", "gabls3", "-o", str(single_out)],
        cwd=repo_root,
        check=True,
    )
    subprocess.run(
        [sys.executable, "utahlsm_offline.py", "-c", "gabls3_2x2", "-o", str(multi_out)],
        cwd=repo_root,
        check=True,
    )

    # Absolute tolerance is set well below any physically meaningful
    # signal; it just absorbs the ULP-level FP reordering between the
    # scalar-like (ncol=1) and broadcast (ncol=4) code paths in the
    # canopy partition arithmetic.
    results = cmp.compare_outputs(
        single_path=single_out,
        multi_path=multi_out,
        rtol=0.0,
        atol=1e-10,
    )

    failures = [res for res in results if not res.matches]
    if failures:
        details = ", ".join(
            f"{res.name} (max diff {res.max_abs_diff:.3e})"
            for res in failures
        )
        raise AssertionError(f"GABLS3 2x2 mismatch: {details}")

    for res in results:
        if res.max_abs_col_diff is not None:
            assert res.max_abs_col_diff <= 1e-12, (
                f"{res.name} columns differ by {res.max_abs_col_diff:.3e}"
            )
