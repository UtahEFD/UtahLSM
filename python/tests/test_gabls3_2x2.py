"""Integration test for temporary 2x2 GABLS3 regression coverage."""

from __future__ import annotations

from pathlib import Path

import pytest

from tests import compare_gabls3_2x2 as cmp
from tests.gabls3_case_factory import run_case, write_2x2_case


@pytest.mark.integration
@pytest.mark.slow
def test_gabls3_2x2_matches_single(tmp_path: Path) -> None:
    """Ensures the vectorized 2x2 path matches a tiled single-column case."""
    repo_root = Path(__file__).resolve().parents[1]
    cases_root = (repo_root / ".." / "cases").resolve()
    case_single = cases_root / "gabls3"

    if not case_single.exists():
        pytest.skip("GABLS3 case files not available.")

    case_multi = tmp_path / "gabls3_2x2_generated"
    single_out = tmp_path / "lsm_gabls3_py.nc"
    multi_out = tmp_path / "lsm_gabls3_2x2_py.nc"

    write_2x2_case(case_single, case_multi)
    run_case(case_single, single_out)
    run_case(case_multi, multi_out)

    # Absolute tolerance is set well below any physically meaningful
    # signal; it just absorbs the ULP-level FP reordering between the
    # scalar-like (ncol=1) and broadcast (ncol=4) code paths in the
    # canopy partition arithmetic.
    results = cmp.compare_outputs(
        single_path=single_out,
        multi_path=multi_out,
        rtol=0.0,
        atol=1e-9,
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
