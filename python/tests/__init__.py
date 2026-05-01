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
"""UtahLSM Test Suite.

This package contains the comprehensive test suite for UtahLSM, organized
into unit tests for different model components.

Test Organization:
- test_solvers.py: Numerical solver unit tests (tridiagonal, root_brent)
- test_soil_models.py: Soil physics model tests
- test_seb_solver.py: Surface energy budget solver tests
- test_smb_solver.py: Surface moisture budget solver tests
- test_data_models.py: Configuration and state dataclass tests

To run all tests:
    pytest

To run tests by marker:
    pytest -m solver        # Run solver tests
    pytest -m soil          # Run soil physics tests
    pytest -m seb           # Run SEB tests
    pytest -m smb           # Run SMB tests
    pytest -m datamodel     # Run data model tests

To run with verbose output:
    pytest -v

To run a single test file:
    pytest tests/test_solvers.py
"""
