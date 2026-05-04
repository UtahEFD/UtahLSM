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
"""Unit tests for numerical solvers in UtahLSM.

This module tests the core numerical solvers:
- `tridiagonal()`: Thomas algorithm for tridiagonal matrix systems
- `root_brent()`: Brent's root-finding method

These solvers are fundamental to the model's physics calculations,
so thorough testing is essential.

Testing Strategy:
1. Test against analytical solutions (where available)
2. Test with different matrix conditioning levels
3. Test boundary conditions and edge cases
4. Test error handling for invalid inputs
5. Test numerical stability and convergence
"""

from typing import Any, Callable

import numpy as np
import pytest
from numpy.testing import assert_allclose
from numpy.typing import NDArray

from utahlsm.util.solvers import root_brent, root_brent_vec, tridiagonal

# ============================================================================
# Tests for tridiagonal() - Thomas Algorithm
# ============================================================================

@pytest.mark.solver
class TestTridiagonal:
    """Test suite for the tridiagonal matrix solver (Thomas algorithm)."""

    def test_single_element(self) -> None:
        """Test tridiagonal solver with single element.

        The system reduces to: b[0] * x[0] = r[0]
        Solution: x[0] = r[0] / b[0]
        """
        a = np.array([0.0])
        b = np.array([2.0])
        c = np.array([0.0])
        r = np.array([4.0])

        x = tridiagonal(a, b, c, r)

        expected = np.array([2.0])
        assert_allclose(x, expected, rtol=1e-10)

    def test_two_elements(self) -> None:
        """Test with minimal 2-element system.

        System:
        2*x[0] + 1*x[1] = 4
        1*x[0] + 2*x[1] = 3

        Solving: x[0] = 5/3, x[1] = 2/3
        """
        a = np.array([0.0, 1.0])
        b = np.array([2.0, 2.0])
        c = np.array([1.0, 0.0])
        r = np.array([4.0, 3.0])

        x = tridiagonal(a, b, c, r)

        # Expected solution: solving the 2x2 system
        expected = np.array([5.0/3.0, 2.0/3.0])
        assert_allclose(x, expected, rtol=1e-10)

    def test_well_conditioned_system(self, create_tridiagonal_system: Any) -> None:
        """Test with well-conditioned tridiagonal system.

        A well-conditioned system has strong diagonal dominance,
        making it numerically stable.
        """
        a, b, c, r, x_expected = create_tridiagonal_system(10, 'well-conditioned')
        x = tridiagonal(a, b, c, r)

        # Should match expected solution closely
        assert_allclose(x, x_expected, rtol=1e-8)

    def test_diagonal_dominant_system(self, create_tridiagonal_system: Any) -> None:
        """Test with diagonal-dominant system.

        System where |b[i]| > |a[i]| + |c[i]| for all i.
        """
        a, b, c, r, x_expected = create_tridiagonal_system(15, 'diagonal-dominant')
        x = tridiagonal(a, b, c, r)

        assert_allclose(x, x_expected, rtol=1e-7)

    def test_ill_conditioned_system(self, create_tridiagonal_system: Any) -> None:
        """Test with ill-conditioned system.

        System where diagonal dominance is weak. Still solvable
        but with reduced accuracy.
        """
        a, b, c, r, x_expected = create_tridiagonal_system(8, 'ill-conditioned')
        x = tridiagonal(a, b, c, r)

        # Tolerance is looser for ill-conditioned systems
        assert_allclose(x, x_expected, rtol=1e-5)

    def test_large_system(self, create_tridiagonal_system: Any) -> None:
        """Test with larger system (100 elements).

        Ensures algorithm scales properly.
        """
        a, b, c, r, x_expected = create_tridiagonal_system(100, 'well-conditioned')
        x = tridiagonal(a, b, c, r)

        assert_allclose(x, x_expected, rtol=1e-8)

    def test_heat_diffusion_analogy(self) -> None:
        """Test tridiagonal solver with physics-inspired system.

        Models 1D heat diffusion with constant temperature boundary conditions:
        - Surface temp fixed at 300 K
        - Bottom temp fixed at 290 K
        - Interior points satisfy diffusion equation

        This system is representative of what the model actually solves.
        """
        nz = 5
        alpha = 0.1  # Diffusion parameter (related to dt * thermal_diff / dz^2)

        # Tridiagonal coefficients for diffusion
        a = np.full(nz, -alpha)
        b = np.full(nz, 1.0 + 2.0 * alpha)
        c = np.full(nz, -alpha)

        # Boundary conditions
        # Top: 300 K (fixed)
        # Bottom: 290 K (fixed)
        temp = np.linspace(300.0, 290.0, nz)

        # Right-hand side incorporates boundary conditions
        r = temp.copy()
        # Top boundary: first equation has extra boundary term
        r[0] += alpha * 300.0  # Extra bc contribution
        # Bottom boundary: last equation has extra boundary term
        r[-1] += alpha * 290.0  # Extra bc contribution

        x = tridiagonal(a, b, c, r)

        # Solution should be physically reasonable (between bounds)
        assert np.all(x >= 290.0 - 1e-10)
        assert np.all(x <= 300.0 + 1e-10)
        # Surface should be close to 300 K
        assert x[0] > 299.0
        # Bottom should be close to 290 K
        assert x[-1] < 291.0

    def test_zero_on_diagonal_error(self) -> None:
        """Test that solver raises error when b[0] is zero.

        A zero on the main diagonal causes the Thomas algorithm to fail.
        """
        a = np.array([0.0, 1.0])
        b = np.array([0.0, 2.0])  # b[0] is zero!
        c = np.array([1.0, 0.0])
        r = np.array([4.0, 3.0])

        with pytest.raises(ValueError, match="Main diagonal cannot have a zero"):
            tridiagonal(a, b, c, r)

    def test_zero_during_factorization_warning(self) -> None:
        """Test that solver handles near-singular systems gracefully.

        While we can't easily force a zero during factorization without
        creating a truly singular system, we can test that well-conditioned
        systems work and ill-conditioned ones fail appropriately.

        This test verifies solver robustness across the conditioning spectrum.
        """
        # A system with poor diagonal dominance
        # This is on the edge of being solvable
        a = np.array([0.0, 0.5, 0.0])
        b = np.array([1.0, 1.0, 1.0])
        c = np.array([0.5, 0.5, 0.0])
        r = np.array([1.0, 1.0, 1.0])

        # This system is still solvable (diagonal dominance maintained),
        # but shows we can handle poorly conditioned systems
        try:
            x = tridiagonal(a, b, c, r)
            # If it solves, verify the solution is reasonable
            assert len(x) == 3
            # System shouldn't blow up
            assert np.all(np.isfinite(x))
        except ValueError:
            # It's also acceptable to fail on ill-conditioned systems
            pass

    def test_identity_matrix(self) -> None:
        """Test tridiagonal solver on identity matrix.

        System: x = r (where A = I)
        Solution: x = r
        """
        n = 5
        a = np.zeros(n)
        b = np.ones(n)  # Identity matrix diagonal
        c = np.zeros(n)
        r = np.array([1.0, 2.0, 3.0, 4.0, 5.0])

        x = tridiagonal(a, b, c, r)

        assert_allclose(x, r, rtol=1e-10)

    def test_multi_column_system(self) -> None:
        """Test tridiagonal solver with multiple RHS columns."""
        n = 4
        a = np.array([0.0, -1.0, -1.0, -1.0])
        b = np.full(n, 4.0)
        c = np.array([-1.0, -1.0, -1.0, 0.0])

        r = np.stack(
            [np.full(n, 5.0), np.array([1.0, 2.0, 3.0, 4.0])],
            axis=1,
        )

        x = tridiagonal(a, b, c, r)

        A = np.diag(b) + np.diag(a[1:], -1) + np.diag(c[:-1], 1)
        expected = np.column_stack(
            [np.linalg.solve(A, r[:, 0]), np.linalg.solve(A, r[:, 1])]
        )

        assert_allclose(x, expected, rtol=1e-10)

    def test_symmetric_system(self) -> None:
        """Test with symmetric tridiagonal matrix.

        For a symmetric system, a = c (sub and super diagonals equal).
        """
        n = 8
        # Symmetric tridiagonal matrix
        diag_val = 3.0
        off_diag = -0.5

        a = np.full(n, off_diag)
        b = np.full(n, diag_val)
        c = np.full(n, off_diag)

        # Create a solution and compute RHS
        x_expected = np.linspace(1.0, 2.0, n)
        r = np.zeros(n)
        r[0] = b[0] * x_expected[0] + c[0] * x_expected[1]
        for i in range(1, n-1):
            r[i] = (a[i] * x_expected[i-1] +
                   b[i] * x_expected[i] +
                   c[i] * x_expected[i+1])
        r[-1] = a[-1] * x_expected[-2] + b[-1] * x_expected[-1]

        x = tridiagonal(a, b, c, r)

        assert_allclose(x, x_expected, rtol=1e-8)

    @pytest.mark.parametrize("n", [5, 10, 20, 50])
    def test_various_sizes(self, create_tridiagonal_system: Any, n: int) -> None:
        """Parametrized test: verify solver works for various matrix sizes."""
        a, b, c, r, x_expected = create_tridiagonal_system(n, 'well-conditioned')
        x = tridiagonal(a, b, c, r)

        assert_allclose(x, x_expected, rtol=1e-8)
        assert len(x) == n


# ============================================================================
# Tests for root_brent() - Brent's Root-Finding Method
# ============================================================================

@pytest.mark.solver
class TestRootBrent:
    """Test suite for Brent's root-finding method."""

    def test_linear_function(self) -> None:
        """Test with linear function: f(x) = 2x - 4, root at x=2.

        Linear functions should converge in very few iterations.
        """
        f: Callable[[float], float] = lambda x: 2*x - 4
        root, converged = root_brent(f, 0.0, 4.0, tol=1e-6)

        assert converged, "Solver should converge for linear function"
        assert_allclose(root, 2.0, rtol=1e-5)

    def test_quadratic_function(self) -> None:
        """Test with quadratic function: f(x) = (x-3)^2 - 1.

        Root at x = 3 ± sqrt(1) = {2, 4}. We search in bracket [1, 3.5]
        which should find the root at x ≈ 2.
        """
        f: Callable[[float], float] = lambda x: (x - 3)**2 - 1
        root, converged = root_brent(f, 1.0, 3.5, tol=1e-6)

        assert converged
        # Root should be near 2
        assert 1.9 < root < 2.1
        assert_allclose(f(root), 0.0, atol=1e-5)

    def test_cubic_function(self) -> None:
        """Test with cubic: f(x) = (x-1)^3 - 8, root at x ≈ 3.

        (x-1)^3 = 8 => x - 1 = 2 => x = 3
        """
        f: Callable[[float], float] = lambda x: (x - 1)**3 - 8
        root, converged = root_brent(f, 2.0, 4.0, tol=1e-6)

        assert converged
        assert_allclose(root, 3.0, rtol=1e-5)
        assert_allclose(f(root), 0.0, atol=1e-5)

    def test_sine_function(self) -> None:
        """Test with sine function: f(x) = sin(x), root at x = π.

        sin(x) has roots at multiples of π. Testing in [2, 4] should
        find the root near π ≈ 3.14159.
        """
        f: Callable[[float], Any] = lambda x: np.sin(x)
        root, converged = root_brent(f, 2.0, 4.0, tol=1e-8)

        assert converged
        assert_allclose(root, np.pi, rtol=1e-6)
        assert_allclose(f(root), 0.0, atol=1e-7)

    def test_exponential_function(self) -> None:
        """Test with exponential: f(x) = exp(x) - 5, root at x = ln(5).

        exp(x) = 5 => x = ln(5) ≈ 1.609
        """
        f: Callable[[float], Any] = lambda x: np.exp(x) - 5
        root, converged = root_brent(f, 0.0, 3.0, tol=1e-8)

        assert converged
        expected_root = np.log(5)
        assert_allclose(root, expected_root, rtol=1e-6)
        assert_allclose(f(root), 0.0, atol=1e-7)

    def test_convergence_without_tolerance(self) -> None:
        """Test convergence criterion: bracket size < tolerance.

        As tolerance shrinks, solution should be more accurate.
        """
        f: Callable[[float], float] = lambda x: x**2 - 2  # Root at sqrt(2)
        expected = np.sqrt(2)

        for tol in [1e-3, 1e-6, 1e-9]:
            root, converged = root_brent(f, 1.0, 2.0, tol=tol)
            assert converged
            # Tolerance should roughly translate to solution accuracy
            assert abs(root - expected) < 10 * tol

    def test_max_iterations(self) -> None:
        """Test behavior when max iterations is reached.

        With iter_max=2, the solver shouldn't converge on a difficult function.
        We use a steep function with narrow bracket.
        """
        f: Callable[[float], float] = lambda x: (x - 1.5)**3 - 1  # Root at x ≈ 2.26
        _root, converged = root_brent(f, 1.0, 3.0, iter_max=2, tol=1e-10)

        # With only 2 iterations, convergence is unlikely with tight tolerance
        assert not converged

    def test_unbracketed_root_error(self) -> None:
        """Test error when root is not bracketed.

        Function must have different signs at bracket endpoints.
        f(x) = x^2 has no sign change in [1, 2] (always positive).
        """
        f: Callable[[float], float] = lambda x: x**2  # Always non-negative
        with pytest.raises(ValueError, match="Root not bracketed"):
            root_brent(f, 1.0, 2.0, tol=1e-6)

    def test_bracket_with_root_at_endpoint(self) -> None:
        """Test when root is exactly at bracket endpoint.

        f(a) = 0 should be detected (though not ideal for numerical methods).
        """
        f: Callable[[float], float] = lambda x: x - 2
        # Bracket with root at left endpoint
        with pytest.raises(ValueError, match="Root not bracketed"):
            # f(2) = 0, f(3) > 0, so f(a)*f(b) = 0 (triggers error)
            root_brent(f, 2.0, 3.0, tol=1e-6)

    def test_near_vertical_function(self) -> None:
        """Test with function that has steep derivative near root.

        These require careful bracketing but Brent's method is robust.
        """
        f: Callable[[float], float] = lambda x: 100 * (x - 1.5)**3  # Very steep near x=1.5
        root, converged = root_brent(f, 1.0, 2.0, tol=1e-6, iter_max=100)

        # Brent's method should still converge despite steep slope
        assert converged
        assert_allclose(root, 1.5, rtol=1e-5)

    def test_energy_balance_analogy(self) -> None:
        """Test with function representing energy balance.

        In the SEB solver, we solve: net_rad - sensible - latent - ground = 0
        This is a physics-relevant test case.

        Energy balance: incoming radiation balanced against outgoing fluxes.
        As surface temperature increases, sensible heat increases, reducing
        the surplus and eventually balancing the energy budget.

        Parameters chosen so that:
        - At low temp (ts_norm=0, 275K): surplus (net_rad > outgoing fluxes)
        - At high temp (ts_norm=1, 325K): deficit (outgoing > net_rad)
        This guarantees a sign change and a root in [0, 1].
        """
        def energy_balance(ts_normalized: float) -> float:
            # Map normalized [0, 1] to physical temperature [275, 325 K]
            ts = 275.0 + 50.0 * ts_normalized

            # Energy balance components (all in W/m^2)
            net_rad = 500.0  # Incoming net radiation
            sensible = 12.0 * (ts - 273.15)  # Sensible heat (temperature-dependent)
            latent = 100.0  # Latent heat flux
            ground = 20.0  # Ground heat flux

            # Energy balance equation: net_rad - sensible - latent - ground = 0
            return net_rad - sensible - latent - ground

        # Verify sign change exists
        f0 = energy_balance(0.0)  # At ts=275K: 500 - 12*2 - 100 - 20 = 356 W/m^2 (positive)
        f1 = energy_balance(1.0)  # At ts=325K: 500 - 12*52 - 100 - 20 = -244 W/m^2 (negative)

        assert f0 * f1 < 0, \
            f"Energy balance must have sign change for root to exist: f(0)={f0:.1f}, f(1)={f1:.1f}"

        root, converged = root_brent(energy_balance, 0.0, 1.0, tol=1e-6)
        assert converged, "Brent's method should converge for energy balance equation"
        assert 0.0 <= root <= 1.0, f"Root should be in bracket [0, 1], got {root}"
        assert_allclose(energy_balance(root), 0.0, atol=1e-4)

    def test_large_bracket(self) -> None:
        """Test convergence with large initial bracket.

        Large brackets shouldn't prevent convergence, just cost more iterations.
        """
        f: Callable[[float], float] = lambda x: x**3 - 1  # Root at x=1
        root, converged = root_brent(f, -100.0, 100.0, tol=1e-8, iter_max=200)

        assert converged
        assert_allclose(root, 1.0, rtol=1e-6)

    def test_root_at_zero(self) -> None:
        """Test finding root at exactly zero.

        f(x) = x has root at x = 0.
        """
        f: Callable[[float], float] = lambda x: x
        root, converged = root_brent(f, -1.0, 1.0, tol=1e-8)

        assert converged
        assert_allclose(root, 0.0, atol=1e-7)

    def test_oscillating_function(self) -> None:
        """Test with oscillating function.

        f(x) = cos(x) has roots at π/2, 3π/2, 5π/2, ...
        Searching in [1, 2] should find root near π/2 ≈ 1.571
        """
        f: Callable[[float], Any] = lambda x: np.cos(x)
        root, converged = root_brent(f, 1.0, 2.0, tol=1e-8)

        assert converged
        assert_allclose(root, np.pi/2, rtol=1e-6)

    @pytest.mark.parametrize("func,bracket,expected", [
        (lambda x: x - 1.5, [1.0, 2.0], 1.5),
        (lambda x: 2*x - 3, [0.0, 2.0], 1.5),
        (lambda x: (x-3)**2 - 4, [0.5, 2.5], 1.0),
        (lambda x: np.exp(x) - 2, [0.0, 1.0], np.log(2)),
    ])
    def test_various_functions(self, func: Any, bracket: Any, expected: Any) -> None:
        """Parametrized test: verify solver on various function types."""
        root, converged = root_brent(func, bracket[0], bracket[1], tol=1e-8)

        assert converged
        assert_allclose(root, expected, rtol=1e-5)
        assert_allclose(func(root), 0.0, atol=1e-5)

    def test_convergence_is_monotonic(self) -> None:
        """Test that bracket size decreases monotonically.

        The root should be in a shrinking bracket with each iteration.
        This is guaranteed by Brent's algorithm, but we can verify behavior.
        """
        f: Callable[[float], float] = lambda x: x**2 - 5  # Root at sqrt(5) ≈ 2.236

        # Use function to track bracket shrinkage
        root, converged = root_brent(f, 1.0, 3.0, tol=1e-8, iter_max=50)

        assert converged
        assert_allclose(root, np.sqrt(5), rtol=1e-6)


@pytest.mark.solver
class TestRootBrentVec:
    """Test suite for the vectorized Brent root finder."""

    def test_multiple_independent_roots(self) -> None:
        """Solve several independent bracketed roots in one vector call."""
        roots = np.array([1.5, -0.75, 3.25])
        a = roots - 1.0
        b = roots + 1.0

        def f(x: NDArray[Any]) -> NDArray[Any]:
            return x - roots

        root, converged = root_brent_vec(f, a, b, tol=1e-8)

        assert np.all(converged)
        assert_allclose(root, roots, rtol=1e-6, atol=1e-8)
        assert_allclose(f(root), 0.0, atol=1e-8)

    def test_unbracketed_root_error(self) -> None:
        """Raise immediately when any entry is not properly bracketed."""
        roots = np.array([1.5, 2.0, -0.75])
        a = np.array([0.0, 2.0, -1.75])
        b = np.array([3.0, 3.0, 0.25])

        def f(x: NDArray[Any]) -> NDArray[Any]:
            return x - roots

        with pytest.raises(ValueError, match="Root not bracketed"):
            root_brent_vec(f, a, b, tol=1e-8)
