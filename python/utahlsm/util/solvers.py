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
"""A collection of numerical solver functions.

This module provides robust and efficient numerical solvers for common
mathematical problems encountered in the land-surface model, such as
solving systems of linear equations and finding roots of functions.
"""
from typing import Callable

import numpy as np
from numpy.typing import NDArray

from ..util.io import logging_helper

logger = logging_helper.get_logger('UTIL: Solvers')

# Tolerance for floating-point comparisons in numerical solvers
_SOLVER_TOL = 1e-12

def tridiagonal(
    a: NDArray[np.float64],
    b: NDArray[np.float64],
    c: NDArray[np.float64],
    r: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Solves a tridiagonal system of equations using the Thomas algorithm.

    This function efficiently solves the equation Ax = r, where A is a
    tridiagonal matrix defined by its sub-diagonal (a), main diagonal (b),
    and super-diagonal (c).

    Handles both 1D (single column) and 2D (multi-column) inputs via a
    unified code path. 1D inputs are internally promoted to 2D with ncol=1
    and the result is squeezed back before returning.

    Args:
        a: The sub-diagonal of the matrix (size n) or (n, ncol). a[0] is
            ignored.
        b: The main diagonal of the matrix (size n) or (n, ncol).
        c: The super-diagonal of the matrix (size n) or (n, ncol). c[n-1] is
            ignored.
        r: The right-hand side vector (size n) or (n, ncol).

    Returns:
        The solution vector u (size n) or (n, ncol).

    Raises:
        ValueError: If an element on the main diagonal is effectively zero
            during factorization.
    """
    r = np.asarray(r, dtype=np.float64)

    # Promote 1D to 2D so a single code path handles both cases
    squeeze = False
    if r.ndim == 1:
        r = r[:, None]
        squeeze = True
    elif r.ndim != 2:
        raise ValueError(
            f'Right-hand side has unsupported dimensions: {r.ndim}.')

    n, ncol = r.shape

    def _broadcast(vec: NDArray[np.float64]) -> NDArray[np.float64]:
        vec = np.asarray(vec, dtype=np.float64)
        if vec.ndim == 1:
            if vec.shape[0] != n:
                raise ValueError(
                    f'Tridiagonal coefficient has length {vec.shape[0]} '
                    f'but expected {n}.')
            return vec[:, None]
        if vec.ndim == 2:
            if vec.shape != (n, ncol):
                raise ValueError(
                    f'Tridiagonal coefficient has shape {vec.shape} but '
                    f'expected {(n, ncol)}.')
            return vec
        raise ValueError(
            f'Tridiagonal coefficient has unsupported dimensions: '
            f'{vec.ndim}.')

    a = _broadcast(a)
    b = _broadcast(b)
    c = _broadcast(c)

    u = np.zeros_like(r)
    gam = np.zeros_like(r)

    bet = b[0].copy()
    if np.any(np.abs(bet) < _SOLVER_TOL):
        logger.error(
            'Error in solve_tridiagonal: b[0] has entries effectively '
            'zero.')
        raise ValueError(
            'Main diagonal cannot have a zero on the first element.')

    u[0] = r[0] / bet

    for j in range(1, n):
        gam[j] = c[j-1] / bet
        bet = b[j] - a[j] * gam[j]
        if np.any(np.abs(bet) < _SOLVER_TOL):
            logger.error(
                'Error in solve_tridiagonal: effective zero on main '
                'diagonal at index %d.', j)
            raise ValueError(
                f'Effective zero on main diagonal at index {j} during '
                f'factorization.')
        u[j] = (r[j] - a[j] * u[j-1]) / bet

    for j in range(n-2, -1, -1):
        u[j] -= gam[j+1] * u[j+1]

    if squeeze:
        return u[:, 0]
    return u

def root_brent(f: Callable[[float], float], a: float, b: float,
    iter_max: int = 100, tol: float = 1e-6) -> tuple[float, bool]:
    """Finds the root of a function using Brent's method.

    This is a robust and fast root-finding algorithm that combines bisection,
    the secant method, and inverse quadratic interpolation. It is guaranteed
    to find a root if one exists within the given bracket.

    Args:
        f: The function for which to find a root, f(x) = 0.
        a: The lower bound of the bracket [a, b].
        b: The upper bound of the bracket [a, b].
        iter_max: The maximum number of iterations. Defaults to 100.
        tol: The desired tolerance for the root. Defaults to 1e-6.

    Returns:
        A tuple containing:
            - The approximate root of the function.
            - A boolean indicating whether the solver converged.

    Raises:
        ValueError: If the root is not bracketed (i.e., f(a) * f(b) >= 0).
    """
    fa = f(a)
    fb = f(b)

    if fa * fb >= 0:
        raise ValueError(
            'Root not bracketed in solve_root_brent (f(a) * f(b) >= 0).')

    # Ensure 'b' is the best current guess (the one with the function
    # value closer to zero)
    if abs(fa) < abs(fb):
        a, b = b, a
        fa, fb = fb, fa

    c, fc = a, fa  # c is the previous best approximation
    d: float = a   # d is the second to last best guess
    mflag = True   # mflag is true if the last step was a bisection

    for _ in range(iter_max):
        # Check for convergence: if the bracket is smaller than the tolerance
        if abs(b - a) < tol:
            return b, True

        # Use fast inverse quadratic interpolation if the three points
        # are distinct
        if (abs(fa) > tol and abs(fb) > tol and abs(fc) > tol and
                fa != fc and fb != fc):
            s = (a * fb * fc / ((fa - fb) * (fa - fc)) +
                 b * fa * fc / ((fb - fa) * (fb - fc)) +
                 c * fa * fb / ((fc - fa) * (fc - fb)))
        # Otherwise, fall back to the secant method
        else:
            s = b - fb * (b - a) / (fb - fa)

        # Condition 1: Is the new point outside the desired range?
        cond1 = (s < (3 * a + b) / 4.0) or (s > b)
        # Condition 2: Is the step not decreasing fast enough
        # (bisection was last step)?
        cond2 = mflag and (abs(s - b) >= abs(b - c) / 2.0)
        # Condition 3: Is the step not decreasing fast enough
        # (interpolation was last step)?
        cond3 = (not mflag) and (abs(s - b) >= abs(c - d) / 2.0)
        # Condition 4: Is the bracket shrinking too slowly
        # (bisection was last step)?
        cond4 = mflag and (abs(b - c) < tol)
        # Condition 5: Is the bracket shrinking too slowly
        # (interpolation was last step)?
        cond5 = (not mflag) and (abs(c - d) < tol)

        if cond1 or cond2 or cond3 or cond4 or cond5:
            # Fallback to bisection
            s = (a + b) / 2.0
            mflag = True
        else:
            mflag = False

        fs = f(s)
        d = c          # d is now the second to last best guess
        c, fc = b, fb  # The last best guess becomes the second to last

        # Move the bounds to keep the root bracketed
        if fa * fs < 0:
            b, fb = s, fs
        else:
            a, fa = s, fs

        # Ensure 'b' is always the best current root estimate
        if abs(fa) < abs(fb):
            a, b = b, a
            fa, fb = fb, fa

        # Check for convergence
        if abs(b - a) < tol:
            return b, True

    return b, False


def root_brent_vec(
    f: Callable[[NDArray[np.float64]], NDArray[np.float64]],
    a: NDArray[np.float64],
    b: NDArray[np.float64],
    iter_max: int = 100,
    tol: float = 1e-6
) -> tuple[NDArray[np.float64], NDArray[np.bool_]]:
    """Vectorized Brent's method for finding roots of multiple functions.

    This function efficiently solves f(x) = 0 for multiple independent problems
    simultaneously using array operations. It uses Brent's method which combines
    bisection, secant, and inverse quadratic interpolation.

    Based on the vectorization approach from https://github.com/adonath/array-brentq

    Args:
        f: A vectorized function that takes an array of x values and returns
            an array of function values f(x).
        a: Array of lower bracket bounds (size n).
        b: Array of upper bracket bounds (size n).
        iter_max: Maximum number of iterations. Defaults to 100.
        tol: Desired tolerance for convergence. Defaults to 1e-6.

    Returns:
        A tuple containing:
            - Array of approximate roots (size n).
            - Boolean array indicating convergence for each root.

    Raises:
        ValueError: If any root is not bracketed; i.e., if ``f(a) * f(b) >= 0``
            for one or more entries.

    Note:
        All problems are iterated together until all converge or iter_max is
        reached. This is efficient when problems have similar convergence rates.
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)

    # Handle scalar case
    if a.ndim == 0:
        a = np.atleast_1d(a)
        b = np.atleast_1d(b)
        scalar_input = True
    else:
        scalar_input = False

    n = len(a)

    fa = f(a)
    fb = f(b)

    # Check bracketing
    invalid = fa * fb >= 0
    if np.any(invalid):
        invalid_idx = np.where(invalid)[0]
        raise ValueError(
            'Root not bracketed in solve_root_brent_vec for '
            f'{invalid_idx.size} of {n} problems at indices '
            f'{invalid_idx.tolist()}.'
        )

    # Ensure b has the smaller function value (best guess)
    swap = np.abs(fa) < np.abs(fb)
    a, b = np.where(swap, b, a), np.where(swap, a, b)
    fa, fb = np.where(swap, fb, fa), np.where(swap, fa, fb)

    c = a.copy()
    fc = fa.copy()
    d = a.copy()

    mflag = np.ones(n, dtype=bool)
    converged = np.abs(b - a) < tol

    for _ in range(iter_max):
        if converged.all():
            break

        # Inverse quadratic interpolation conditions
        use_iqi = (
            (np.abs(fa) > tol) & (np.abs(fb) > tol) & (np.abs(fc) > tol)
            & (fa != fc) & (fb != fc) & ~converged
        )

        # Inverse quadratic interpolation
        # Protect denominators to avoid division by zero (np.where evaluates both branches)
        denom1 = (fa - fb) * (fa - fc)
        denom2 = (fb - fa) * (fb - fc)
        denom3 = (fc - fa) * (fc - fb)
        denom1 = np.where(
            np.abs(denom1) < _SOLVER_TOL,
            np.copysign(_SOLVER_TOL, denom1), denom1)
        denom2 = np.where(
            np.abs(denom2) < _SOLVER_TOL,
            np.copysign(_SOLVER_TOL, denom2), denom2)
        denom3 = np.where(
            np.abs(denom3) < _SOLVER_TOL,
            np.copysign(_SOLVER_TOL, denom3), denom3)

        s_iqi = np.where(
            use_iqi,
            (a * fb * fc / denom1 + b * fa * fc / denom2 + c * fa * fb / denom3),
            0.0
        )

        # Secant method (fallback)
        denom = fb - fa
        denom = np.where(
            np.abs(denom) < _SOLVER_TOL,
            np.copysign(_SOLVER_TOL, denom), denom)
        s_sec = b - fb * (b - a) / denom

        s = np.where(use_iqi, s_iqi, s_sec)

        # Conditions for falling back to bisection
        bound_lo = (3 * a + b) / 4.0
        bound_hi = b
        # Ensure bound_lo < bound_hi for the comparison
        bound_lo, bound_hi = np.minimum(bound_lo, bound_hi), np.maximum(bound_lo, bound_hi)

        cond1 = (s < bound_lo) | (s > bound_hi)
        cond2 = mflag & (np.abs(s - b) >= np.abs(b - c) / 2.0)
        cond3 = ~mflag & (np.abs(s - b) >= np.abs(c - d) / 2.0)
        cond4 = mflag & (np.abs(b - c) < tol)
        cond5 = ~mflag & (np.abs(c - d) < tol)

        use_bisect = cond1 | cond2 | cond3 | cond4 | cond5
        s = np.where(use_bisect & ~converged, (a + b) / 2.0, s)
        mflag = np.where(~converged, use_bisect, mflag)

        # Evaluate function at new point
        fs = f(s)

        # Update history
        d = np.where(~converged, c, d)
        c = np.where(~converged, b, c)
        fc = np.where(~converged, fb, fc)

        # Update bracket
        update_a = (fa * fs < 0) & ~converged
        b = np.where(update_a, s, b)
        fb = np.where(update_a, fs, fb)
        a = np.where(~update_a & ~converged, s, a)
        fa = np.where(~update_a & ~converged, fs, fa)

        # Ensure b is the best guess
        swap = (np.abs(fa) < np.abs(fb)) & ~converged
        a, b = np.where(swap, b, a), np.where(swap, a, b)
        fa, fb = np.where(swap, fb, fa), np.where(swap, fa, fb)

        # Check convergence
        converged |= np.abs(b - a) < tol

    result = b
    if scalar_input:
        return result[0], converged[0]
    return result, converged
