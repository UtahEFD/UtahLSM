# 
# UtahLSM
# 
# Copyright (c) 2017–2025 Jeremy A. Gibbs
# Copyright (c) 2017–2025 Rob Stoll
# Copyright (c) 2017–2025 Eric Pardyjak
# Copyright (c) 2017–2025 Pete Willemsen
# 
# This file is part of UtahLSM.
# 
# This software is free and is distributed under the MIT License.
# See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
# 

import logging
import numpy as np

# Module-level logger
logger = logging.getLogger(__name__)

def tridiagonal(a: np.ndarray, b: np.ndarray, c: np.ndarray, r: np.ndarray) -> np.ndarray:
    """
    Solves a tridiagonal system of equations using the Thomas algorithm.
    Solves the equation Ax = r, where A is a tridiagonal matrix.

    :param a: The sub-diagonal of the matrix (size n). a[0] is ignored.
    :param b: The main diagonal of the matrix (size n).
    :param c: The super-diagonal of the matrix (size n). c[n-1] is ignored.
    :param r: The right-hand side vector (size n).
    :return: The solution vector u (size n).
    """
    n = len(b)
    u = np.zeros(n)
    gam = np.zeros(n)
    
    if b[0] == 0.0:
        logger.error("Error in solve_tridiagonal: b[0] is zero.")
        raise ValueError("Main diagonal cannot have a zero on the first element.")

    bet = b[0]
    u[0] = r[0] / bet

    for j in range(1, n):
        gam[j] = c[j-1] / bet
        bet = b[j] - a[j] * gam[j]
        if bet == 0.0:
            logger.error(f"Error in solve_tridiagonal: zero on main diagonal at index {j}.")
            raise ValueError(f"Zero on main diagonal at index {j} during factorization.")
        u[j] = (r[j] - a[j] * u[j-1]) / bet

    for j in range(n-2, -1, -1):
        u[j] -= gam[j+1] * u[j+1]
        
    return u

def root_brent(f, a, b, tol=1e-6, max_iter=100) -> float:
    """
    Finds the root of a function within a bracketed interval using Brent's method.
    This is a robust and fast hybrid of bisection, secant, and inverse quadratic interpolation.
    
    :param f: The function for which to find a root, f(x) = 0.
    :param a: The lower bound of the bracket [a, b].
    :param b: The upper bound of the bracket [a, b].
    :param tol: The desired tolerance for the root.
    :param max_iter: The maximum number of iterations to perform.
    :return: The root of the function.
    :raises ValueError: If the root is not bracketed.
    """
    fa = f(a)
    fb = f(b)
    
    if fa * fb >= 0:
        raise ValueError("Root not bracketed in solve_root_brent (f(a) * f(b) >= 0).")
    
    # Ensure 'b' is the best current guess (the one with the function value closer to zero)
    if abs(fa) < abs(fb):
        a, b = b, a
        fa, fb = fb, fa
    
    c, fc = a, fa  # c is the previous best approximation
    mflag = True   # mflag is true if the last step was a bisection
    
    for i in range(max_iter):
        # Check for convergence: if the bracket is smaller than the tolerance
        if abs(b - a) < tol:
            return b
    
        # Use fast inverse quadratic interpolation if the three points are distinct
        if abs(fa) > tol and abs(fb) > tol and abs(fc) > tol and fa != fc and fb != fc:
            s = (a * fb * fc / ((fa - fb) * (fa - fc)) +
                 b * fa * fc / ((fb - fa) * (fb - fc)) +
                 c * fa * fb / ((fc - fa) * (fc - fb)))
        # Otherwise, fall back to the secant method
        else:
            s = b - fb * (b - a) / (fb - fa)
             
        # This block is the core of Brent's method's robustness.
        # It decides whether to accept the interpolated point 's' or fall back to bisection.
        # Condition 1: Is the new point outside the desired range?
        cond1 = (s < (3 * a + b) / 4.0) or (s > b)
        # Condition 2: Is the step not decreasing fast enough (bisection was last step)?
        cond2 = mflag and (abs(s - b) >= abs(b - c) / 2.0)
        # Condition 3: Is the step not decreasing fast enough (interpolation was last step)?
        cond3 = (not mflag) and (abs(s - b) >= abs(c - d) / 2.0)
        # Condition 4: Is the bracket shrinking too slowly (bisection was last step)?
        cond4 = mflag and (abs(b - c) < tol)
        # Condition 5: Is the bracket shrinking too slowly (interpolation was last step)?
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
             