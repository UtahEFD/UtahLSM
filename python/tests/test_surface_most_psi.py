"""Tests for selectable stable MOST ψ functions."""

import numpy as np

from utahlsm.physics.surface.sfc_most import SurfaceMOST


def test_default_stable_psi_matches_dyer_hicks():
    """Default stable ψ uses the legacy Dyer/Hicks linear form."""
    sfc = SurfaceMOST()
    zeta = np.array([0.0, 0.1, 0.5, 1.0, 5.0])
    expected = -5.0 * zeta
    actual = np.array([float(sfc.psim_stable(float(z))) for z in zeta])
    assert np.allclose(actual, expected)


def test_beljaars_holtslag_1991_is_well_behaved_for_stable_zeta():
    """BH91 stable ψ returns finite, non-positive values for zeta >= 0."""
    sfc = SurfaceMOST(psi_stable="beljaars-holtslag")
    zeta = np.array([0.0, 1e-6, 0.1, 1.0, 5.0, 20.0])

    psim = np.array([float(sfc.psim_stable(float(z))) for z in zeta])
    psih = np.array([float(sfc.psih_stable(float(z))) for z in zeta])

    assert np.isfinite(psim).all()
    assert np.isfinite(psih).all()
    assert (psim <= 1e-12).all()
    assert (psih <= 1e-12).all()
    assert psim[0] == 0.0
    assert psih[0] == 0.0


def test_cheng_brutsaert_2005_is_well_behaved_for_stable_zeta():
    """CB05 stable ψ returns finite, non-positive values for zeta >= 0."""
    sfc = SurfaceMOST(psi_stable="cheng-brutsaert")
    zeta = np.array([0.0, 1e-6, 0.1, 1.0, 5.0, 20.0])

    psim = np.array([float(sfc.psim_stable(float(z))) for z in zeta])
    psih = np.array([float(sfc.psih_stable(float(z))) for z in zeta])

    assert np.isfinite(psim).all()
    assert np.isfinite(psih).all()
    assert (psim <= 1e-12).all()
    assert (psih <= 1e-12).all()
    assert psim[0] == 0.0
    assert psih[0] == 0.0
