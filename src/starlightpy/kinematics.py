"""Line-of-sight Gaussian in velocity, applied on a log-wavelength grid."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.ndimage import gaussian_filter1d

C_KMS = 299792.458


def apply_losvd(
    wavelength: ArrayLike,
    flux: ArrayLike,
    v0_kms: float,
    sigma_kms: float,
) -> NDArray[np.float64]:
    """Shift by v0 and broaden by sigma (km/s). Identity if both are ~0."""
    wave = np.asarray(wavelength, dtype=float)
    y = np.asarray(flux, dtype=float)
    if abs(v0_kms) < 1e-9 and sigma_kms < 1e-9:
        return y.copy()
    if np.any(wave <= 0):
        raise ValueError("Wavelengths must be positive for log-lambda LOSVD.")

    ln = np.log(wave)
    ln_out = ln.copy()
    ln_shifted = ln - v0_kms / C_KMS
    y_shift = np.interp(ln_out, ln_shifted, y)

    if sigma_kms < 1e-9:
        return y_shift

    d_ln = float(np.median(np.diff(ln)))
    if d_ln <= 0:
        raise ValueError("Wavelength array must be strictly increasing.")
    sigma_pix = sigma_kms / C_KMS / d_ln
    if sigma_pix < 0.05:
        return y_shift
    return gaussian_filter1d(y_shift, sigma=sigma_pix, mode="nearest")
