"""Line-of-sight Gaussian in velocity space (uniform log-λ), then resample back."""

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
    """Shift and broaden in velocity.

    ``v0_kms > 0`` is recession: absorption moves to longer wavelength
    (``λ → λ exp(v/c)``). Broadening is a Gaussian in ln(λ), not in Å pixels.
    """
    wave = np.asarray(wavelength, dtype=float)
    y = np.asarray(flux, dtype=float)
    if wave.shape != y.shape:
        raise ValueError("wavelength and flux must have the same shape.")
    if abs(v0_kms) < 1e-9 and sigma_kms < 1e-9:
        return y.copy()
    if np.any(wave <= 0):
        raise ValueError("Wavelengths must be positive for log-lambda LOSVD.")
    if np.any(np.diff(wave) <= 0):
        raise ValueError("Wavelength array must be strictly increasing.")

    ln = np.log(wave)
    n = wave.size
    ln_u = np.linspace(ln[0], ln[-1], n)
    y_u = np.interp(ln_u, ln, y)
    y_u = np.interp(ln_u - v0_kms / C_KMS, ln_u, y_u)

    if sigma_kms >= 1e-9:
        d_ln = float(ln_u[1] - ln_u[0])
        sigma_pix = sigma_kms / C_KMS / d_ln
        if sigma_pix >= 0.05:
            y_u = gaussian_filter1d(y_u, sigma=sigma_pix, mode="nearest")

    return np.interp(ln, ln_u, y_u)
