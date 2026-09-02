"""Line-of-sight Gaussian in velocity space (uniform log-λ), then resample back."""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.ndimage import gaussian_filter1d

C_KMS = 299792.458


def _ln_grid(
    ln: NDArray[np.float64],
    v0_kms: float,
    sigma_kms: float,
    pad: bool,
    oversample: int,
) -> NDArray[np.float64]:
    oversample = int(oversample)
    if oversample < 1:
        raise ValueError("losvd_oversample must be >= 1.")
    d_ln = float(ln[-1] - ln[0]) / (ln.size - 1) / oversample
    pad_ln = 0.0
    if pad:
        pad_ln = 5.0 * max(abs(v0_kms), sigma_kms, 40.0) / C_KMS
    return np.arange(ln[0] - pad_ln, ln[-1] + pad_ln + 0.5 * d_ln, d_ln)


def apply_losvd(
    wavelength: ArrayLike,
    flux: ArrayLike,
    v0_kms: float,
    sigma_kms: float,
    *,
    pad: bool = False,
    oversample: int = 1,
) -> NDArray[np.float64]:
    """Shift and broaden in velocity.

    ``v0_kms > 0`` is recession: absorption moves to longer wavelength
    (``λ → λ exp(v/c)``). Broadening is a Gaussian in ln(λ), not in Å pixels.
    ``pad`` extends the lnλ grid before the Gaussian so spectrum ends are less
    corrupted. ``oversample`` densifies the uniform lnλ grid.
    """
    wave = np.asarray(wavelength, dtype=float)
    y = np.asarray(flux, dtype=float)
    if wave.shape != y.shape:
        raise ValueError("wavelength and flux must have the same shape.")
    if abs(v0_kms) < 1e-9 and sigma_kms < 1e-9 and not pad and oversample == 1:
        return y.copy()
    if np.any(wave <= 0):
        raise ValueError("Wavelengths must be positive for log-lambda LOSVD.")
    if np.any(np.diff(wave) <= 0):
        raise ValueError("Wavelength array must be strictly increasing.")

    ln = np.log(wave)
    if not pad and oversample == 1:
        ln_u = np.linspace(ln[0], ln[-1], ln.size)
    else:
        ln_u = _ln_grid(ln, v0_kms, sigma_kms, pad, oversample)
    y_u = np.interp(ln_u, ln, y)
    y_u = np.interp(ln_u - v0_kms / C_KMS, ln_u, y_u)

    if sigma_kms >= 1e-9:
        d_ln = float(ln_u[1] - ln_u[0])
        sigma_pix = sigma_kms / C_KMS / d_ln
        if sigma_pix >= 0.05:
            y_u = gaussian_filter1d(y_u, sigma=sigma_pix, mode="nearest")

    return np.interp(ln, ln_u, y_u)


def losvd_design(
    wavelength: ArrayLike,
    bases: ArrayLike,
    v0_kms: float,
    sigma_kms: float,
    *,
    pad: bool = False,
    oversample: int = 1,
) -> NDArray[np.float64]:
    """Apply LOSVD to each template column. Used to reuse one (v, σ) design."""
    wave = np.asarray(wavelength, dtype=float)
    matrix = np.asarray(bases, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != wave.size:
        raise ValueError("bases must be (n_wave, n_comp) on the wavelength grid.")
    out = np.empty_like(matrix, dtype=float)
    for j in range(matrix.shape[1]):
        out[:, j] = apply_losvd(
            wave, matrix[:, j], v0_kms, sigma_kms, pad=pad, oversample=oversample
        )
    return out


class LosvdDesignCache:
    """Fingerprint-keyed cache for a single (v, σ) design matrix."""

    def __init__(self) -> None:
        self._store = {}

    def _key(
        self,
        bases: NDArray[np.float64],
        v0_kms: float,
        sigma_kms: float,
        pad: bool,
        oversample: int,
    ) -> Tuple:
        return (
            round(float(v0_kms), 6),
            round(float(sigma_kms), 6),
            bool(pad),
            int(oversample),
            bases.shape,
            round(float(bases[0, 0]), 12),
            round(float(bases[-1, -1]), 12),
            round(float(np.sum(bases)), 8),
        )

    def design(
        self,
        wavelength: ArrayLike,
        bases: ArrayLike,
        v0_kms: float,
        sigma_kms: float,
        *,
        pad: bool = False,
        oversample: int = 1,
    ) -> NDArray[np.float64]:
        matrix = np.asarray(bases, dtype=float)
        key = self._key(matrix, v0_kms, sigma_kms, pad, oversample)
        hit = self._store.get(key)
        if hit is not None:
            return hit
        computed = losvd_design(
            wavelength, matrix, v0_kms, sigma_kms, pad=pad, oversample=oversample
        )
        self._store[key] = computed
        return computed
