"""Linear mix of reddened SSP templates. No interpolation, no LOSVD."""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
from numpy.typing import ArrayLike, NDArray


def normalize_at(
    wavelengths: ArrayLike,
    flux: ArrayLike,
    window: Tuple[float, float] = (4010.0, 4060.0),
) -> Tuple[NDArray[np.float64], float]:
    wave = np.asarray(wavelengths, dtype=float)
    y = np.asarray(flux, dtype=float)
    mask = (wave >= window[0]) & (wave <= window[1])
    if not np.any(mask):
        raise ValueError(f"Normalization window {window} does not overlap the wavelength array.")
    scale = float(np.median(y[mask]))
    if scale == 0:
        raise ValueError("Normalization scale is zero.")
    return y / scale, scale


def normalize_bases(wavelengths: ArrayLike, bases: ArrayLike, window: Tuple[float, float]) -> NDArray[np.float64]:
    wave = np.asarray(wavelengths, dtype=float)
    matrix = np.asarray(bases, dtype=float)
    out = np.empty_like(matrix)
    for j in range(matrix.shape[1]):
        out[:, j], _ = normalize_at(wave, matrix[:, j], window)
    return out


def reddening_factor(q_lambda: ArrayLike, a_v: float, q_lambda0: float) -> NDArray[np.float64]:
    q = np.asarray(q_lambda, dtype=float)
    return 10.0 ** (-0.4 * (q - q_lambda0) * a_v)


def build_model(
    x_j: ArrayLike,
    base_matrix: ArrayLike,
    extinction_curve: ArrayLike,
    a_v: float,
    q_lambda0: Optional[float] = None,
) -> NDArray[np.float64]:
    """M_λ = Σ_j x_j b_{λ,j} 10^{-0.4 (q_λ - q_λ0) A_V} (no kinematics)."""
    bases = np.asarray(base_matrix, dtype=float)
    x = np.asarray(x_j, dtype=float)
    q = np.asarray(extinction_curve, dtype=float)
    n_wave, n_comp = bases.shape
    if x.shape != (n_comp,):
        raise ValueError(f"x_j has shape {x.shape}, expected ({n_comp},)")
    if q.shape != (n_wave,):
        raise ValueError("extinction_curve must match the wavelength axis of base_matrix.")
    if q_lambda0 is None:
        q_lambda0 = float(np.median(q))
    factor = reddening_factor(q, a_v, q_lambda0)
    return (bases * factor[:, None]) @ x
