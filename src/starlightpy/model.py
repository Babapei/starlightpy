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


def reddening_columns(
    q_lambda: ArrayLike,
    a_v: float,
    q_lambda0: float,
    n_comp: int,
    a_yv: float = 0.0,
    young_flags: Optional[ArrayLike] = None,
) -> NDArray[np.float64]:
    """Per-template reddening. Young templates get A_V + A_YV when flagged."""
    q = np.asarray(q_lambda, dtype=float)
    if young_flags is None or abs(float(a_yv)) < 1e-12:
        return reddening_factor(q, a_v, q_lambda0)[:, None]
    young = np.asarray(young_flags, dtype=float)
    if young.shape != (n_comp,):
        raise ValueError("young_flags must have one value per template.")
    av_j = float(a_v) + float(a_yv) * young
    return 10.0 ** (-0.4 * (q[:, None] - q_lambda0) * av_j[None, :])


def build_model(
    x_j: ArrayLike,
    base_matrix: ArrayLike,
    extinction_curve: ArrayLike,
    a_v: float,
    q_lambda0: Optional[float] = None,
    a_yv: float = 0.0,
    young_flags: Optional[ArrayLike] = None,
) -> NDArray[np.float64]:
    """M_λ = Σ_j x_j b_{λ,j} r_{λ,j}(A_V, A_YV) (no kinematics)."""
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
    rmat = reddening_columns(q, a_v, q_lambda0, n_comp, a_yv=a_yv, young_flags=young_flags)
    return (bases * rmat) @ x
