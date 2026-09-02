"""Coarse uncertainties. Profile χ² slices or repeat fits; not covariances."""

from __future__ import annotations

from typing import Callable, Dict, Optional

import numpy as np
from numpy.typing import NDArray


def halfwidth_chi2_plus_one(
    theta: NDArray[np.float64],
    chi2: NDArray[np.float64],
    chi2_min: float,
    fallback: float,
) -> float:
    """Half-width of {θ : χ²(θ) ≤ χ²_min + 1}. Not a Bayesian interval."""
    ok = np.asarray(chi2) <= float(chi2_min) + 1.0 + 1e-9
    th = np.asarray(theta, dtype=float)
    if th.size == 0 or not np.any(ok):
        return float(abs(fallback))
    lo = float(th[ok].min())
    hi = float(th[ok].max())
    if hi <= lo:
        return 0.5 * float(abs(fallback))
    return 0.5 * (hi - lo)


def scan_slice(
    lo: float,
    hi: float,
    n: int,
    chi2_at: Callable[[float], float],
) -> tuple:
    if hi < lo:
        lo, hi = hi, lo
    if n < 3:
        n = 3
    grid = np.linspace(lo, hi, n)
    chi2 = np.empty(n)
    for i, t in enumerate(grid):
        chi2[i] = float(chi2_at(t))
    return grid, chi2


def x_fraction_std(samples: NDArray[np.float64]) -> NDArray[np.float64]:
    if samples.ndim != 2 or samples.shape[0] < 2:
        raise ValueError("Need at least two repeat fits to estimate x scatter.")
    return np.std(samples, axis=0, ddof=1)
