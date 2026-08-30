"""Prototype fitter: same physical model as STARLIGHT, different optimizer.

STARLIGHT (Cid Fernandes et al. 2005) minimises χ² with simulated annealing
plus Metropolis. This prototype instead grids A_V and solves the linear
mixture with non-negative least squares. Results are therefore *not*
numerically interchangeable with the Fortran binary.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence, Tuple

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import nnls

from .extinction import get_extinction_curve
from .model import build_model, normalize_at, reddening_factor


@dataclass
class FitResult:
    x: NDArray[np.float64]
    a_v: float
    model: NDArray[np.float64]
    chi2: float
    chi2_reduced: float
    n_good: int
    a_v_grid: NDArray[np.float64]
    chi2_grid: NDArray[np.float64]


def _weights(error: NDArray[np.float64], good: NDArray[np.bool_]) -> NDArray[np.float64]:
    w = np.zeros_like(error, dtype=float)
    usable = good & np.isfinite(error) & (error > 0)
    w[usable] = 1.0 / error[usable]
    return w


def fit_spectrum(
    wavelengths: ArrayLike,
    flux: ArrayLike,
    error: ArrayLike,
    base_matrix: ArrayLike,
    mask: Optional[ArrayLike] = None,
    *,
    law: str = "CCM",
    r_v: float = 3.1,
    a_v_bounds: Tuple[float, float] = (0.0, 3.0),
    a_v_step: float = 0.05,
    norm_window: Tuple[float, float] = (4010.0, 4060.0),
    yav_flags: Optional[Sequence[int]] = None,
) -> FitResult:
    """Fit luminosity fractions x_j and a global A_V.

    Parameters
    ----------
    wavelengths, flux, error
        Observed spectrum. ``error`` is used as 1/weight.
    base_matrix
        Template matrix (n_wave, n_components), already sampled on the
        same wavelength grid as the observation.
    mask
        Optional boolean array, True = keep pixel.
    """
    wave = np.asarray(wavelengths, dtype=float)
    obs = np.asarray(flux, dtype=float)
    err = np.asarray(error, dtype=float)
    bases = np.asarray(base_matrix, dtype=float)
    n_wave, n_comp = bases.shape
    if obs.shape != (n_wave,) or err.shape != (n_wave,):
        raise ValueError("Observation, error and base wavelength axes must match.")

    good = np.ones(n_wave, dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
    obs_n, obs_scale = normalize_at(wave, obs, norm_window)
    bases_n = np.empty_like(bases)
    for j in range(n_comp):
        bases_n[:, j], _ = normalize_at(wave, bases[:, j], norm_window)
    err_n = err / obs_scale

    q = get_extinction_curve(wave, law=law, r_v=r_v)
    q0 = float(np.median(q[(wave >= norm_window[0]) & (wave <= norm_window[1])]))
    weights = _weights(err_n, good)

    a_v_grid = np.arange(a_v_bounds[0], a_v_bounds[1] + 0.5 * a_v_step, a_v_step)
    chi2_grid = np.empty_like(a_v_grid)
    x_grid = np.empty((a_v_grid.size, n_comp))

    for i, a_v in enumerate(a_v_grid):
        r = reddening_factor(q, float(a_v), q0)
        design = (bases_n * r[:, None]) * weights[:, None]
        y = obs_n * weights
        x_hat, _ = nnls(design, y)
        residual = (obs_n - (bases_n * r[:, None]) @ x_hat) * weights
        chi2_grid[i] = float(np.sum(residual**2))
        x_grid[i] = x_hat

    best = int(np.argmin(chi2_grid))
    x_opt = x_grid[best]
    a_v_opt = float(a_v_grid[best])
    model_n = build_model(x_opt, bases_n, q, a_v_opt, q_lambda0=q0, yav_flags=yav_flags)
    n_good = int(np.sum(weights > 0))
    dof = max(n_good - n_comp - 1, 1)
    return FitResult(
        x=x_opt,
        a_v=a_v_opt,
        model=model_n * obs_scale,
        chi2=chi2_grid[best],
        chi2_reduced=chi2_grid[best] / dof,
        n_good=n_good,
        a_v_grid=a_v_grid,
        chi2_grid=chi2_grid,
    )
