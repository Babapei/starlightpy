"""Phase A optimizer: A_V grid + NNLS. Not STARLIGHT annealing (phase C)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import nnls

from .config import FitConfig
from .extinction import get_extinction_curve
from .kinematics import apply_losvd
from .model import build_model, normalize_at, normalize_bases, reddening_factor


@dataclass
class FitResult:
    x: NDArray[np.float64]
    x_fraction: NDArray[np.float64]
    a_v: float
    v0_kms: float
    sigma_kms: float
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
    config: Optional[FitConfig] = None,
) -> FitResult:
    """Fit x_j and A_V. Kinematics from config are applied as fixed values in phase A.

    Observation, error and templates must share the same wavelength grid.
    """
    if config is None:
        config = FitConfig()
    if config.search_kinematics:
        raise NotImplementedError("Kinematic search is PLAN phase B. Pass fixed v0_kms/sigma_kms.")
    if config.clip_nsigma is not None:
        raise NotImplementedError("Clipping is PLAN phase D.")

    wave = np.asarray(wavelengths, dtype=float)
    obs = np.asarray(flux, dtype=float)
    err = np.asarray(error, dtype=float)
    bases = np.asarray(base_matrix, dtype=float)
    n_wave, n_comp = bases.shape
    if obs.shape != (n_wave,) or err.shape != (n_wave,):
        raise ValueError("Observation, error and base wavelength axes must match.")

    good = np.ones(n_wave, dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
    obs_n, obs_scale = normalize_at(wave, obs, config.norm_window)
    bases_n = normalize_bases(wave, bases, config.norm_window)
    err_n = err / obs_scale

    q = get_extinction_curve(wave, law=config.law, r_v=config.r_v)
    win = (wave >= config.norm_window[0]) & (wave <= config.norm_window[1])
    q0 = float(np.median(q[win]))
    weights = _weights(err_n, good)

    a_v_grid = np.arange(config.a_v_bounds[0], config.a_v_bounds[1] + 0.5 * config.a_v_step, config.a_v_step)
    chi2_grid = np.empty_like(a_v_grid)
    x_grid = np.empty((a_v_grid.size, n_comp))

    for i, a_v in enumerate(a_v_grid):
        r = reddening_factor(q, float(a_v), q0)
        reddened = bases_n * r[:, None]
        design = np.empty_like(reddened)
        for j in range(n_comp):
            design[:, j] = apply_losvd(wave, reddened[:, j], config.v0_kms, config.sigma_kms)
        y = obs_n * weights
        x_hat, _ = nnls(design * weights[:, None], y)
        model = apply_losvd(wave, build_model(x_hat, bases_n, q, float(a_v), q_lambda0=q0), config.v0_kms, config.sigma_kms)
        chi2_grid[i] = float(np.sum(((obs_n - model) * weights) ** 2))
        x_grid[i] = x_hat

    best = int(np.argmin(chi2_grid))
    x_opt = x_grid[best]
    a_v_opt = float(a_v_grid[best])
    model_n = apply_losvd(
        wave,
        build_model(x_opt, bases_n, q, a_v_opt, q_lambda0=q0),
        config.v0_kms,
        config.sigma_kms,
    )
    n_good = int(np.sum(weights > 0))
    dof = max(n_good - n_comp - 1, 1)
    x_sum = float(np.sum(x_opt))
    x_frac = x_opt / x_sum if x_sum > 0 else x_opt
    return FitResult(
        x=x_opt,
        x_fraction=x_frac,
        a_v=a_v_opt,
        v0_kms=config.v0_kms,
        sigma_kms=config.sigma_kms,
        model=model_n * obs_scale,
        chi2=chi2_grid[best],
        chi2_reduced=chi2_grid[best] / dof,
        n_good=n_good,
        a_v_grid=a_v_grid,
        chi2_grid=chi2_grid,
    )
