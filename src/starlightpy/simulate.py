"""Build mock observations with the same forward model as the fitter.

Used for tests (PLAN B1). Templates must have absorption lines if the mock
will later be used to recover velocity dispersion.
"""

from __future__ import annotations

from typing import Optional, Sequence, Tuple

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .config import FitConfig
from .extinction import get_extinction_curve
from .kinematics import apply_losvd
from .model import build_model, normalize_at, normalize_bases


def absorption_template(
    wavelength: ArrayLike,
    line_centers: Sequence[float],
    line_width: float = 6.0,
    line_depth: float = 0.45,
    continuum: float = 1.0,
) -> NDArray[np.float64]:
    """Smooth continuum minus Gaussian absorption. Lines should sit outside the fit norm window."""
    wave = np.asarray(wavelength, dtype=float)
    spec = np.full_like(wave, continuum, dtype=float)
    for center in line_centers:
        spec -= line_depth * np.exp(-0.5 * ((wave - center) / line_width) ** 2)
    return np.clip(spec, 0.15, None)


def default_absorption_bases(wavelength: ArrayLike) -> NDArray[np.float64]:
    """Three distinguishable bases; lines avoid 4010–4060 Å so normalization stays on continuum."""
    wave = np.asarray(wavelength, dtype=float)
    return np.column_stack(
        [
            absorption_template(wave, [4861.0, 5172.0]),
            absorption_template(wave, [4920.0, 5270.0], line_depth=0.35),
            absorption_template(wave, [4959.0, 5320.0], line_width=8.0, line_depth=0.40),
        ]
    )


def mock_observation(
    wavelength: ArrayLike,
    base_matrix: ArrayLike,
    x_j: ArrayLike,
    a_v: float,
    v0_kms: float = 0.0,
    sigma_kms: float = 0.0,
    *,
    config: Optional[FitConfig] = None,
    snr: Optional[float] = None,
    rng: Optional[np.random.Generator] = None,
    a_yv: float = 0.0,
    young_flags: Optional[ArrayLike] = None,
) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return flux and error on the same grid.

    Scaling matches ``fit_spectrum``: templates and the noise-free model are
    normalized in ``config.norm_window``, then LOSVD is applied. ``x_j`` should
    be light fractions (non-negative). If ``snr`` is None, error is 1% of the
    median flux and no noise is added.
    """
    if config is None:
        config = FitConfig()
    wave = np.asarray(wavelength, dtype=float)
    bases = np.asarray(base_matrix, dtype=float)
    x = np.asarray(x_j, dtype=float)
    bases_n = normalize_bases(wave, bases, config.norm_window)
    q = get_extinction_curve(wave, law=config.law, r_v=config.r_v)
    win = (wave >= config.norm_window[0]) & (wave <= config.norm_window[1])
    if not np.any(win):
        raise ValueError("norm_window does not overlap the wavelength array.")
    q0 = float(np.median(q[win]))
    mixed = build_model(
        x, bases_n, q, a_v, q_lambda0=q0, a_yv=a_yv, young_flags=young_flags
    )
    model = apply_losvd(wave, mixed, v0_kms, sigma_kms)
    _, model_scale = normalize_at(wave, model, config.norm_window)
    flux = model / model_scale
    median = float(np.median(np.abs(flux[win])))
    if snr is None:
        err = np.full_like(flux, max(0.01 * median, 1e-8))
        return flux, err
    noise_std = median / float(snr)
    if rng is None:
        rng = np.random.default_rng()
    flux = flux + rng.normal(0.0, noise_std, size=flux.shape)
    err = np.full_like(flux, max(noise_std, 1e-8))
    return flux, err
