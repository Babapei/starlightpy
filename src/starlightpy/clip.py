"""NSIGMA residual clip and EX0-style dropping of tiny x_j (PLAN phase D)."""

from __future__ import annotations

from typing import Optional

import numpy as np
from numpy.typing import ArrayLike, NDArray


def clip_outliers(
    observed: ArrayLike,
    model: ArrayLike,
    error: ArrayLike,
    nsigma: float,
    good: Optional[ArrayLike] = None,
) -> NDArray[np.bool_]:
    """Pixels with |O-M| > nsigma * e among currently good pixels become False."""
    if nsigma <= 0:
        raise ValueError("clip nsigma must be positive.")
    obs = np.asarray(observed, dtype=float)
    mdl = np.asarray(model, dtype=float)
    err = np.asarray(error, dtype=float)
    if obs.shape != mdl.shape or obs.shape != err.shape:
        raise ValueError("observed, model and error must have the same shape.")
    if good is None:
        mask = np.ones(obs.shape, dtype=bool)
    else:
        mask = np.asarray(good, dtype=bool)
        if mask.shape != obs.shape:
            raise ValueError("good mask must match the observation shape.")
    out = mask.copy()
    usable = mask & np.isfinite(err) & (err > 0) & np.isfinite(obs) & np.isfinite(mdl)
    out[usable] = np.abs(obs[usable] - mdl[usable]) <= nsigma * err[usable]
    return out


def keep_components(x_fraction: ArrayLike, x_min_keep: float) -> NDArray[np.bool_]:
    """True for light fractions at or above the EX0 threshold. Always keeps at least one."""
    if x_min_keep < 0:
        raise ValueError("x_min_keep must be >= 0.")
    frac = np.asarray(x_fraction, dtype=float)
    keep = frac >= x_min_keep
    if not np.any(keep):
        keep = np.zeros(frac.size, dtype=bool)
        keep[int(np.argmax(frac))] = True
    return keep


def clip_and_refit(
    observed: ArrayLike,
    model: ArrayLike,
    error: ArrayLike,
    nsigma: float,
    good: Optional[ArrayLike] = None,
) -> NDArray[np.bool_]:
    """Return the good-pixel mask after NSIGMA clipping. ``fit_spectrum`` does the refit."""
    return clip_outliers(observed, model, error, nsigma, good)
