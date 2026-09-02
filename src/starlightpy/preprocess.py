"""Observation-frame helpers applied before fitting. Do not search redshift."""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .io import resample_to


def to_rest_frame(
    wavelength: ArrayLike,
    flux: ArrayLike,
    redshift: float,
    error: Optional[ArrayLike] = None,
):
    """Move F_λ from observed to rest frame. ``redshift`` is applied, never fitted."""
    if redshift < 0:
        raise ValueError("redshift must be >= 0.")
    factor = 1.0 + float(redshift)
    wave = np.asarray(wavelength, dtype=float) / factor
    y = np.asarray(flux, dtype=float) * factor
    if error is None:
        return wave, y
    err = np.asarray(error, dtype=float) * factor
    return wave, y, err


def vacuum_to_air(wavelength: ArrayLike) -> NDArray[np.float64]:
    """Morton / SDSS-style conversion; wavelength in Å."""
    vac = np.asarray(wavelength, dtype=float)
    if np.any(vac <= 0):
        raise ValueError("Wavelengths must be positive.")
    return vac / (1.0 + 2.735182e-4 + 131.4182 / vac**2 + 2.76249e8 / vac**4)


def air_to_vacuum(wavelength: ArrayLike) -> NDArray[np.float64]:
    """Approximate inverse of ``vacuum_to_air``; wavelength in Å."""
    air = np.asarray(wavelength, dtype=float)
    if np.any(air <= 0):
        raise ValueError("Wavelengths must be positive.")
    return air * (1.0 + 2.735182e-4 + 131.4182 / air**2 + 2.76249e8 / air**4)


def align_observation(
    template_wave: ArrayLike,
    flux: ArrayLike,
    error: ArrayLike,
    redshift: float = 0.0,
    wave_frame: str = "as_is",
) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Put observed F_λ onto the rest-frame vacuum template grid."""
    frame = wave_frame.lower()
    if frame not in ("as_is", "air", "vacuum"):
        raise ValueError("wave_frame must be 'as_is', 'air', or 'vacuum'.")
    target = np.asarray(template_wave, dtype=float)
    wave_data = air_to_vacuum(target) if frame == "air" else target.copy()
    if redshift == 0.0 and frame != "air":
        return np.asarray(flux, dtype=float), np.asarray(error, dtype=float)
    wave_rest, flux_rest, err_rest = to_rest_frame(wave_data, flux, redshift, error)
    return resample_to(wave_rest, flux_rest, target), resample_to(wave_rest, err_rest, target)
