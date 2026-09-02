"""Post-fit products. Refuse mass or age/Z summaries without metadata."""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
from numpy.typing import ArrayLike, NDArray


def light_to_mass(
    x: ArrayLike,
    mass_to_light: Optional[ArrayLike],
) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Convert light weights ``x`` to mass weights and mass fractions.

    ``mass_to_light`` is Υ_j = M/L at the same norm window as ``x``.
    Missing M/L is not invented.
    """
    if mass_to_light is None:
        raise ValueError("mass_to_light is required; refuse to invent M/L.")
    x_arr = np.asarray(x, dtype=float)
    ml = np.asarray(mass_to_light, dtype=float)
    if ml.shape != x_arr.shape:
        raise ValueError("mass_to_light must match x.")
    if np.any(ml < 0) or not np.all(np.isfinite(ml)):
        raise ValueError("mass_to_light must be finite and non-negative.")
    if not np.all(np.isfinite(x_arr)):
        raise ValueError("x must be finite.")
    mass = x_arr * ml
    total = float(np.sum(mass))
    if total <= 0:
        raise ValueError("Total mass weight is zero.")
    return mass, mass / total


def light_weighted_mean(x: ArrayLike, values: Optional[ArrayLike]) -> float:
    """Light-weighted mean of per-template metadata (age, Z, ...)."""
    if values is None:
        raise ValueError("values metadata is required.")
    x_arr = np.asarray(x, dtype=float)
    v = np.asarray(values, dtype=float)
    if v.shape != x_arr.shape:
        raise ValueError("values must match x.")
    if not np.all(np.isfinite(v)):
        raise ValueError("values must be finite.")
    if not np.all(np.isfinite(x_arr)):
        raise ValueError("x must be finite.")
    wsum = float(np.sum(x_arr))
    if wsum <= 0:
        raise ValueError("Total light weight is zero.")
    return float(np.sum(x_arr * v) / wsum)


def light_weighted_age(x: ArrayLike, ages: Optional[ArrayLike]) -> float:
    if ages is None:
        raise ValueError("ages metadata is required.")
    return light_weighted_mean(x, ages)


def light_weighted_metallicity(x: ArrayLike, metallicities: Optional[ArrayLike]) -> float:
    if metallicities is None:
        raise ValueError("metallicities metadata is required.")
    return light_weighted_mean(x, metallicities)
