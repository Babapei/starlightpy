"""Read a simple ASCII spectrum: wavelength, flux, optional error."""

from __future__ import annotations

from pathlib import Path
from typing import Tuple, Union

import numpy as np
from numpy.typing import NDArray

PathLike = Union[str, Path]


def load_spectrum(
    filename: PathLike,
) -> Tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Load λ, flux, error from 2 or 3 column ASCII (comments start with #)."""
    wavelengths: list[float] = []
    fluxes: list[float] = []
    errors: list[float] = []

    with open(filename, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                raise ValueError(f"Invalid spectrum line: {line}")
            wavelengths.append(float(parts[0]))
            fluxes.append(float(parts[1]))
            errors.append(float(parts[2]) if len(parts) >= 3 else 1.0)

    return np.asarray(wavelengths), np.asarray(fluxes), np.asarray(errors)
