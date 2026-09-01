"""Read ASCII spectra, plus an optional SDSS-like FITS loader.

Fitting still goes through pPXF (Cappellari). This module only reads files.
"""

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


def load_sdss_fits(
    filename: PathLike,
) -> Tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.bool_]]:
    """Read a 1D SDSS-like FITS spectrum for pPXF. Requires astropy.

    Independent copy of the starlightpy loader (no shared import). Fitting
    itself is pPXF — cite Cappellari, not this wrapper.
    """
    try:
        from astropy.io import fits
    except ImportError as exc:
        raise ImportError("load_sdss_fits requires astropy.") from exc

    path = Path(filename)
    with fits.open(path) as hdul:
        table_hdu = None
        for hdu in hdul:
            data = getattr(hdu, "data", None)
            columns = getattr(hdu, "columns", None)
            if data is None or columns is None:
                continue
            names = [c.name.lower() for c in columns]
            if "flux" in names:
                table_hdu = hdu
                break
        if table_hdu is None:
            raise ValueError("No FITS table HDU with a flux column.")
        data = table_hdu.data
        names = {n.lower(): n for n in data.dtype.names}
        flux = np.asarray(data[names["flux"]], dtype=float).reshape(-1)
        if "loglam" in names:
            wave = 10.0 ** np.asarray(data[names["loglam"]], dtype=float).reshape(-1)
        else:
            wave = None
            for key in ("wavelength", "lambda", "wave", "lam"):
                if key in names:
                    wave = np.asarray(data[names[key]], dtype=float).reshape(-1)
                    break
            if wave is None:
                raise ValueError("FITS table has flux but no wavelength/loglam column.")
        if "ivar" in names:
            ivar = np.asarray(data[names["ivar"]], dtype=float).reshape(-1)
            err = np.full_like(flux, np.inf)
            positive = ivar > 0
            err[positive] = 1.0 / np.sqrt(ivar[positive])
        elif "error" in names:
            err = np.asarray(data[names["error"]], dtype=float).reshape(-1)
        else:
            median = float(np.median(np.abs(flux[np.isfinite(flux)]))) if flux.size else 1.0
            err = np.full_like(flux, max(0.01 * median, 1e-8))
        good = np.isfinite(flux) & np.isfinite(err) & (err < np.inf)
        if "ivar" in names:
            good &= np.asarray(data[names["ivar"]], dtype=float).reshape(-1) > 0
        for key in ("and_mask", "mask", "or_mask"):
            if key in names:
                good &= np.asarray(data[names[key]]).reshape(-1) == 0
                break
        return wave, flux, err, good
