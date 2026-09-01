"""Convenience wrapper around pPXF. This is not STARLIGHT and does not fit itself."""

from .fit import FitResult, fit_spectrum
from .io import load_sdss_fits, load_spectrum

__all__ = ["FitResult", "fit_spectrum", "load_sdss_fits", "load_spectrum"]
__version__ = "0.1.0"
