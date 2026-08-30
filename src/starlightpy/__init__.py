"""STARLIGHT-style spectral synthesis prototype (not a Fortran clone)."""

from .extinction import get_extinction_curve
from .fit import FitResult, fit_spectrum
from .io import apply_mask, load_base_master, load_base_spectra, load_mask, load_spectrum
from .model import build_model, normalize_at

__all__ = [
    "FitResult",
    "apply_mask",
    "build_model",
    "fit_spectrum",
    "get_extinction_curve",
    "load_base_master",
    "load_base_spectra",
    "load_mask",
    "load_spectrum",
    "normalize_at",
]

__version__ = "0.1.0"
