"""STARLIGHT-style spectral fitting (not a Fortran clone). See docs/PLAN.md."""

from .config import FitConfig
from .extinction import get_extinction_curve
from .fit import FitResult, fit_spectrum
from .io import apply_mask, load_base_master, load_base_spectra, load_mask, load_spectrum, resample_to
from .kinematics import apply_losvd
from .model import build_model, normalize_at
from .simulate import default_absorption_bases, mock_observation

__all__ = [
    "FitConfig",
    "FitResult",
    "apply_losvd",
    "apply_mask",
    "build_model",
    "fit_spectrum",
    "get_extinction_curve",
    "load_base_master",
    "load_base_spectra",
    "load_mask",
    "load_spectrum",
    "default_absorption_bases",
    "mock_observation",
    "normalize_at",
    "resample_to",
]

__version__ = "0.1.0"
