"""STARLIGHT-style spectral fitting (not a Fortran clone). See docs/PLAN.md."""

from .clip import clip_and_refit, clip_outliers, keep_components
from .config import FitConfig
from .extinction import get_extinction_curve
from .fit import FitResult, fit_spectrum
from .io import (
    apply_mask,
    combine_good,
    good_from_flags,
    iter_ascii_spectra,
    load_base_master,
    load_base_spectra,
    load_cxt,
    load_mask,
    load_sdss_fits,
    load_spectrum,
    resample_to,
)
from .kinematics import apply_losvd
from .model import build_model, normalize_at
from .preprocess import (
    air_to_vacuum,
    align_observation,
    match_instrumental_fwhm,
    optical_emission_mask_regions,
    to_rest_frame,
    vacuum_to_air,
)
from .simulate import default_absorption_bases, mock_observation

__all__ = [
    "FitConfig",
    "FitResult",
    "apply_losvd",
    "apply_mask",
    "build_model",
    "fit_spectrum",
    "get_extinction_curve",
    "combine_good",
    "good_from_flags",
    "iter_ascii_spectra",
    "load_base_master",
    "load_base_spectra",
    "load_cxt",
    "load_mask",
    "load_sdss_fits",
    "load_spectrum",
    "clip_and_refit",
    "clip_outliers",
    "keep_components",
    "default_absorption_bases",
    "mock_observation",
    "air_to_vacuum",
    "align_observation",
    "match_instrumental_fwhm",
    "optical_emission_mask_regions",
    "to_rest_frame",
    "vacuum_to_air",
    "normalize_at",
    "resample_to",
]

__version__ = "0.1.0"
