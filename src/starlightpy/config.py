from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple


@dataclass
class FitConfig:
    """Whitelist of options. New fields must be added in docs/PLAN.md first."""

    norm_window: Tuple[float, float] = (4010.0, 4060.0)
    law: str = "CCM"
    r_v: float = 3.1
    a_v_bounds: Tuple[float, float] = (0.0, 3.0)
    a_v_step: float = 0.05
    v0_kms: float = 0.0
    sigma_kms: float = 0.0
    search_kinematics: bool = False
    v_bounds: Tuple[float, float] = (-300.0, 300.0)
    v_step: float = 50.0
    sigma_bounds: Tuple[float, float] = (40.0, 250.0)
    sigma_step: float = 50.0
    clip_nsigma: Optional[float] = None
    x_min_keep: float = 0.0
    redshift: float = 0.0
    wave_frame: str = "as_is"
    fwhm_data: Optional[float] = None
    fwhm_template: Optional[float] = None
    estimate_error: bool = False
    refine_kinematics: bool = False
    regularize_x: Optional[str] = None
    regularize_strength: float = 1.0
    age_bin_edges: Optional[Tuple[float, ...]] = None
    error_method: Optional[str] = None
    n_repeat: int = 8
    repeat_seed: int = 0
    pad_losvd: bool = False
    losvd_oversample: int = 1
    fit_ayv: bool = False
    a_yv_bounds: Tuple[float, float] = (0.0, 1.5)
    a_yv_step: float = 0.2
