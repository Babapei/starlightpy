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
    clip_nsigma: Optional[float] = None
    x_min_keep: float = 0.0
