"""Local refine of (A_V, v, sigma). x_j stays NNLS; never Metropolis."""

from __future__ import annotations

from typing import Callable, Tuple

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import minimize

EvalFn = Callable[
    [float, float, float],
    Tuple[NDArray[np.float64], NDArray[np.float64], float],
]


def local_axis_grid(
    center: float,
    coarse_step: float,
    low: float,
    high: float,
    factor: float = 5.0,
) -> NDArray[np.float64]:
    """Fine samples spanning about ± one coarse step around ``center``."""
    if coarse_step <= 0:
        return np.array([center], dtype=float)
    fine = coarse_step / factor
    n = int(round(factor))
    vals = center + fine * np.arange(-n, n + 1, dtype=float)
    vals = vals[(vals >= low - 1e-12) & (vals <= high + 1e-12)]
    if vals.size == 0:
        clipped = min(max(center, low), high)
        return np.array([clipped], dtype=float)
    return np.unique(np.round(vals, decimals=10))


def _clip_params(
    a_v: float,
    v0: float,
    sigma: float,
    a_v_bounds: Tuple[float, float],
    v_bounds: Tuple[float, float],
    sigma_bounds: Tuple[float, float],
    vary_kinematics: bool,
    v_fixed: float,
    s_fixed: float,
) -> Tuple[float, float, float]:
    a_v = float(min(max(a_v, a_v_bounds[0]), a_v_bounds[1]))
    if vary_kinematics:
        v0 = float(min(max(v0, v_bounds[0]), v_bounds[1]))
        sigma = float(min(max(sigma, sigma_bounds[0]), sigma_bounds[1]))
    else:
        v0 = v_fixed
        sigma = s_fixed
    return a_v, v0, sigma


def refine_av_v_sigma(
    eval_chi2: EvalFn,
    a_v: float,
    v0_kms: float,
    sigma_kms: float,
    *,
    a_v_step: float,
    v_step: float,
    sigma_step: float,
    a_v_bounds: Tuple[float, float],
    v_bounds: Tuple[float, float],
    sigma_bounds: Tuple[float, float],
    vary_kinematics: bool = True,
) -> Tuple[NDArray[np.float64], float, NDArray[np.float64], float, float, float]:
    """Coordinate-wise fine grids, then Nelder-Mead. Each trial is NNLS."""
    v_fixed = v0_kms
    s_fixed = sigma_kms
    x_hat, model, chi2 = eval_chi2(a_v, v0_kms, sigma_kms)
    best = [chi2, x_hat, a_v, model, v0_kms, sigma_kms]

    def consider(av: float, vv: float, ss: float) -> None:
        av, vv, ss = _clip_params(
            av, vv, ss, a_v_bounds, v_bounds, sigma_bounds, vary_kinematics, v_fixed, s_fixed
        )
        xx, mm, cc = eval_chi2(av, vv, ss)
        if cc < best[0]:
            best[0], best[1], best[2], best[3], best[4], best[5] = cc, xx, av, mm, vv, ss

    for _ in range(2):
        _, _, a_now, _, v_now, s_now = best
        if vary_kinematics:
            for vv in local_axis_grid(v_now, v_step, v_bounds[0], v_bounds[1]):
                consider(a_now, float(vv), s_now)
            _, _, a_now, _, v_now, s_now = best
            for ss in local_axis_grid(s_now, sigma_step, sigma_bounds[0], sigma_bounds[1]):
                consider(a_now, v_now, float(ss))
            _, _, a_now, _, v_now, s_now = best
        for av in local_axis_grid(a_now, a_v_step, a_v_bounds[0], a_v_bounds[1]):
            consider(float(av), v_now, s_now)

    _, _, a_now, _, v_now, s_now = best
    a_loc = local_axis_grid(a_now, a_v_step, a_v_bounds[0], a_v_bounds[1])
    if vary_kinematics:
        v_loc = local_axis_grid(v_now, v_step, v_bounds[0], v_bounds[1])
        s_loc = local_axis_grid(s_now, sigma_step, sigma_bounds[0], sigma_bounds[1])
        # Keep the 3-D neighbourhood cheap: every other fine point (±2 fine steps).
        a_loc = a_loc[::2]
        v_loc = v_loc[::2]
        s_loc = s_loc[::2]
        for vv in v_loc:
            for ss in s_loc:
                for av in a_loc:
                    consider(float(av), float(vv), float(ss))
    else:
        for av in a_loc:
            consider(float(av), v_now, s_now)

    _, _, a_now, _, v_now, s_now = best

    def objective(p: NDArray[np.float64]) -> float:
        av, vv, ss = _clip_params(
            float(p[0]),
            float(p[1]),
            float(p[2]),
            a_v_bounds,
            v_bounds,
            sigma_bounds,
            vary_kinematics,
            v_fixed,
            s_fixed,
        )
        _, _, cc = eval_chi2(av, vv, ss)
        return cc

    polished = minimize(
        objective,
        np.array([a_now, v_now, s_now], dtype=float),
        method="Nelder-Mead",
        options={"maxiter": 40, "xatol": 0.25, "fatol": 1e-8, "disp": False},
    )
    if polished.success or np.isfinite(polished.fun):
        av, vv, ss = _clip_params(
            float(polished.x[0]),
            float(polished.x[1]),
            float(polished.x[2]),
            a_v_bounds,
            v_bounds,
            sigma_bounds,
            vary_kinematics,
            v_fixed,
            s_fixed,
        )
        consider(av, vv, ss)

    chi2, x_hat, a_v, model, v0_kms, sigma_kms = best
    return x_hat, float(a_v), model, float(chi2), float(v0_kms), float(sigma_kms)
