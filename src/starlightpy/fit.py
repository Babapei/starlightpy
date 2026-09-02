"""A_V grid + NNLS; optional coarse grid over (v, sigma). x_j is never annealed."""

from __future__ import annotations

import warnings
from dataclasses import dataclass, replace
from typing import Optional, Tuple, Dict

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import nnls

from .clip import clip_outliers, keep_components
from .config import FitConfig
from .extinction import get_extinction_curve
from .kinematics import apply_losvd
from .model import build_model, normalize_at, normalize_bases, reddening_factor
from .preprocess import align_observation, estimate_rms_error, match_instrumental_fwhm
from .errors import halfwidth_chi2_plus_one, scan_slice, x_fraction_std
from .refine import refine_av_v_sigma
from .regularize import XRegularizer


@dataclass
class FitResult:
    x: NDArray[np.float64]
    x_fraction: NDArray[np.float64]
    a_v: float
    v0_kms: float
    sigma_kms: float
    model: NDArray[np.float64]
    chi2: float
    chi2_reduced: float
    n_good: int
    a_v_grid: NDArray[np.float64]
    chi2_grid: NDArray[np.float64]
    n_clipped: int = 0
    good: Optional[NDArray[np.bool_]] = None
    obs_scale: float = 1.0
    config: Optional[FitConfig] = None
    dropped: Optional[NDArray[np.bool_]] = None
    errors: Optional[dict] = None


def _weights(error: NDArray[np.float64], good: NDArray[np.bool_]) -> NDArray[np.float64]:
    w = np.zeros_like(error, dtype=float)
    usable = good & np.isfinite(error) & (error > 0)
    w[usable] = 1.0 / error[usable]
    return w


def _closed_grid(low: float, high: float, step: float) -> NDArray[np.float64]:
    if step <= 0:
        raise ValueError("Grid step must be positive.")
    if high < low:
        raise ValueError("Grid high must be >= low.")
    return np.arange(low, high + 0.5 * step, step)


def _design_and_nnls(
    wave: NDArray[np.float64],
    obs_n: NDArray[np.float64],
    weights: NDArray[np.float64],
    bases_n: NDArray[np.float64],
    q: NDArray[np.float64],
    q0: float,
    a_v: float,
    v0_kms: float,
    sigma_kms: float,
    x_reg: Optional[XRegularizer] = None,
) -> Tuple[NDArray[np.float64], NDArray[np.float64], float]:
    r = reddening_factor(q, a_v, q0)
    reddened = bases_n * r[:, None]
    n_comp = bases_n.shape[1]
    design = np.empty_like(reddened)
    for j in range(n_comp):
        design[:, j] = apply_losvd(wave, reddened[:, j], v0_kms, sigma_kms)
    weighted_design = design * weights[:, None]
    weighted_obs = obs_n * weights
    if x_reg is None:
        x_hat, _ = nnls(weighted_design, weighted_obs)
    else:
        x_hat = x_reg.solve(weighted_design, weighted_obs)
    model = apply_losvd(wave, build_model(x_hat, bases_n, q, a_v, q_lambda0=q0), v0_kms, sigma_kms)
    chi2 = float(np.sum(((obs_n - model) * weights) ** 2))
    return x_hat, model, chi2


def _best_av_for_kinematics(
    wave: NDArray[np.float64],
    obs_n: NDArray[np.float64],
    weights: NDArray[np.float64],
    bases_n: NDArray[np.float64],
    q: NDArray[np.float64],
    q0: float,
    a_v_grid: NDArray[np.float64],
    v0_kms: float,
    sigma_kms: float,
    x_reg: Optional[XRegularizer] = None,
) -> Tuple[NDArray[np.float64], float, NDArray[np.float64], float, NDArray[np.float64]]:
    n_comp = bases_n.shape[1]
    chi2_grid = np.empty(a_v_grid.size)
    x_grid = np.empty((a_v_grid.size, n_comp))
    models = np.empty((a_v_grid.size, wave.size))
    for i, a_v in enumerate(a_v_grid):
        x_hat, model, chi2 = _design_and_nnls(
            wave, obs_n, weights, bases_n, q, q0, float(a_v), v0_kms, sigma_kms, x_reg=x_reg
        )
        x_grid[i] = x_hat
        models[i] = model
        chi2_grid[i] = chi2
    best = int(np.argmin(chi2_grid))
    return x_grid[best], float(a_v_grid[best]), models[best], float(chi2_grid[best]), chi2_grid


def _errors_chi2_slice(
    config: FitConfig,
    wave: NDArray[np.float64],
    obs_n: NDArray[np.float64],
    weights: NDArray[np.float64],
    bases_n: NDArray[np.float64],
    q: NDArray[np.float64],
    q0: float,
    a_v_opt: float,
    v_opt: float,
    s_opt: float,
    chi2: float,
    x_reg: Optional[XRegularizer],
) -> dict:
    def chi2_at_av(av: float) -> float:
        _, _, c = _design_and_nnls(
            wave, obs_n, weights, bases_n, q, q0, av, v_opt, s_opt, x_reg=x_reg
        )
        return c

    lo = max(config.a_v_bounds[0], a_v_opt - 2.0 * config.a_v_step)
    hi = min(config.a_v_bounds[1], a_v_opt + 2.0 * config.a_v_step)
    grid, chi2s = scan_slice(lo, hi, 21, chi2_at_av)
    out: Dict[str, float] = {
        "a_v": halfwidth_chi2_plus_one(grid, chi2s, chi2, config.a_v_step),
    }
    if config.search_kinematics:
        def chi2_at_v(vv: float) -> float:
            _, _, c = _design_and_nnls(
                wave, obs_n, weights, bases_n, q, q0, a_v_opt, vv, s_opt, x_reg=x_reg
            )
            return c

        def chi2_at_s(ss: float) -> float:
            _, _, c = _design_and_nnls(
                wave, obs_n, weights, bases_n, q, q0, a_v_opt, v_opt, ss, x_reg=x_reg
            )
            return c

        vlo = max(config.v_bounds[0], v_opt - 2.0 * config.v_step)
        vhi = min(config.v_bounds[1], v_opt + 2.0 * config.v_step)
        vg, vc = scan_slice(vlo, vhi, 15, chi2_at_v)
        out["v0_kms"] = halfwidth_chi2_plus_one(vg, vc, chi2, config.v_step)
        slo = max(config.sigma_bounds[0], s_opt - 2.0 * config.sigma_step)
        shi = min(config.sigma_bounds[1], s_opt + 2.0 * config.sigma_step)
        sg, sc = scan_slice(slo, shi, 15, chi2_at_s)
        out["sigma_kms"] = halfwidth_chi2_plus_one(sg, sc, chi2, config.sigma_step)
    return out


def fit_spectrum(
    wavelengths: ArrayLike,
    flux: ArrayLike,
    error: Optional[ArrayLike],
    base_matrix: ArrayLike,
    mask: Optional[ArrayLike] = None,
    config: Optional[FitConfig] = None,
    template_ages: Optional[ArrayLike] = None,
) -> FitResult:
    """Fit x_j and A_V; optionally grid-search v and sigma; clip/EX0 if configured."""
    if config is None:
        config = FitConfig()
    config = replace(config)

    wave = np.asarray(wavelengths, dtype=float)
    obs = np.asarray(flux, dtype=float)
    bases = np.asarray(base_matrix, dtype=float)
    n_wave, n_comp = bases.shape
    if obs.shape != (n_wave,):
        raise ValueError("Observation and base wavelength axes must match.")
    err_in: Optional[NDArray[np.float64]]
    if error is None:
        err_in = None
    else:
        err_in = np.asarray(error, dtype=float)
        if err_in.shape != (n_wave,):
            raise ValueError("Observation, error and base wavelength axes must match.")
    if np.any(np.diff(wave) <= 0):
        raise ValueError("wavelengths must be strictly increasing.")
    if not np.all(np.isfinite(wave)):
        raise ValueError("wavelengths must be finite.")
    if not np.all(np.isfinite(obs)):
        raise ValueError("flux contains non-finite values.")
    if not np.all(np.isfinite(bases)):
        raise ValueError("base_matrix contains non-finite values.")

    x_reg = XRegularizer.from_config(
        config.regularize_x,
        template_ages,
        n_comp,
        config.regularize_strength,
        config.age_bin_edges,
    )
    if config.error_method is not None:
        method = config.error_method.lower()
        if method not in ("chi2_slice", "repeat"):
            raise ValueError("error_method must be None, 'chi2_slice', or 'repeat'.")
        if method == "repeat" and config.n_repeat < 2:
            raise ValueError("n_repeat must be >= 2.")
    else:
        method = None

    if config.redshift < 0:
        raise ValueError("redshift must be >= 0.")
    if config.wave_frame.lower() not in ("as_is", "air", "vacuum"):
        raise ValueError("wave_frame must be 'as_is', 'air', or 'vacuum'.")
    if config.redshift != 0.0 or config.wave_frame.lower() == "air":
        obs, err_in = align_observation(
            wave, obs, err_in, redshift=config.redshift, wave_frame=config.wave_frame
        )

    if (config.fwhm_data is None) != (config.fwhm_template is None):
        raise ValueError("fwhm_data and fwhm_template must both be set or both omitted.")
    if config.fwhm_data is not None:
        if config.fwhm_data < 0.0 or config.fwhm_template < 0.0:
            raise ValueError("fwhm_data and fwhm_template must be >= 0.")
        if config.fwhm_data > config.fwhm_template:
            bases = match_instrumental_fwhm(wave, bases, config.fwhm_data, config.fwhm_template)
        elif config.fwhm_data < config.fwhm_template:
            warnings.warn(
                "fwhm_data < fwhm_template; cannot deconvolve. Skipping LSF match.",
                UserWarning,
                stacklevel=2,
            )

    good = np.ones(n_wave, dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
    if good.shape != (n_wave,):
        raise ValueError("mask must match the wavelength axis.")

    if err_in is None:
        if not config.estimate_error:
            raise ValueError("error is required unless estimate_error is True.")
        warnings.warn(
            "Estimated a constant RMS error; this is not a true χ².",
            UserWarning,
            stacklevel=2,
        )
        err = estimate_rms_error(obs, good)
    else:
        err = err_in
    obs_n, obs_scale = normalize_at(wave, obs, config.norm_window)
    bases_n = normalize_bases(wave, bases, config.norm_window)
    err_n = err / obs_scale

    if not config.search_kinematics and config.sigma_kms == 0.0:
        warnings.warn(
            "search_kinematics is False and sigma_kms is 0; velocity dispersion is not fitted.",
            UserWarning,
            stacklevel=2,
        )

    q = get_extinction_curve(wave, law=config.law, r_v=config.r_v)
    win = (wave >= config.norm_window[0]) & (wave <= config.norm_window[1])
    q0 = float(np.median(q[win]))
    weights = _weights(err_n, good)
    a_v_grid = _closed_grid(*config.a_v_bounds, config.a_v_step)

    if config.search_kinematics:
        v_grid = _closed_grid(*config.v_bounds, config.v_step)
        s_grid = _closed_grid(*config.sigma_bounds, config.sigma_step)
        best = None
        for v0 in v_grid:
            for sig in s_grid:
                x_hat, a_v, model_n, chi2, chi2_av = _best_av_for_kinematics(
                    wave, obs_n, weights, bases_n, q, q0, a_v_grid, float(v0), float(sig), x_reg=x_reg
                )
                if best is None or chi2 < best[0]:
                    best = (chi2, x_hat, a_v, model_n, chi2_av, float(v0), float(sig))
        if best is None:
            raise ValueError("Kinematic grids are empty.")
        chi2, x_opt, a_v_opt, model_n, chi2_av, v_opt, s_opt = best
    else:
        x_opt, a_v_opt, model_n, chi2, chi2_av = _best_av_for_kinematics(
            wave,
            obs_n,
            weights,
            bases_n,
            q,
            q0,
            a_v_grid,
            config.v0_kms,
            config.sigma_kms,
            x_reg=x_reg,
        )
        v_opt = config.v0_kms
        s_opt = config.sigma_kms

    if config.refine_kinematics:
        def _eval(av: float, vv: float, ss: float):
            return _design_and_nnls(wave, obs_n, weights, bases_n, q, q0, av, vv, ss, x_reg=x_reg)

        x_opt, a_v_opt, model_n, chi2, v_opt, s_opt = refine_av_v_sigma(
            _eval,
            a_v_opt,
            v_opt,
            s_opt,
            a_v_step=config.a_v_step,
            v_step=config.v_step,
            sigma_step=config.sigma_step,
            a_v_bounds=config.a_v_bounds,
            v_bounds=config.v_bounds,
            sigma_bounds=config.sigma_bounds,
            vary_kinematics=config.search_kinematics,
        )
        _, _, _, _, chi2_av = _best_av_for_kinematics(
            wave, obs_n, weights, bases_n, q, q0, a_v_grid, v_opt, s_opt, x_reg=x_reg
        )

    n_clipped = 0
    if config.clip_nsigma is not None:
        clipped_good = clip_outliers(obs_n, model_n, err_n, config.clip_nsigma, good)
        n_clipped = int(np.sum(good & ~clipped_good))
        if n_clipped > 0:
            good = clipped_good
            weights = _weights(err_n, good)
            x_opt, a_v_opt, model_n, chi2, chi2_av = _best_av_for_kinematics(
                wave,
                obs_n,
                weights,
                bases_n,
                q,
                q0,
                a_v_grid,
                v_opt,
                s_opt,
                x_reg=x_reg,
            )

    n_used = n_comp
    dropped = np.zeros(n_comp, dtype=bool)
    if config.x_min_keep > 0:
        x_sum = float(np.sum(x_opt))
        x_frac_now = x_opt / x_sum if x_sum > 0 else x_opt
        keep = keep_components(x_frac_now, config.x_min_keep)
        if np.any(~keep):
            x_sub, a_v_opt, model_n, chi2, chi2_av = _best_av_for_kinematics(
                wave,
                obs_n,
                weights,
                bases_n[:, keep],
                q,
                q0,
                a_v_grid,
                v_opt,
                s_opt,
                x_reg=x_reg.restrict(keep) if x_reg is not None else None,
            )
            x_full = np.zeros(n_comp, dtype=float)
            x_full[keep] = x_sub
            x_opt = x_full
            dropped = ~keep
            n_used = int(np.sum(keep))

    n_good = int(np.sum(weights > 0))
    dof = max(n_good - n_used - 1, 1)
    x_sum = float(np.sum(x_opt))
    x_frac = x_opt / x_sum if x_sum > 0 else x_opt

    extra_errors = None
    if method == "chi2_slice":
        extra_errors = _errors_chi2_slice(
            config,
            wave,
            obs_n,
            weights,
            bases_n,
            q,
            q0,
            a_v_opt,
            v_opt,
            s_opt,
            chi2,
            x_reg,
        )
    elif method == "repeat":
        inner = replace(
            config,
            error_method=None,
            redshift=0.0,
            wave_frame="as_is",
            fwhm_data=None,
            fwhm_template=None,
            clip_nsigma=None,
        )
        rng = np.random.default_rng(config.repeat_seed)
        model_obs = model_n * obs_scale
        samples = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            for _ in range(config.n_repeat):
                noisy = model_obs + rng.normal(0.0, err)
                rep = fit_spectrum(
                    wave,
                    noisy,
                    err,
                    bases,
                    mask=good,
                    config=inner,
                    template_ages=template_ages,
                )
                samples.append(rep.x_fraction)
        extra_errors = {"x": x_fraction_std(np.vstack(samples))}

    return FitResult(
        x=x_opt,
        x_fraction=x_frac,
        a_v=a_v_opt,
        v0_kms=v_opt,
        sigma_kms=s_opt,
        model=model_n * obs_scale,
        chi2=chi2,
        chi2_reduced=chi2 / dof,
        n_good=n_good,
        a_v_grid=a_v_grid,
        chi2_grid=chi2_av,
        n_clipped=n_clipped,
        good=good.copy(),
        obs_scale=float(obs_scale),
        config=config,
        dropped=dropped,
        errors=extra_errors,
    )
