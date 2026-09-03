"""Call pPXF after the boilerplate astronomers usually copy from the examples.

pPXF wants log-rebinned spectra and templates on the same velocity scale.
This module takes linearly sampled wavelength + flux and does that step,
then returns kinematics. Cite Cappellari (see FitResult.cite); do not cite this
wrapper as a fitting method.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from ppxf.ppxf import ppxf
from ppxf.ppxf_util import determine_goodpixels, log_rebin

CITE = (
    "Cappellari (2023, MNRAS, 526, 3273); "
    "Cappellari (2017, MNRAS, 466, 798); "
    "Cappellari & Emsellem (2004, PASP, 116, 138)"
)

# pPXF ``set_lam_input`` default velocity bounds when ``bounds`` is omitted.
_C_KMS = 299792.458
_PPXF_DEFAULT_V_MARGIN_KMS = 2900.0


def _goodpixels_within_ppxf_template_margin(
    lam_gal: NDArray[np.float64],
    lam_temp: NDArray[np.float64],
    goodpixels: Optional[NDArray[np.int_]],
    v_margin_kms: float = _PPXF_DEFAULT_V_MARGIN_KMS,
) -> NDArray[np.int_]:
    """Keep galaxy pixels pPXF can cover with the default ±2900 km/s bounds.

    Passing ``lam`` and ``lam_temp`` makes pPXF require extra template
    wavelength coverage. Same-grid calls used to hit an ``AssertionError``.
    """
    factor = float(np.exp(v_margin_kms / _C_KMS))
    lo = float(np.min(lam_temp)) * factor
    hi = float(np.max(lam_temp)) / factor
    if (not np.isfinite(lo)) or (not np.isfinite(hi)) or lo >= hi:
        raise ValueError(
            "templates wavelength range is too short for pPXF "
            f"(need about ±{int(v_margin_kms)} km/s extra coverage beyond the galaxy)."
        )
    cover = np.flatnonzero((lam_gal > lo) & (lam_gal < hi))
    if goodpixels is None:
        selected = cover
    else:
        selected = np.intersect1d(np.asarray(goodpixels, dtype=int), cover)
    selected = np.asarray(np.sort(selected), dtype=int)
    if selected.size < 10:
        raise ValueError(
            "Too few galaxy pixels remain after leaving the pPXF velocity margin "
            f"(about ±{int(v_margin_kms)} km/s). "
            "Give templates a wider wavelength range than the galaxy."
        )
    return selected


@dataclass
class FitResult:
    """Kinematics from pPXF. ``bestfit`` and ``wavelength`` are on the log grid."""

    velocity: float
    sigma: float
    sol: NDArray[np.float64]
    error: NDArray[np.float64]
    chi2: float
    bestfit: NDArray[np.float64]
    wavelength: NDArray[np.float64]
    galaxy: NDArray[np.float64]
    weights: NDArray[np.float64]
    velscale: float
    cite: str = CITE


def _as_2d_templates(spec: NDArray[np.float64]) -> NDArray[np.float64]:
    if spec.ndim == 1:
        return spec[:, np.newaxis]
    if spec.ndim == 2:
        return spec
    raise ValueError("templates must be (n_wave,) or (n_wave, n_templates)")


def fit_spectrum(
    wavelength: ArrayLike,
    flux: ArrayLike,
    templates: ArrayLike,
    template_wavelength: ArrayLike,
    error: Optional[ArrayLike] = None,
    *,
    redshift: float = 0.0,
    start: Sequence[float] = (0.0, 150.0),
    moments: int = 2,
    degree: int = 4,
    mask_emission: bool = True,
    quiet: bool = True,
) -> FitResult:
    """Fit one rest-frame-ish spectrum with user-supplied templates via pPXF.

    Parameters
    ----------
    wavelength, flux
        Observed spectrum, linear or irregular sampling, Angstroms.
    templates
        Template fluxes ``(n_wave,)`` or ``(n_wave, n_templates)`` on
        ``template_wavelength``.
    error
        Optional 1σ error, same sampling as ``flux``. If omitted, a constant
        1% of the median flux is used (only for a first look).
    redshift
        Approximate redshift; the wavelength axis is divided by ``1+z``
        before the fit (pPXF example convention).
    """
    wave = np.asarray(wavelength, dtype=float)
    flux = np.asarray(flux, dtype=float)
    t_wave = np.asarray(template_wavelength, dtype=float)
    t_lin = _as_2d_templates(np.asarray(templates, dtype=float))

    if wave.shape != flux.shape:
        raise ValueError("wavelength and flux must have the same shape.")
    if t_lin.shape[0] != t_wave.size:
        raise ValueError("templates first axis must match template_wavelength.")

    if redshift != 0.0:
        wave = wave / (1.0 + redshift)

    galaxy_log, ln_lam, velscale = log_rebin(wave, flux)
    velscale = float(np.squeeze(velscale))
    scale = float(np.median(galaxy_log))
    if scale == 0:
        raise ValueError("Spectrum median is zero after log-rebin.")
    galaxy = galaxy_log / scale

    if error is None:
        noise = np.full_like(galaxy, 0.01)
    else:
        err = np.asarray(error, dtype=float)
        noise_log, _, _ = log_rebin(np.asarray(wavelength, dtype=float) / (1.0 + redshift), err, velscale=velscale)
        noise = np.clip(noise_log / scale, 1e-8, None)

    templates_log, ln_temp, _ = log_rebin(t_wave, t_lin, velscale=velscale)
    templates_log = _as_2d_templates(np.asarray(templates_log, dtype=float))
    templates_log = templates_log / np.median(templates_log, axis=0, keepdims=True)

    lam_gal = np.exp(ln_lam)
    lam_temp = np.exp(ln_temp)

    goodpixels = None
    if mask_emission:
        goodpixels = determine_goodpixels(ln_lam, [lam_temp.min(), lam_temp.max()])
    goodpixels = _goodpixels_within_ppxf_template_margin(lam_gal, lam_temp, goodpixels)

    # When both wavelength vectors are passed, pPXF forbids vsyst.
    pp = ppxf(
        templates_log,
        galaxy,
        noise,
        velscale,
        list(start),
        goodpixels=goodpixels,
        moments=moments,
        degree=degree,
        lam=lam_gal,
        lam_temp=lam_temp,
        quiet=quiet,
    )
    sol = np.atleast_1d(np.asarray(pp.sol, dtype=float))
    err = np.atleast_1d(np.asarray(pp.error, dtype=float)) * np.sqrt(max(pp.chi2, 1e-12))
    return FitResult(
        velocity=float(sol[0]),
        sigma=float(sol[1]) if sol.size > 1 else float("nan"),
        sol=sol,
        error=err,
        chi2=float(pp.chi2),
        bestfit=np.asarray(pp.bestfit, dtype=float),
        wavelength=lam_gal,
        galaxy=galaxy,
        weights=np.asarray(pp.weights, dtype=float),
        velscale=velscale,
    )
