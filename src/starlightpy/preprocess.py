"""Observation-frame helpers applied before fitting. Do not search redshift."""

from __future__ import annotations

from typing import List, Optional, Tuple

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.ndimage import gaussian_filter1d

from .io import resample_to

# FWHM = 2 sqrt(2 ln 2) σ
_FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))

# STARLIGHT-style (λ_ini, λ_fin, weight); weight 0 means fully masked.
# Optical nebular lines only; not a complete rest-UV / NIR catalogue.
_OPTICAL_EMISSION_MASK_REGIONS: Tuple[Tuple[float, float, float], ...] = (
    (3722.0, 3736.0, 0.0),  # [O II]
    (3865.0, 3885.0, 0.0),  # [Ne III] / He I
    (3964.0, 3978.0, 0.0),  # [Ne III]
    (4098.0, 4110.0, 0.0),  # Hδ
    (4336.0, 4348.0, 0.0),  # Hγ
    (4360.0, 4366.0, 0.0),  # [O III] 4363
    (4684.0, 4690.0, 0.0),  # He II
    (4856.0, 4868.0, 0.0),  # Hβ
    (4954.0, 4964.0, 0.0),  # [O III] 4959
    (5002.0, 5014.0, 0.0),  # [O III] 5007
    (5870.0, 5882.0, 0.0),  # He I
    (6296.0, 6308.0, 0.0),  # [O I]
    (6543.0, 6553.0, 0.0),  # [N II] 6548
    (6558.0, 6572.0, 0.0),  # Hα 6562.8
    (6578.0, 6588.0, 0.0),  # [N II] 6584
    (6711.0, 6724.0, 0.0),  # [S II] 6716
    (6728.0, 6742.0, 0.0),  # [S II] 6731
)


def to_rest_frame(
    wavelength: ArrayLike,
    flux: ArrayLike,
    redshift: float,
    error: Optional[ArrayLike] = None,
):
    """Move F_λ from observed to rest frame. ``redshift`` is applied, never fitted."""
    if redshift < 0:
        raise ValueError("redshift must be >= 0.")
    factor = 1.0 + float(redshift)
    wave = np.asarray(wavelength, dtype=float) / factor
    y = np.asarray(flux, dtype=float) * factor
    if error is None:
        return wave, y
    err = np.asarray(error, dtype=float) * factor
    return wave, y, err


def vacuum_to_air(wavelength: ArrayLike) -> NDArray[np.float64]:
    """Morton / SDSS-style conversion; wavelength in Å."""
    vac = np.asarray(wavelength, dtype=float)
    if np.any(vac <= 0):
        raise ValueError("Wavelengths must be positive.")
    return vac / (1.0 + 2.735182e-4 + 131.4182 / vac**2 + 2.76249e8 / vac**4)


def air_to_vacuum(wavelength: ArrayLike) -> NDArray[np.float64]:
    """Approximate inverse of ``vacuum_to_air``; wavelength in Å."""
    air = np.asarray(wavelength, dtype=float)
    if np.any(air <= 0):
        raise ValueError("Wavelengths must be positive.")
    return air * (1.0 + 2.735182e-4 + 131.4182 / air**2 + 2.76249e8 / air**4)


def align_observation(
    template_wave: ArrayLike,
    flux: ArrayLike,
    error: Optional[ArrayLike] = None,
    redshift: float = 0.0,
    wave_frame: str = "as_is",
) -> Tuple[NDArray[np.float64], Optional[NDArray[np.float64]]]:
    """Put observed F_λ onto the rest-frame vacuum template grid."""
    frame = wave_frame.lower()
    if frame not in ("as_is", "air", "vacuum"):
        raise ValueError("wave_frame must be 'as_is', 'air', or 'vacuum'.")
    target = np.asarray(template_wave, dtype=float)
    wave_data = air_to_vacuum(target) if frame == "air" else target.copy()
    if redshift == 0.0 and frame != "air":
        flux_out = np.asarray(flux, dtype=float)
        if error is None:
            return flux_out, None
        return flux_out, np.asarray(error, dtype=float)
    if error is None:
        wave_rest, flux_rest = to_rest_frame(wave_data, flux, redshift)
        return resample_to(wave_rest, flux_rest, target), None
    wave_rest, flux_rest, err_rest = to_rest_frame(wave_data, flux, redshift, error)
    return resample_to(wave_rest, flux_rest, target), resample_to(wave_rest, err_rest, target)


def extra_instrumental_fwhm(fwhm_data: float, fwhm_template: float) -> float:
    """Quadratic extra FWHM in Å. Zero if data is not broader (no deconvolution)."""
    if fwhm_data < 0 or fwhm_template < 0:
        raise ValueError("FWHM must be >= 0.")
    if fwhm_data <= fwhm_template:
        return 0.0
    return float(np.sqrt(fwhm_data**2 - fwhm_template**2))


def _broaden_linear_angstrom(
    wave: NDArray[np.float64],
    flux_1d: NDArray[np.float64],
    sigma_angstrom: float,
) -> NDArray[np.float64]:
    """Gaussian in linear Å (instrumental LSF), not velocity-space LOSVD."""
    if sigma_angstrom < 1e-9:
        return flux_1d.copy()
    if wave.size != flux_1d.size:
        raise ValueError("wavelength and flux must have the same length.")
    if wave.size < 2:
        raise ValueError("Need at least two wavelength samples to broaden.")
    dw = np.diff(wave)
    if np.any(dw <= 0):
        raise ValueError("wavelengths must be strictly increasing.")
    if np.allclose(dw, dw[0], rtol=1e-6, atol=1e-8):
        sigma_pix = sigma_angstrom / float(dw[0])
        if sigma_pix < 0.05:
            return flux_1d.copy()
        return gaussian_filter1d(flux_1d, sigma=sigma_pix, mode="nearest")
    dw_u = float(np.min(dw))
    n = int(np.ceil((wave[-1] - wave[0]) / dw_u)) + 1
    wave_u = np.linspace(wave[0], wave[-1], n)
    y_u = np.interp(wave_u, wave, flux_1d)
    sigma_pix = sigma_angstrom / float(wave_u[1] - wave_u[0])
    if sigma_pix >= 0.05:
        y_u = gaussian_filter1d(y_u, sigma=sigma_pix, mode="nearest")
    return np.interp(wave, wave_u, y_u)


def match_instrumental_fwhm(
    wavelength: ArrayLike,
    flux: ArrayLike,
    fwhm_data: float,
    fwhm_template: float,
) -> NDArray[np.float64]:
    """Broaden templates (or a 1-D spectrum) from ``fwhm_template`` to ``fwhm_data``.

    Units are Å FWHM. If ``fwhm_data <= fwhm_template`` the array is copied:
    this helper never sharpens.
    """
    extra = extra_instrumental_fwhm(fwhm_data, fwhm_template)
    sigma = extra * _FWHM_TO_SIGMA
    wave = np.asarray(wavelength, dtype=float)
    y = np.asarray(flux, dtype=float)
    if y.ndim == 1:
        return _broaden_linear_angstrom(wave, y, sigma)
    if y.ndim == 2:
        if y.shape[0] != wave.size:
            raise ValueError("base wavelength axis must match wavelength.")
        out = np.empty_like(y, dtype=float)
        for j in range(y.shape[1]):
            out[:, j] = _broaden_linear_angstrom(wave, y[:, j], sigma)
        return out
    raise ValueError("flux must be 1-D or a (n_wave, n_comp) matrix.")


def optical_emission_mask_regions() -> List[Tuple[float, float, float]]:
    """STARLIGHT-style optical emission windows for ``apply_mask``."""
    return list(_OPTICAL_EMISSION_MASK_REGIONS)


def estimate_rms_error(
    flux: ArrayLike,
    good: Optional[ArrayLike] = None,
) -> NDArray[np.float64]:
    """Constant RMS of unmasked flux. This is not a true per-pixel error spectrum."""
    y = np.asarray(flux, dtype=float)
    used = np.isfinite(y)
    if good is not None:
        used &= np.asarray(good, dtype=bool)
    if int(np.sum(used)) < 2:
        raise ValueError("Need at least two unmasked finite pixels to estimate RMS error.")
    rms = float(np.std(y[used], ddof=1))
    if not np.isfinite(rms) or rms <= 0.0:
        med = float(np.median(np.abs(y[used])))
        rms = 0.01 * med if med > 0.0 else 0.0
    if rms <= 0.0:
        raise ValueError("Cannot estimate a positive RMS error from the flux.")
    return np.full(y.shape, rms, dtype=float)
