"""q(λ) = A(λ)/A_V. Wavelengths in Angstroms."""

from __future__ import annotations

from typing import Optional

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.interpolate import interp1d

_GORDON_TABLES = {
    "GD1": [
        (0.125, 7.75), (0.150, 5.91), (0.175, 4.89), (0.200, 4.16),
        (0.250, 3.20), (0.300, 2.62), (0.440, 1.76), (0.550, 1.00),
        (0.700, 0.71), (0.900, 0.43),
    ],
    "GD2": [
        (0.125, 6.95), (0.150, 5.58), (0.175, 4.67), (0.200, 4.01),
        (0.250, 3.12), (0.300, 2.55), (0.440, 1.75), (0.550, 1.00),
        (0.700, 0.70), (0.900, 0.46),
    ],
    "GD3": [
        (0.125, 5.61), (0.150, 4.54), (0.175, 3.80), (0.200, 3.29),
        (0.250, 2.53), (0.300, 2.05), (0.440, 1.58), (0.550, 1.00),
        (0.700, 0.72), (0.900, 0.51),
    ],
}


def extinction_ccm89(wavelengths: ArrayLike, r_v: float = 3.1) -> NDArray[np.float64]:
    x = 1.0e4 / np.asarray(wavelengths, dtype=float)
    a = np.zeros_like(x)
    b = np.zeros_like(x)

    ir = (x >= 0.3) & (x < 1.1)
    a[ir] = 0.574 * x[ir] ** 1.61
    b[ir] = -0.527 * x[ir] ** 1.61

    opt = (x >= 1.1) & (x < 3.3)
    y = x[opt] - 1.82
    a[opt] = (
        1.0 + 0.17699 * y - 0.50447 * y**2 - 0.02427 * y**3
        + 0.72085 * y**4 + 0.01979 * y**5 - 0.77530 * y**6 + 0.32999 * y**7
    )
    b[opt] = (
        1.41338 * y + 2.28305 * y**2 + 1.07233 * y**3 - 5.38434 * y**4
        - 0.62251 * y**5 + 5.30260 * y**6 - 2.09002 * y**7
    )

    uv = (x >= 3.3) & (x <= 8.0)
    a[uv] = 1.752 - 0.316 * x[uv] - 0.104 / ((x[uv] - 4.67) ** 2 + 0.341)
    b[uv] = -3.090 + 1.825 * x[uv] + 1.206 / ((x[uv] - 4.62) ** 2 + 0.263)
    bump = uv & (x >= 5.9)
    y_uv = x[bump] - 5.9
    a[bump] += -0.04473 * y_uv**2 - 0.009779 * y_uv**3
    b[bump] += 0.2130 * y_uv**2 + 0.1207 * y_uv**3

    fuv = (x > 8.0) & (x <= 10.0)
    y_fuv = x[fuv] - 8.0
    a[fuv] = -1.073 - 0.628 * y_fuv + 0.137 * y_fuv**2 - 0.070 * y_fuv**3
    b[fuv] = 13.670 + 4.257 * y_fuv - 0.420 * y_fuv**2 + 0.374 * y_fuv**3
    return a + b / r_v


def extinction_calzetti(wavelengths: ArrayLike, r_v: float = 4.05) -> NDArray[np.float64]:
    wl = np.asarray(wavelengths, dtype=float) / 1e4
    k = np.zeros_like(wl)
    blue = (wl >= 0.12) & (wl <= 0.63)
    k[blue] = 2.659 * (-2.156 + 1.509 / wl[blue] - 0.198 / wl[blue] ** 2 + 0.011 / wl[blue] ** 3) + r_v
    red = (wl > 0.63) & (wl <= 2.2)
    k[red] = 2.659 * (-1.857 + 1.040 / wl[red]) + r_v
    return k / r_v


def extinction_gordon(wavelengths: ArrayLike, table: str = "GD1") -> NDArray[np.float64]:
    data = _GORDON_TABLES.get(table.upper())
    if data is None:
        raise ValueError(f"Unknown Gordon table {table!r}.")
    lam_um, q = zip(*data)
    interpolator = interp1d(lam_um, q, kind="linear", bounds_error=False, fill_value="extrapolate")
    return interpolator(np.asarray(wavelengths, dtype=float) * 1e-4)


def _dust_extinction_curve(
    wavelengths: ArrayLike, model_name: str, r_v: float
) -> NDArray[np.float64]:
    name = model_name.strip()
    if not name or not name.isidentifier() or name.startswith("_") or name.startswith("Base"):
        raise ValueError(f"Unknown dust_extinction model {model_name!r}.")
    try:
        from astropy import units as u
        from dust_extinction import parameter_averages
    except ImportError as exc:
        raise ImportError(
            "law='dust:...' requires the dust_extinction package. "
            "Install with: pip install 'starlightpy[dust]'."
        ) from exc
    cls = getattr(parameter_averages, name, None)
    if cls is None:
        matches = [
            item
            for item in dir(parameter_averages)
            if item.lower() == name.lower() and item.isidentifier() and not item.startswith("_")
        ]
        cls = getattr(parameter_averages, matches[0], None) if len(matches) == 1 else None
    if cls is None or not isinstance(cls, type) or not hasattr(cls, "evaluate"):
        raise ValueError(
            f"Unknown dust_extinction model {model_name!r}. "
            "Use a class name from dust_extinction.parameter_averages (e.g. F99)."
        )
    wave = np.asarray(wavelengths, dtype=float)
    try:
        model = cls(Rv=float(r_v))
    except TypeError:
        model = cls()
    return np.asarray(model(wave * u.AA), dtype=float)


def get_extinction_curve(
    wavelengths: ArrayLike,
    law: str = "CCM",
    r_v: float = 3.1,
    custom_law_path: Optional[str] = None,
) -> NDArray[np.float64]:
    raw = str(law).strip()
    if raw.lower().startswith("dust:"):
        return _dust_extinction_curve(wavelengths, raw.split(":", 1)[1], r_v)
    law = raw.upper()
    if law == "CCM":
        return extinction_ccm89(wavelengths, r_v)
    if law == "CAL":
        return extinction_calzetti(wavelengths, r_v if r_v != 3.1 else 4.05)
    if law in _GORDON_TABLES:
        return extinction_gordon(wavelengths, law)
    if custom_law_path is not None:
        data = np.loadtxt(custom_law_path)
        interpolator = interp1d(data[:, 0], data[:, 1], kind="linear", bounds_error=False, fill_value="extrapolate")
        return interpolator(np.asarray(wavelengths, dtype=float))
    raise ValueError(f"Reddening law {law!r} is not implemented.")
