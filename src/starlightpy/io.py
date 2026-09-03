"""STARLIGHT-style ASCII readers and FitResult npz/json save (not Fortran .out)."""

from __future__ import annotations

import gzip
import json
from dataclasses import asdict, fields
from pathlib import Path
from typing import Any, Iterator, List, Optional, Sequence, Tuple, Union

import numpy as np
from numpy.typing import ArrayLike, NDArray

PathLike = Union[str, Path]


def _open_text(path: Path):
    name = path.name.lower()
    if name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def _resolve_existing(path: Path) -> Path:
    if path.exists():
        return path
    gz = Path(str(path) + ".gz")
    if gz.exists():
        return gz
    raise FileNotFoundError(path)


def load_spectrum(
    filename: PathLike,
    skip_header: bool = False,
) -> Tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.int_]]:
    """Read λ, flux, error, flag from STARLIGHT-style ASCII (optionally gzipped).

    After comments (``#``) and an optional skipped header line, a single-token
    line is treated as ``Npix`` and ignored. Remaining rows are 2–4 columns:
    wavelength, flux, optional error, optional integer flag.
    """
    path = _resolve_existing(Path(filename))
    with _open_text(path) as handle:
        lines = handle.readlines()
    if skip_header and lines:
        lines = lines[1:]

    data_lines: List[str] = []
    for raw in lines:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        data_lines.append(line)
    if data_lines and len(data_lines[0].split()) == 1:
        data_lines = data_lines[1:]

    wavelengths: List[float] = []
    fluxes: List[float] = []
    errors: List[float] = []
    flags: List[int] = []
    for line in data_lines:
        parts = line.split()
        if len(parts) == 4:
            wl, fl, err, flag = parts
            wavelengths.append(float(wl))
            fluxes.append(float(fl))
            errors.append(float(err))
            flags.append(int(float(flag)))
        elif len(parts) == 3:
            wl, fl, err = parts
            wavelengths.append(float(wl))
            fluxes.append(float(fl))
            errors.append(float(err))
            flags.append(0)
        elif len(parts) == 2:
            wl, fl = parts
            wavelengths.append(float(wl))
            fluxes.append(float(fl))
            errors.append(1.0)
            flags.append(0)
        else:
            raise ValueError(f"Invalid spectrum line: {line}")

    return (
        np.asarray(wavelengths),
        np.asarray(fluxes),
        np.asarray(errors),
        np.asarray(flags, dtype=int),
    )


def load_cxt(
    filename: PathLike,
    skip_header: bool = False,
) -> Tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.int_]]:
    """Alias of ``load_spectrum`` for STARLIGHT ``.cxt`` names."""
    return load_spectrum(filename, skip_header=skip_header)


def good_from_flags(flags: ArrayLike) -> NDArray[np.bool_]:
    """STARLIGHT: ``flag >= 2`` is ignored in the fit."""
    return np.asarray(flags, dtype=int) < 2


def load_base_master(master_path: PathLike) -> List[dict]:
    path = _resolve_existing(Path(master_path))
    with _open_text(path) as handle:
        lines = [line for line in handle.readlines() if line.strip() and not line.startswith("#")]
    if not lines:
        raise ValueError("Base master file is empty.")
    n_base = int(lines[0].split()[0])
    base_list = []
    for line in lines[1 : n_base + 1]:
        parts = line.split()
        if len(parts) < 7:
            raise ValueError(f"Base line malformed: {line}")
        base_list.append(
            {
                "filename": parts[0],
                "age": float(parts[1]),
                "Z": float(parts[2]),
                "nickname": parts[3],
                "f_star": float(parts[4]),
                "YAV": int(parts[5]),
                "alpha_Fe": float(parts[6]),
            }
        )
    return base_list


def load_base_spectra(base_dir: PathLike, base_list: Sequence[dict]) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    base_dir = Path(base_dir)
    spectra = []
    wl_ref = None
    for base in base_list:
        path = _resolve_existing(base_dir / base["filename"])
        wavelengths: List[float] = []
        fluxes: List[float] = []
        with _open_text(path) as handle:
            for raw in handle:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 2:
                    continue
                wavelengths.append(float(parts[0]))
                fluxes.append(float(parts[1]))
        wave = np.asarray(wavelengths)
        flux = np.asarray(fluxes)
        if wl_ref is None:
            wl_ref = wave
        elif not np.allclose(wl_ref, wave):
            raise ValueError(f"Base file {path} has a different wavelength sampling.")
        spectra.append(flux)
    if wl_ref is None:
        raise ValueError("No base spectra were loaded.")
    return wl_ref, np.column_stack(spectra)


def load_mask(filename: PathLike) -> List[Tuple[float, float, float]]:
    path = _resolve_existing(Path(filename))
    with _open_text(path) as handle:
        lines = [line for line in handle.readlines() if line.strip() and not line.startswith("#")]
    if not lines:
        raise ValueError("Mask file is empty.")
    n_masks = int(lines[0].split()[0])
    regions: List[Tuple[float, float, float]] = []
    for line in lines[1 : n_masks + 1]:
        parts = line.split()
        if len(parts) < 3:
            continue
        regions.append((float(parts[0]), float(parts[1]), float(parts[2])))
    return regions


def apply_mask(
    wavelengths: NDArray[np.float64],
    mask_regions: Sequence[Tuple[float, float, float]],
) -> NDArray[np.bool_]:
    good = np.ones_like(wavelengths, dtype=bool)
    for lam_ini, lam_fin, _weight in mask_regions:
        good &= ~((wavelengths >= lam_ini) & (wavelengths <= lam_fin))
    return good


def combine_good(
    wavelengths: ArrayLike,
    flags: Optional[ArrayLike] = None,
    mask_regions: Optional[Sequence[Tuple[float, float, float]]] = None,
) -> NDArray[np.bool_]:
    """AND of STARLIGHT flags and mask-file windows."""
    wave = np.asarray(wavelengths, dtype=float)
    good = np.ones(wave.shape, dtype=bool)
    if flags is not None:
        good &= good_from_flags(flags)
    if mask_regions is not None:
        good &= apply_mask(wave, mask_regions)
    return good


def resample_to(wavelength, flux, wave_out):
    """Linear interpolation onto ``wave_out``. Call this before ``fit_spectrum``."""
    wave = np.asarray(wavelength, dtype=float)
    y = np.asarray(flux, dtype=float)
    target = np.asarray(wave_out, dtype=float)
    if wave.size != y.size:
        raise ValueError("wavelength and flux must have the same length.")
    if wave.size < 2:
        raise ValueError("Need at least two wavelength samples to resample.")
    if np.any(np.diff(wave) <= 0) or np.any(np.diff(target) <= 0):
        raise ValueError("Wavelength arrays must be strictly increasing.")
    return np.interp(target, wave, y)


def iter_ascii_spectra(
    directory: PathLike,
    pattern: str = "*.cxt",
    skip_header: bool = False,
) -> Iterator[Tuple[Path, Tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.int_]]]]:
    """Yield ``(path, load_spectrum(...))`` for files in ``directory``. Not a job queue."""
    folder = Path(directory)
    paths = sorted(folder.glob(pattern))
    for path in paths:
        yield path, load_spectrum(path, skip_header=skip_header)


def load_sdss_fits(
    filename: PathLike,
) -> Tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.bool_]]:
    """Read a 1D SDSS-like FITS spectrum. Requires astropy (optional extra ``fits``).

    Accepts a table HDU with ``flux`` plus ``loglam`` or ``wavelength``/``wave``,
    optional ``ivar``/``error``, optional ``and_mask``/``mask``; or a primary 1D
    array with ``COEFF0``/``COEFF1`` or ``CRVAL1``.
    """
    try:
        from astropy.io import fits
    except ImportError as exc:
        raise ImportError("load_sdss_fits requires astropy. pip install 'starlightpy[fits]'") from exc

    path = Path(filename)
    with fits.open(path) as hdul:
        table_hdu = None
        for hdu in hdul:
            data = getattr(hdu, "data", None)
            columns = getattr(hdu, "columns", None)
            if data is None or columns is None:
                continue
            names = [c.name.lower() for c in columns]
            if "flux" in names:
                table_hdu = hdu
                break
        if table_hdu is not None:
            data = table_hdu.data
            names = {n.lower(): n for n in data.dtype.names}
            flux = np.asarray(data[names["flux"]], dtype=float).reshape(-1)
            if "loglam" in names:
                wave = 10.0 ** np.asarray(data[names["loglam"]], dtype=float).reshape(-1)
            else:
                wave = None
                for key in ("wavelength", "lambda", "wave", "lam"):
                    if key in names:
                        wave = np.asarray(data[names[key]], dtype=float).reshape(-1)
                        break
                if wave is None:
                    raise ValueError("FITS table has flux but no wavelength/loglam column.")
            if "ivar" in names:
                ivar = np.asarray(data[names["ivar"]], dtype=float).reshape(-1)
                err = np.full_like(flux, np.inf)
                positive = ivar > 0
                err[positive] = 1.0 / np.sqrt(ivar[positive])
            else:
                if "error" in names:
                    err = np.asarray(data[names["error"]], dtype=float).reshape(-1)
                else:
                    median = float(np.median(np.abs(flux[np.isfinite(flux)]))) if flux.size else 1.0
                    err = np.full_like(flux, max(0.01 * median, 1e-8))
            good = np.isfinite(flux) & np.isfinite(err) & (err < np.inf)
            if "ivar" in names:
                good &= np.asarray(data[names["ivar"]], dtype=float).reshape(-1) > 0
            for key in ("and_mask", "mask", "or_mask"):
                if key in names:
                    good &= np.asarray(data[names[key]]).reshape(-1) == 0
                    break
            return wave, flux, err, good

        primary = hdul[0]
        if primary.data is None:
            raise ValueError("FITS file has neither a flux table nor a primary array.")
        flux = np.asarray(primary.data, dtype=float)
        if flux.ndim > 1:
            flux = np.asarray(flux[0], dtype=float)
        hdr = primary.header
        if "COEFF0" in hdr and "COEFF1" in hdr:
            wave = 10.0 ** (float(hdr["COEFF0"]) + float(hdr["COEFF1"]) * np.arange(flux.size))
        elif "CRVAL1" in hdr:
            dx = float(hdr.get("CD1_1", hdr.get("CDELT1", 1.0)))
            wave = float(hdr["CRVAL1"]) + dx * np.arange(flux.size)
            ctype = str(hdr.get("CTYPE1", "")).upper()
            if ctype.startswith("LOG"):
                wave = 10.0 ** wave
        else:
            raise ValueError("No wavelength calibration (COEFF0 or CRVAL1) in the primary HDU.")
        median = float(np.median(np.abs(flux[np.isfinite(flux)]))) if flux.size else 1.0
        err = np.full_like(flux, max(0.01 * median, 1e-8))
        good = np.isfinite(flux)
        return wave, flux, err, good


_RESULT_FORMAT = "starlightpy-fitresult"
_RESULT_VERSION = 1
_CONFIG_TUPLE_FIELDS = {
    "norm_window",
    "a_v_bounds",
    "v_bounds",
    "sigma_bounds",
    "a_yv_bounds",
    "age_bin_edges",
}


def _jsonable(obj: Any) -> Any:
    if obj is None:
        return None
    if isinstance(obj, dict):
        return {str(key): _jsonable(value) for key, value in obj.items()}
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer, np.bool_)):
        return obj.item()
    if isinstance(obj, (list, tuple)):
        return [_jsonable(item) for item in obj]
    return obj


def _errors_from_json(data: Any) -> Optional[dict]:
    if data is None:
        return None
    out: dict = {}
    for key, value in data.items():
        if isinstance(value, list):
            out[key] = np.asarray(value, dtype=float)
        else:
            out[key] = value
    return out


def _config_to_json(config: Any) -> Optional[dict]:
    if config is None:
        return None
    return _jsonable(asdict(config))


def _config_from_json(data: Any):
    from .config import FitConfig

    if data is None:
        return None
    allowed = {item.name for item in fields(FitConfig)}
    kwargs = {}
    for key, value in data.items():
        if key not in allowed:
            continue
        if key in _CONFIG_TUPLE_FIELDS and value is not None:
            kwargs[key] = tuple(value)
        else:
            kwargs[key] = value
    return FitConfig(**kwargs)


def _result_to_payload(result: Any) -> dict:
    return {
        "format": _RESULT_FORMAT,
        "version": _RESULT_VERSION,
        "x": np.asarray(result.x, dtype=float).tolist(),
        "x_fraction": np.asarray(result.x_fraction, dtype=float).tolist(),
        "a_v": float(result.a_v),
        "a_yv": float(result.a_yv),
        "v0_kms": float(result.v0_kms),
        "sigma_kms": float(result.sigma_kms),
        "model": np.asarray(result.model, dtype=float).tolist(),
        "chi2": float(result.chi2),
        "chi2_reduced": float(result.chi2_reduced),
        "n_good": int(result.n_good),
        "a_v_grid": np.asarray(result.a_v_grid, dtype=float).tolist(),
        "chi2_grid": np.asarray(result.chi2_grid, dtype=float).tolist(),
        "n_clipped": int(result.n_clipped),
        "good": None if result.good is None else np.asarray(result.good, dtype=bool).tolist(),
        "obs_scale": float(result.obs_scale),
        "dropped": None
        if result.dropped is None
        else np.asarray(result.dropped, dtype=bool).tolist(),
        "errors": _jsonable(result.errors),
        "config": _config_to_json(result.config),
    }


def _payload_to_result(data: dict):
    from .fit import FitResult

    if data.get("format") != _RESULT_FORMAT:
        raise ValueError("Not a starlightpy FitResult file.")
    good = data.get("good")
    dropped = data.get("dropped")
    return FitResult(
        x=np.asarray(data["x"], dtype=float),
        x_fraction=np.asarray(data["x_fraction"], dtype=float),
        a_v=float(data["a_v"]),
        v0_kms=float(data["v0_kms"]),
        sigma_kms=float(data["sigma_kms"]),
        model=np.asarray(data["model"], dtype=float),
        chi2=float(data["chi2"]),
        chi2_reduced=float(data["chi2_reduced"]),
        n_good=int(data["n_good"]),
        a_v_grid=np.asarray(data["a_v_grid"], dtype=float),
        chi2_grid=np.asarray(data["chi2_grid"], dtype=float),
        n_clipped=int(data.get("n_clipped", 0)),
        good=None if good is None else np.asarray(good, dtype=bool),
        obs_scale=float(data.get("obs_scale", 1.0)),
        config=_config_from_json(data.get("config")),
        dropped=None if dropped is None else np.asarray(dropped, dtype=bool),
        errors=_errors_from_json(data.get("errors")),
        a_yv=float(data.get("a_yv", 0.0)),
    )


def save_fit_result(path: PathLike, result: Any) -> None:
    """Write a FitResult to ``.npz`` or ``.json``. Fortran ``.out`` is not supported."""
    dest = Path(path)
    suffix = dest.suffix.lower()
    payload = _result_to_payload(result)
    if suffix == ".json":
        dest.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        return
    if suffix == ".npz":
        np.savez_compressed(
            dest,
            payload=np.frombuffer(json.dumps(payload).encode("utf-8"), dtype=np.uint8),
        )
        return
    raise ValueError(
        "save_fit_result only writes .npz or .json; Fortran .out is not a product."
    )


def load_fit_result(path: PathLike):
    """Read a FitResult saved by ``save_fit_result``."""
    src = Path(path)
    suffix = src.suffix.lower()
    if suffix == ".json":
        data = json.loads(src.read_text(encoding="utf-8"))
        return _payload_to_result(data)
    if suffix == ".npz":
        with np.load(src, allow_pickle=False) as npz:
            if "payload" not in npz:
                raise ValueError("Not a starlightpy FitResult file.")
            data = json.loads(npz["payload"].tobytes().decode("utf-8"))
        return _payload_to_result(data)
    raise ValueError(
        "load_fit_result only reads .npz or .json; Fortran .out is not a product."
    )
