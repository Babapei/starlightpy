"""STARLIGHT-style ASCII readers. Interpolation onto a common grid is the caller's job."""

from __future__ import annotations

from pathlib import Path
from typing import List, Sequence, Tuple, Union

import numpy as np
from numpy.typing import NDArray

PathLike = Union[str, Path]


def load_spectrum(
    filename: PathLike,
    skip_header: bool = False,
) -> Tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.int_]]:
    wavelengths: List[float] = []
    fluxes: List[float] = []
    errors: List[float] = []
    flags: List[int] = []

    with open(filename, "r", encoding="utf-8") as handle:
        lines = handle.readlines()
    if skip_header and lines:
        lines = lines[1:]

    for raw in lines:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
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


def load_base_master(master_path: PathLike) -> List[dict]:
    with open(master_path, "r", encoding="utf-8") as handle:
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
        path = base_dir / base["filename"]
        wavelengths: List[float] = []
        fluxes: List[float] = []
        with open(path, "r", encoding="utf-8") as handle:
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
    with open(filename, "r", encoding="utf-8") as handle:
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
