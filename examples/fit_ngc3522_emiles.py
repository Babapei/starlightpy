#!/usr/bin/env python3
"""Fit public NGC 3522 (SDSS) with a user-downloaded E-MILES subset.

This is a workflow demo, not a standard stellar-population solution.
Templates are not bundled in starlightpy; they are fetched like a pPXF user
would (Cappellari's ppxf_data, Vazdekis et al. 2016 E-MILES).
"""
from __future__ import annotations

import argparse
import warnings
from pathlib import Path
from urllib.request import urlretrieve

import numpy as np

from starlightpy import (
    FitConfig,
    air_to_vacuum,
    apply_mask,
    fit_spectrum,
    light_to_mass,
    light_weighted_age,
    load_sdss_fits,
    optical_emission_mask_regions,
    resample_to,
    save_fit_result,
    to_rest_frame,
)

SDSS_URL = (
    "https://data.sdss.org/sas/dr17/sdss/spectro/redux/26/spectra/lite/"
    "2488/spec-2488-54149-0001.fits"
)
EMILES_URL = "https://raw.githubusercontent.com/micappe/ppxf_data/main/spectra_emiles_9.0.npz"
AGE_IDX = (0, 6, 12, 17, 22, 24)
MET_IDX = (2, 4, 5)


def cache_dir() -> Path:
    return Path.home() / ".cache" / "starlightpy"


def fetch(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        print(f"downloading {url} -> {dest}")
        urlretrieve(url, dest)
    return dest


def emiles_subset(path: Path):
    npz = np.load(path)
    wave_vac = air_to_vacuum(np.asarray(npz["lam"], dtype=float))
    ages_gyr = np.asarray(npz["ages"], dtype=float)
    metals = np.asarray(npz["metals"], dtype=float)
    masses = np.asarray(npz["masses"], dtype=float)
    tpl = np.asarray(npz["templates"], dtype=float)
    cols, age_yr, z_list, mass_list = [], [], [], []
    for ia in AGE_IDX:
        for iz in MET_IDX:
            cols.append(tpl[:, ia, iz])
            age_yr.append(ages_gyr[ia] * 1e9)
            z_list.append(metals[iz])
            mass_list.append(masses[ia, iz])
    return (
        wave_vac,
        np.column_stack(cols),
        np.asarray(age_yr),
        np.asarray(z_list),
        np.asarray(mass_list),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=cache_dir())
    parser.add_argument("--out", type=Path, default=Path("ngc3522_starlightpy.npz"))
    parser.add_argument("--plot", type=Path, default=None, help="optional PNG path")
    args = parser.parse_args()

    spec = fetch(SDSS_URL, args.cache / "spec-2488-54149-0001.fits")
    emiles = fetch(EMILES_URL, args.cache / "spectra_emiles_9.0.npz")

    from astropy.io import fits

    wave_obs, flux_obs, err_obs, good_obs = load_sdss_fits(spec)
    with fits.open(spec) as hdul:
        z = float(hdul["SPECOBJ"].data["z"].item())

    wave_tpl, bases_full, ages, _metals, mass_star = emiles_subset(emiles)
    grid = np.arange(3800.0, 7200.1, 1.0)
    bases = np.column_stack(
        [resample_to(wave_tpl, bases_full[:, j], grid) for j in range(bases_full.shape[1])]
    )
    wave_rest, flux_rest, err_rest = to_rest_frame(wave_obs, flux_obs, z, err_obs)
    flux = resample_to(wave_rest, flux_rest, grid)
    error = resample_to(wave_rest, err_rest, grid)
    good = resample_to(wave_rest, good_obs.astype(float), grid) > 0.5
    good &= apply_mask(grid, optical_emission_mask_regions())
    good &= np.isfinite(flux) & np.isfinite(error) & (error > 0)

    win = (grid >= 4010.0) & (grid <= 4060.0)
    light = np.median(bases[win], axis=0)
    ml = mass_star / np.clip(light, 1e-30, None)
    ml = ml / np.median(ml)

    config = FitConfig(
        redshift=0.0,
        wave_frame="vacuum",
        fwhm_data=2.51,
        fwhm_template=2.51,
        search_kinematics=True,
        refine_kinematics=True,
        clip_nsigma=3.0,
        a_v_bounds=(0.0, 1.2),
        a_v_step=0.1,
        v_bounds=(-200.0, 200.0),
        v_step=50.0,
        sigma_bounds=(40.0, 220.0),
        sigma_step=40.0,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = fit_spectrum(grid, flux, error, bases, mask=good, config=config)

    age_L = light_weighted_age(result.x_fraction, ages)
    _, mu = light_to_mass(result.x_fraction, ml)
    print(
        f"NGC 3522  z={z:.5f}  A_V={result.a_v:.3f}  v={result.v0_kms:.1f}  "
        f"sigma={result.sigma_kms:.1f}  chi2/nu={result.chi2_reduced:.3f}  "
        f"age_L={age_L/1e9:.2f} Gyr"
    )
    print("x_fraction", np.round(result.x_fraction, 3))
    print("mu        ", np.round(mu, 3))
    print("Not a standard SFH solution; templates and grid are a demo subset.")
    save_fit_result(args.out, result)
    print("wrote", args.out)

    if args.plot is not None:
        from matplotlib import pyplot as plt

        fig, axes = plt.subplots(
            2, 1, figsize=(11, 6), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]}
        )
        ax, axr = axes
        ax.plot(result.wavelength, result.flux, color="0.45", lw=0.6, label="observed")
        ax.plot(result.wavelength, result.model, color="C3", lw=0.8, label="model")
        if result.good is not None:
            bad = ~result.good
            ax.plot(result.wavelength[bad], result.flux[bad], "o", ms=2.0, color="C0", alpha=0.65, label="masked")
        ax.set_ylabel("flux")
        ax.set_title(
            f"NGC 3522  demo fit (not a standard SFH)  "
            f"A_V={result.a_v:.2f}  v={result.v0_kms:.0f}  sig={result.sigma_kms:.0f}"
        )
        ax.legend(loc="upper right", fontsize=8)
        axr.plot(result.wavelength, result.flux - result.model, color="0.25", lw=0.55)
        axr.axhline(0.0, color="C3", lw=0.6)
        axr.set_xlabel("rest-frame vacuum wavelength (A)")
        axr.set_ylabel("residual")
        fig.tight_layout()
        args.plot.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.plot, dpi=140)
        plt.close(fig)
        print("wrote", args.plot)


if __name__ == "__main__":
    main()
