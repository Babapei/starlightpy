import numpy as np

from starlightpy import FitConfig, fit_spectrum
from starlightpy.simulate import default_absorption_bases, mock_observation


def test_mock_with_known_kinematics_still_recovers_x_and_av():
    """B1: fitter is told the true v, σ; it must still recover mix and dust."""
    wave = np.arange(3800.0, 5601.0, 1.0)
    bases = default_absorption_bases(wave)
    true_x = np.array([0.50, 0.35, 0.15])
    true_av = 0.30
    true_v = 100.0
    true_sig = 120.0
    config = FitConfig(
        a_v_bounds=(0.0, 1.0),
        a_v_step=0.05,
        v0_kms=true_v,
        sigma_kms=true_sig,
        search_kinematics=False,
    )
    flux, err = mock_observation(
        wave, bases, true_x, true_av, true_v, true_sig, config=config, snr=None
    )
    result = fit_spectrum(wave, flux, err, bases, config=config)
    assert abs(result.a_v - true_av) <= 0.05
    np.testing.assert_allclose(result.x_fraction, true_x, atol=0.10)
    assert result.chi2 < 1e-4
    assert result.v0_kms == true_v
    assert result.sigma_kms == true_sig
