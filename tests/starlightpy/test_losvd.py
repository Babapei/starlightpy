import numpy as np

from starlightpy import apply_losvd

C_KMS = 299792.458


def test_losvd_identity_when_zero():
    wave = np.linspace(4000, 5000, 50)
    flux = np.linspace(1.0, 2.0, wave.size)
    np.testing.assert_allclose(apply_losvd(wave, flux, 0.0, 0.0), flux)


def test_positive_velocity_moves_absorption_to_the_red():
    wave = np.arange(4800.0, 5200.0, 0.5)
    center = 5000.0
    flux = 1.0 - 0.8 * np.exp(-0.5 * ((wave - center) / 2.0) ** 2)
    v = 300.0
    shifted = apply_losvd(wave, flux, v, 0.0)
    expected = center * np.exp(v / C_KMS)
    assert abs(wave[np.argmin(shifted)] - expected) < 1.5
