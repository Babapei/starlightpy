import numpy as np

from starlightpy import FitConfig, fit_spectrum, get_extinction_curve
from starlightpy.model import build_model, normalize_at, normalize_bases


def _planck_like(wave, temperature):
    wl_cm = np.asarray(wave) * 1e-8
    c2 = 1.4388
    return 1.0 / wl_cm**5 * np.exp(-c2 / (wl_cm * temperature))


def test_noiseless_recovery_of_mixture_and_av():
    wave = np.arange(3800.0, 7001.0, 2.0)
    bases = np.column_stack(
        [
            _planck_like(wave, 4000),
            _planck_like(wave, 6000),
            _planck_like(wave, 9000),
        ]
    )
    true_x = np.array([0.55, 0.30, 0.15])
    true_av = 0.40
    q = get_extinction_curve(wave, law="CCM")
    q0 = float(np.median(q[(wave >= 4010) & (wave <= 4060)]))
    bases_n = normalize_bases(wave, bases, (4010.0, 4060.0))
    obs, _ = normalize_at(wave, build_model(true_x, bases_n, q, true_av, q_lambda0=q0))
    err = np.full_like(obs, 0.01)

    result = fit_spectrum(wave, obs, err, bases, config=FitConfig(a_v_bounds=(0.0, 1.0), a_v_step=0.05))
    assert abs(result.a_v - true_av) <= 0.05
    np.testing.assert_allclose(result.x_fraction, true_x, atol=0.08)
    assert result.chi2 < 1e-6
