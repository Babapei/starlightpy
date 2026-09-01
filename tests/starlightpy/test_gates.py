import numpy as np
import pytest

from starlightpy import FitConfig, fit_spectrum
from starlightpy.optimize import metropolis_anneal


def test_phase_gates_anneal_still_closed_clip_is_open():
    with pytest.raises(NotImplementedError):
        metropolis_anneal()
    wave = np.arange(4010.0, 4065.0)
    result = fit_spectrum(
        wave,
        np.ones(wave.size),
        np.ones(wave.size),
        np.ones((wave.size, 2)),
        config=FitConfig(clip_nsigma=3.0, a_v_bounds=(0.0, 0.1), a_v_step=0.1),
    )
    assert result.n_clipped >= 0
    assert result.a_v >= 0.0
