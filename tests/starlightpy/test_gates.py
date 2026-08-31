import numpy as np
import pytest

from starlightpy import FitConfig, fit_spectrum
from starlightpy.clip import clip_and_refit
from starlightpy.optimize import metropolis_anneal


def test_phase_gates_are_closed():
    with pytest.raises(NotImplementedError):
        metropolis_anneal()
    with pytest.raises(NotImplementedError):
        clip_and_refit()
    with pytest.raises(NotImplementedError):
        fit_spectrum(
            np.arange(4010.0, 4065.0),
            np.ones(55),
            np.ones(55),
            np.ones((55, 2)),
            config=FitConfig(search_kinematics=True),
        )
