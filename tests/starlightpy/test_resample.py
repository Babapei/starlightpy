import numpy as np
import pytest

from starlightpy import resample_to


def test_resample_matches_original_on_same_grid():
    wave = np.array([4000.0, 4010.0, 4020.0])
    flux = np.array([1.0, 2.0, 3.0])
    np.testing.assert_allclose(resample_to(wave, flux, wave), flux)


def test_resample_rejects_unsorted():
    with pytest.raises(ValueError):
        resample_to([4020.0, 4000.0], [1.0, 2.0], [4000.0, 4010.0])
