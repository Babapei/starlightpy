import numpy as np
import pytest

from starlightpy import (
    light_to_mass,
    light_weighted_age,
    light_weighted_metallicity,
)


def test_g6_light_to_mass_fractions():
    x = np.array([0.5, 0.5])
    mass, mu = light_to_mass(x, mass_to_light=np.array([1.0, 3.0]))
    np.testing.assert_allclose(mass, [0.5, 1.5])
    np.testing.assert_allclose(mu, [0.25, 0.75])


def test_g6_missing_ml_raises():
    with pytest.raises(ValueError, match="mass_to_light"):
        light_to_mass(np.array([0.5, 0.5]), mass_to_light=None)


def test_g6_light_weighted_age():
    x = np.array([0.5, 0.5])
    ages = np.array([1.0, 3.0])
    assert light_weighted_age(x, ages) == pytest.approx(2.0)
    with pytest.raises(ValueError, match="ages"):
        light_weighted_age(x, None)


def test_g6_light_weighted_z_requires_metadata():
    x = np.array([0.2, 0.8])
    z = np.array([0.004, 0.02])
    assert light_weighted_metallicity(x, z) == pytest.approx(0.2 * 0.004 + 0.8 * 0.02)
    with pytest.raises(ValueError, match="metallicities"):
        light_weighted_metallicity(x, None)
