import pytest
import numpy as np
from icemeltmodels import ThreeEquationMeltModel

def test_three_equation_melt_model_melt_rate_zero_current_zero_melt():
    """Melt rate should be zero if there is no current (no heat transfer)."""
    model = ThreeEquationMeltModel()
    N = 10
    current_speed = np.zeros(N)
    temperature = np.linspace(-2, 4, N)
    salinity = np.full(N, 35.0)
    pressure = np.zeros(N)

    np.testing.assert_allclose(
        model.melt_rate(current_speed, temperature, salinity, pressure),
        np.zeros(N),
        atol=1e-12
    )