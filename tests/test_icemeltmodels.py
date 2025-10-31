import pytest
import numpy as np

np.random.seed(42)
import icemeltmodels


def test_three_equation_melt_model_melt_rate_zero_current_zero_melt():
    """Melt rate should be zero if there is no current (no heat transfer)."""
    model = icemeltmodels.ThreeEquationMeltModel()
    N = 10
    current_speed = np.zeros(N)
    temperature = np.linspace(-2, 4, N)
    salinity = np.full(N, 35.0)
    pressure = np.zeros(N)

    np.testing.assert_allclose(
        model.melt_rate(current_speed, temperature, salinity, pressure),
        np.zeros(N),
        atol=1e-12,
    )


def test_two_equation_melt_model_neglecting_conduction_consistency():
    """Check that left hand side is equal to right hand side of TwoEquationMeltModelNeglectingConduction()"""
    model = icemeltmodels.TwoEquationMeltModelNeglectingConduction()
    N = 100

    # Generate random fields in appropriate numerical ranges.
    U = np.random.uniform(0, 5, N)
    T = np.random.uniform(-2, 4, N)
    S = np.random.uniform(30, 40, N)
    p = np.random.uniform(0, 5.0e6, N)

    # Compute left and right hand sides of heat conservation equation.
    left_hand_side = model.latent_heat_release(U, T, S, p)
    right_hand_side = model.ocean_heat_flux(U, T, S, p)

    np.testing.assert_allclose(left_hand_side, right_hand_side, atol=1e-12)


def test_two_equation_melt_model_consistency():
    """Check that left hand side is equal to right hand side of TwoEquationMeltModel()"""
    model = icemeltmodels.TwoEquationMeltModel()
    N = 100

    # Generate random fields in appropriate numerical ranges.
    U = np.random.uniform(0, 5, N)
    T = np.random.uniform(-2, 4, N)
    S = np.random.uniform(30, 40, N)
    p = np.random.uniform(0, 5.0e6, N)

    # Compute left and right hand sides of heat conservation equation.
    left_hand_side = model.latent_heat_release(U, T, S, p)
    right_hand_side = model.ocean_heat_flux(U, T, S, p) + model.ice_heat_flux(
        U, T, S, p
    )

    np.testing.assert_allclose(left_hand_side, right_hand_side, atol=1e-12)


def test_three_equation_melt_model_neglecting_conduction_consistency():
    """Check that left hand side is equal to right hand side of ThreeEquationMeltModelNeglectingConduction()"""
    model = icemeltmodels.ThreeEquationMeltModelNeglectingConduction()
    N = 100

    # Generate random fields in appropriate numerical ranges.
    U = np.random.uniform(0, 5, N)
    T = np.random.uniform(-2, 4, N)
    S = np.random.uniform(30, 40, N)
    p = np.random.uniform(0, 5.0e6, N)

    # Compute left and right hand sides of heat conservation equation.
    left_hand_side_heat = model.latent_heat_release(U, T, S, p)
    right_hand_side_heat = model.ocean_heat_flux(U, T, S, p)

    np.testing.assert_allclose(left_hand_side_heat, right_hand_side_heat, atol=1e-12)

    # Compute left and right hand sides of salt conservation equation.
    left_hand_side_salt = model.ice_salt_flux(U, T, S, p)
    right_hand_side_salt = model.ocean_salt_flux(U, T, S, p)

    np.testing.assert_allclose(left_hand_side_salt, right_hand_side_salt, atol=1e-12)


def test_three_equation_melt_model_consistency():
    """Check that left hand side is equal to right hand side of ThreeEquationMeltModel()"""
    model = icemeltmodels.ThreeEquationMeltModel()
    N = 100

    # Generate random fields in appropriate numerical ranges.
    U = np.random.uniform(0, 5, N)
    T = np.random.uniform(-2, 4, N)
    S = np.random.uniform(30, 40, N)
    p = np.random.uniform(0, 5.0e6, N)

    # Compute left and right hand sides of heat conservation equation.
    left_hand_side_heat = model.latent_heat_release(U, T, S, p)
    right_hand_side_heat = model.ocean_heat_flux(U, T, S, p) + model.ice_heat_flux(
        U, T, S, p
    )

    np.testing.assert_allclose(left_hand_side_heat, right_hand_side_heat, atol=1e-12)

    # Compute left and right hand sides of salt conservation equation.
    left_hand_side_salt = model.ice_salt_flux(U, T, S, p)
    right_hand_side_salt = model.ocean_salt_flux(U, T, S, p)

    np.testing.assert_allclose(left_hand_side_salt, right_hand_side_salt, atol=1e-12)
