from types import DynamicClassAttribute
import numpy as np
from scipy.integrate import solve_ivp
import pytest

import matplotlib.pyplot as plt
import sys

from reactor_ import point_kinetics_solver


def test_zero_reactivity_analytical():
    """Regression test: compares current power to a reference."""
    solver = point_kinetics_solver()
    solver.define_insertion()
    solver.load_kinetics_parameters(path="Parameter/Kinetics_parameter.csv")
    solver.solve()
    _, y = solver.get_solution()
    power = y[0]
    assert np.allclose(power, 1.0, rtol=1e-8)


def test_insertion_analytical():
    """Regression test: compares current power to a reference."""

    def P(t, rho, beta, Lambda, lam):
        return rho / (rho - beta) * np.exp((rho - beta) * t / Lambda) - beta / (
            rho - beta
        ) * np.exp(-lam * rho * t / (rho - beta))

    def step_func(t):
        if t < 0:
            return 0
        else:
            return 10

    solver = point_kinetics_solver()
    solver.load_kinetics_parameters(path="Parameter/one_group_test.csv")
    solver.define_insertion(step_func)
    solver.solve()
    t, y = solver.get_solution()

    # generate the analytical solution
    rho = 10
    beta = 670
    Lambda = 0.000432
    lam = 0.01334
    assert y[0] == pytest.approx(P(t, rho, beta, Lambda, lam))


def test_get_power():
    """Regression test: compares power evolution to a reference."""

    def P(t, rho, beta, Lambda, lam):
        return rho / (rho - beta) * np.exp((rho - beta) * t / Lambda) - beta / (
            rho - beta
        ) * np.exp(-lam * rho * t / (rho - beta))

    def step_func(t):
        if t < 0:
            return 0
        else:
            return 10

    solver = point_kinetics_solver()
    solver.load_kinetics_parameters(path="Parameter/one_group_test.csv")
    solver.define_insertion(step_func)
    solver.solve()
    rho = 10
    beta = 670
    Lambda = 0.000432
    lam = 0.01334

    assert solver.get_power() == pytest.approx(P(10, rho, beta, Lambda, lam))


def test_get_outlet_temperature():
    """Regression test: calculate Temperature rise to a solution."""

    solver = point_kinetics_solver()
    solver.set_initial_power(13000000)
    solver.set_inlet_temp(400)
    solver.load_kinetics_parameters(path="Parameter/one_group_test.csv")
    solver.define_insertion()
    solver.solve()
    assert solver.get_outlet_temp() == pytest.approx(400 + 43 + 1 / 3)


def test_get_reactor_temperature():
    """Regression test: calculate reactor Temperature to a solution."""

    solver = point_kinetics_solver()
    solver.set_initial_power(13000000)
    solver.set_inlet_temp(400)
    solver.load_kinetics_parameters(path="Parameter/one_group_test.csv")
    solver.define_insertion()
    solver.solve()
    assert solver.get_reactor_temp() == pytest.approx((400 + 400 + 43 + 1 / 3) / 2)


def test_dynamics():
    """Test: compares current power with dynamics to a value."""

    def insertion_func(t):
        if t < 1:
            return 0
        else:
            return 10

    solver = point_kinetics_solver()
    solver.set_initial_power(13000000)
    solver.load_kinetics_parameters(path="Parameter/one_group_test.csv")
    solver.set_dynamics()
    solver.define_insertion(insertion_func)
    solver.solve()
    assert solver.get_power() == pytest.approx(13213066.438692745)
