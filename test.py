import numpy as np
from scipy.integrate import solve_ivp
import pytest

import matplotlib.pyplot as plt
import sys

from reactor_ import point_kinetics_solver


def test_zero_reactivity_analytical():
    solver = point_kinetics_solver()

    def rho_zero(t):
        return t * 0.0

    solver.load_kinetics_parameters(path="Parameter/Kinetics_parameter.csv")
    solver.define_insertion(rho_zero)
    solver.solve()
    _, y = solver.get_solution()
    power = y[0]
    assert np.allclose(power, 1.0, rtol=1e-8)


def test_insertion_analytical():
    """Regression test: compares current solution to a saved reference."""

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
