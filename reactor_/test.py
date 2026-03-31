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
    solver.load_kinetics_parameters(path="../Parameter/Kinetics_parameter.csv")
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
    solver.load_kinetics_parameters(path="../Parameter/one_group_test.csv")
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
    solver.load_kinetics_parameters(path="../Parameter/one_group_test.csv")
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
    solver.load_kinetics_parameters(path="../Parameter/one_group_test.csv")
    solver.define_insertion()
    solver.solve()
    assert solver.get_outlet_temp() == pytest.approx(400 + 43 + 1 / 3)


def test_get_reactor_temperature():
    """Regression test: calculate reactor Temperature to a solution."""

    solver = point_kinetics_solver()
    solver.set_initial_power(13000000)
    solver.set_inlet_temp(400)
    solver.load_kinetics_parameters(path="../Parameter/one_group_test.csv")
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
    solver.load_kinetics_parameters(path="../Parameter/one_group_test.csv")
    solver.set_dynamics()
    solver.define_insertion(insertion_func)
    solver.solve()
    assert solver.get_power() == pytest.approx(13213066.438692745)


def test_xenon_set_steady_state():
    """Test: compares set steady state Xenon level to a analytical solution."""

    solver = point_kinetics_solver()
    solver.set_initial_power(1)
    solver.define_insertion()
    solver.load_kinetics_parameters(path="../Parameter/one_group_test.csv")
    solver.load_Xenon_parameters(path="../Parameter/Xenon_parameter_test.csv")
    solver.set_initial_X_and_I()
    solver.solve(t_span=(0, 1))
    _, y = solver.get_solution(0)
    assert y[7] == pytest.approx(814.835440)


def test_xenon_steady_state():
    """Test: compares steady state Xenon level to a analytical solution."""

    solver = point_kinetics_solver()
    solver.set_initial_power(1)
    solver.define_insertion()
    solver.load_kinetics_parameters(path="../Parameter/one_group_test.csv")
    solver.load_Xenon_parameters(path="../Parameter/Xenon_parameter_test.csv")
    solver.solve(t_span=(0, 1000000))
    _, y = solver.get_solution(1000000)
    assert y[7] == pytest.approx(814.835440)


def test_xenon_shutdown():
    """Test: Tests Xenon level after shutdown."""

    def insertion_func(t):
        if t < 5:
            return 0
        else:
            return -1000

    solver = point_kinetics_solver()
    solver.set_initial_power(1e9)
    solver.define_insertion(func=insertion_func)
    solver.load_kinetics_parameters(path="../Parameter/one_group_test.csv")
    solver.load_Xenon_parameters(path="../Parameter/Xenon_parameter.csv")
    solver.set_initial_X_and_I()
    # test for 100h
    solver.solve(t_span=(0, 10 * 60 * 60))
    t, y = solver.get_solution(10 * 60 * 60)

    assert y[7] == pytest.approx(5.028635316084359e16)


def test_xenon_feedback():
    """Test: Reactor shutdown due to Xenon level."""
    solver = point_kinetics_solver()
    solver.set_initial_power(13e6)
    solver.define_insertion()
    solver.load_kinetics_parameters(path="../Parameter/one_group_test.csv")
    solver.load_Xenon_parameters(path="../Parameter/Xenon_parameter.csv")
    # solver.set_dynamics()
    solver.set_xenon()
    solver.solve(t_span=(0, 10 * 60 * 60))
    t, y = solver.get_solution(10 * 60 * 60)
    assert y[0] == pytest.approx(17.099465044339095)


def test_xenon_feedback_dynamic():
    """Test: Reactor shutdown due to Xenon level."""
    solver = point_kinetics_solver()
    solver.set_initial_power(3e9)
    solver.define_insertion()
    solver.load_kinetics_parameters(path="../Parameter/one_group_test.csv")
    solver.load_Xenon_parameters(path="../Parameter/Xenon_parameter.csv")
    solver.set_dynamics()
    solver.set_xenon()
    solver.solve(t_span=(0, 20 * 60 * 60))
    t, y = solver.get_solution(20 * 60 * 60)
    assert y[0] == pytest.approx(2.8606867655979525e-05)
