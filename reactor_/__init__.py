import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


class point_kinetics_solver:
    def __init__(
        self,
    ):
        """
        Load kinetics parameters once at initialization.
        """
        # Operational Parameter
        self.Initial_Power = 1  # Watts
        self.Power = self.Initial_Power
        # kinetics parameters
        self.alpha = None  # pcm / Kelvin
        self.beta = None  # pcm
        self.Lambda = None  # second
        self.betas = None  # beta1 to beta6
        self.lambdas = None  # lambda1 to lambda6 # 1/second
        # Temperature
        self.Inlet_temp = 373.15  # Kelvin
        self.Outlet_temp = 373.15
        self.c_pm_R = 300000  # Jules / Kelvin /second
        # opperation PARAMETERS
        self.insertion_func = None  # function
        self.dynamics = False

        # solution
        self.PK_solution = None

    def set_initial_power(self, set_Power):
        self.Initial_Power = set_Power
        self.Power = set_Power

    def load_kinetics_parameters(self, path="../Parameter/Kinetics_parameter.csv"):
        """
        returns alpha, beta, Lambda,
        betas (as an array from 1-N)
        lambdas (as an array from 1-N)
        """
        # beta	beta1	beta2	beta3	beta4	beta5	beta6	Lambda	lksjdföklajsdiklambda1	lambda2	lambda3	lambda4	lambda5	lambda6	alpha
        data = np.loadtxt(path, delimiter=",", skiprows=1)

        self.alpha = data[-1]  # last value
        self.beta = data[0]
        self.Lambda = data[7]
        self.betas = data[1:7]  # beta1 to beta6
        self.lambdas = data[8:14]  # lambda1 to lambda6

    def define_insertion(self, func=None):
        if func is None:

            def zero_rho(t):
                return 0 * t

            self.insertion_func = zero_rho

        else:
            self.insertion_func = func

    def get_rho(self, t):
        # todo add insertion , temp feedback, and xenon
        if self.dynamics:
            return self.insertion_func(t)
        else:
            return self.insertion_func(t)

    def calc_outlet_temperature(self):
        self.Outlet_temp = self.Inlet_temp + self.Power / self.c_pm_R

    def rhs(self, t, y):
        rho = self.get_rho(t)
        self.Power = y[0]
        C = y[1:]
        dPdt = (rho - self.beta) / self.Lambda * self.Power + np.sum(C * self.lambdas)
        dCdt = self.betas / self.Lambda * self.Power - self.lambdas * C
        return np.append(dPdt, dCdt)

    def solve(
        self,
        t_span=(0, 10),
        n_points=1000,
        method="Radau",
        rtol=1e-6,
        atol=1e-8,
    ):
        # initial conditions
        Initial_C = (self.betas * self.Initial_Power) / (self.Lambda * self.lambdas)
        y0 = np.append(self.Initial_Power, Initial_C)

        t_eval = np.linspace(t_span[0], t_span[1], n_points)

        self.PK_solution = solve_ivp(
            fun=self.rhs,
            t_span=t_span,
            y0=y0,
            method=method,
            rtol=rtol,
            atol=atol,
            t_eval=t_eval,
            dense_output=True,
        )
        return self.PK_solution

    def get_solution(self, t=None):
        """Return time and solution (interpolated if t is given)."""
        if self.PK_solution is None:
            raise ValueError("You must call solve() first.")

        if t is None:
            return self.PK_solution.t, self.PK_solution.y
        else:
            return t, self.PK_solution.sol(t)

    def get_power(self):
        """Return time and solution (interpolated if t is given)."""
        if self.PK_solution is None:
            raise ValueError("You must call solve() first.")

        return self.Power

    def get_outlet_temp(self):
        """Return time and solution (interpolated if t is given)."""
        self.calc_outlet_temperature()
        if self.Outlet_temp is None:
            raise ValueError("You must call solve() first.")

        return self.Outlet_temp
