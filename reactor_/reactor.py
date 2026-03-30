import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


class PointKineticsSolver:
    def __init__(self):
        """
        Load kinetics parameters once at initialization.
        """
        self.alpha = None  # last value
        self.beta = None
        self.Lambda = None
        self.betas = None  # beta1 to beta6
        self.lambdas = None  # lambda1 to lambda6
        self.sol = None

    def _load_kinetics_parameters(self, path="../Parameter/Kinetics_parameter.csv"):
        """
        returns alpha, beta, Lambda,
        betas (as an array from 1-N)
        lambdas (as an array from 1-N)
        """
        # beta	beta1	beta2	beta3	beta4	beta5	beta6	Lambda	lambda1	lambda2	lambda3	lambda4	lambda5	lambda6	alpha
        data = np.loadtxt(
            "../Parameter/Kinetics_parameter.csv", delimiter=",", skiprows=1
        )

        self.alpha = data[-1]  # last value
        self.beta = data[0]
        self.Lambda = data[7]
        self.betas = data[1:7]  # beta1 to beta6
        self.lambdas = data[8:14]  # lambda1 to lambda6

    def calc_rho(self, t):
        # todo add insertion , temp feedback, and xenon
        if t < 5:
            return 0
        else:
            return 10

    def rhs(self, t, y):
        rho = self.calc_rho(t)
        P = y[0]
        C = y[1:]
        dPdt = (rho - self.beta) / self.Lambda * P + np.sum(C * self.lambdas)
        dCdt = self.betas / self.Lambda * P - self.lambdas * C
        return np.append(dPdt, dCdt)

    def solve(
        self,
        Initial_Power=1,
        t_span=(0, 10),
        n_points=1000,
        method="Radau",
        rtol=1e-6,
        atol=1e-8,
    ):
        # initial conditions
        Initial_C = (self.betas * Initial_Power) / (self.Lambda * self.lambdas)
        y0 = np.append(Initial_Power, Initial_C)

        t_eval = np.linspace(t_span[0], t_span[1], n_points)

        self.sol = solve_ivp(
            fun=self.rhs,
            t_span=t_span,
            y0=y0,
            method=method,
            rtol=rtol,
            atol=atol,
            t_eval=t_eval,
            dense_output=True,
        )
        if self.sol.success:
            print("Integration successful!")
        else:
            print("Integration failed:", self.sol.message)

        return self.sol

    def get_solution(self, t=None):
        """Return time and solution (interpolated if t is given)."""
        if self.sol is None:
            raise ValueError("You must call solve() first.")

        if t is None:
            return self.sol.t, self.sol.y
        else:
            return t, self.sol.sol(t)


# Solve it
Reactor_module = PointKineticsSolver()

Reactor_module._load_kinetics_parameters()
Reactor_module.solve()

t, y = Reactor_module.get_solution()


plt.plot(t, y[0], label="y1(t)")
# plt.plot(t, y[1], label="y2(t)")

plt.xlabel("t")
plt.ylabel("y")
plt.legend()
plt.grid(True)
plt.show()
