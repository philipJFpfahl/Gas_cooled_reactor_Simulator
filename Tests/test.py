import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


def P(t, rho, beta, Lambda, lam):
    return rho / (rho - beta) * np.exp((rho - beta) * t / Lambda) - beta / (
        rho - beta
    ) * np.exp(-lam * rho * t / (rho - beta))


def C(t, rho, beta, Lambda, lam):
    return rho * beta / (rho - beta) ** 2 * np.exp((rho - beta) * t / Lambda) + beta / (
        lam * Lambda
    ) * np.exp(-lam * rho * t / (rho - beta))
