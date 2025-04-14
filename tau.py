import numpy as np


def tau_mnp(m, n_max, x):
    """Computes τ_mnp using recursion up to n_max."""

    # Array to store tau values (initialize with zeros)
    tau1 = np.zeros(n_max + 2)
    tau2 = np.zeros(n_max + 2)

    # Set initial conditions explicitly
    if m == 0:
        tau2[0] = 0
        tau1[0] = 0
        tau1[1] = np.sqrt(1 - x ** 2)
    elif m == 1:
        tau2[1] = 1
        tau1[1] = -x


    # Compute τ recursively
    for n in range(2, n_max + 1):
        tau2[n] = ((2 * n - 1) / (n - m)) * x * tau2[n - 1] - ((n + m - 1) / (n - m)) * tau2[n - 2]

    for n in range(2, n_max + 1):
        tau1[n] = (((n+m)/m)*tau2[n-1]) - (n*x/m * tau2[n])

    return tau1, tau2  # Returns τ values


