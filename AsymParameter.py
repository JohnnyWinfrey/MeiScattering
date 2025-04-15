import numpy as np
import matplotlib.pyplot as plt
import MeiCoefficient as mc
import scipy.special as sp
import helper as eq
import Q as q

n, e, eta, m, k, a = eq.fuller_values()

C_sca = q.Q_scattering(n, e, eta, m, k, a) * (np.pi*a**2)

# Load Lorenz-Mie coefficients from your module
a_n1, a_n2 = mc.mie_coefficients(n, e, eta, m)

# Precompute arrays for the summation
n_array = np.arange(1, n + 1)

# ---- First Summation Term ---- #
# n(n + 2)/(n + 1) * (a_n1 * conj(a_(n+1)1) + a_n2 * conj(a_(n+1)2))
# note: arrays must be shifted for a_(n+1)
a_n1_shift = np.append(a_n1[1:], 0)  # shift left, pad last with zero
a_n2_shift = np.append(a_n2[1:], 0)

first_sum_terms = (
    n_array * (n_array + 2) / (n_array + 1) *
    (a_n1 * np.conj(a_n1_shift) + a_n2 * np.conj(a_n2_shift))
)

first_sum = np.sum(first_sum_terms)

# ---- Second Summation Term ---- #
# (2n + 1) / [n(n + 1)] * Re(a_n1 * conj(a_n2))
second_sum_terms = (
    (2 * n_array + 1) / (n_array * (n_array + 1)) *
    np.real(a_n1 * np.conj(a_n2))
)

second_sum = np.sum(second_sum_terms)

# ---- Total Summation ---- #
total_sum = first_sum + second_sum

# ---- Final Asymmetry Parameter Calculation ---- #
asymmetry = (4 * np.pi) / (k**2 * C_sca) * total_sum
print(f"C_scat = {C_sca}")
print(f"<cos> = {asymmetry}")
