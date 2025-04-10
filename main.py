import numpy as np
import MeiCoefficient as mc
import helper as eq
import ScatteringMatrix as sm
import Q

"""_________________System Parameters_________________"""
angle = 0
theta = np.cos(np.radians(angle)) # cos of incident angle
n, e, eta, m, k, a = eq.fuller_values()
p = 2

"""_________________Q_________________"""
Q_d = Q.Q_direct(n, e, eta, m, k, a)
Q_a = Q.Q_abs(n, e, eta, m, k, a)

print("Q_abs =", Q_a)
print("Q_direct =", Q_d)
print()

"""_________________Scattering Matrix_________________"""

