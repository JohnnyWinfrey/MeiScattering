import numpy as np
import MeiCoefficient as mc
import helper as eq

n, e, eta, m, k, a = eq.ExercisePart4()
n = n.astype(int)

C_abs_total = 0.0

""" need to plot over "a" with a given value not sum over the values of a"""
ni = n[0][50]
etai = eta[0][50]
ei = e[0][50]

psi_prime = mc.psi_prime(ni, etai)
psi = mc.psi_n(ni, etai)

cn1, cn2 = mc.cn(ni, m[0], ei, etai)


C_abs = 2 * np.pi / (abs(m[0]) ** 2 * k[0] ** 2) * sum((2 * ni + 1) * (1j * psi_prime * np.conj(psi)).real * (m[0] * abs(cn1[i]) ** 2 + np.conj(m[0]) * abs(cn2[i]) ** 2) for i in range(len(cn1)))


