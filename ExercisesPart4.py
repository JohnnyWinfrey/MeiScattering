from string import printable

import numpy as np
import matplotlib.pyplot as plt
import MeiCoefficient as mc
import scipy.special as sp
"""NOTES"""
"""I need to fix that ka is set to it's maximum value for all values of lambda and not for each value of lambda. 
    I also need to figure out how to handle psi prime and psi_n since they do not return arrays.
    There's also the issue of putting an if statement in for when the radius exits the sphere. i.e. setting m
    equal to 1+0j since it would be in air. """

pi = np.pi
r = 5e-5  # radius in meters

def psi_n(n_max, e):
    """Computes psi_n = e * j_n(e) for n from 1 to n_max."""
    psi = []
    for n in range(1, n_max + 1):
        psi.append(e * sp.spherical_jn(n, e))
    return np.array(psi)

def Pz(n_max, z):
    """Computes p_n(z) = j_n(z) / j_{n-1}(z) for n from 1 to n_max."""
    p = []
    for n in range(1, n_max + 1):
        jn = sp.spherical_jn(n, z)
        jn_1 = sp.spherical_jn(n - 1, z)
        if jn_1 == 0:
            p.append(1)  # Avoid division by zero
        else:
            p.append(jn / jn_1)
    return np.array(p)

def A_n(n_max, z):
    """Computes A_n(z) = -n/z + 1/p_n(z) for n from 1 to n_max."""
    pz = Pz(n_max, z)
    A = []
    for n in range(1, n_max + 1):
        if pz[n - 1] == 0:
            A.append(np.inf)  # Avoid division by zero
        else:
            A.append(-(n / z) + (1 / pz[n - 1]))
    return np.array(A)

def psi_prime_n(n_max, e):
    """Computes psi'_n = A_n * psi_n for n from 1 to n_max."""
    A = A_n(n_max, e)
    psi = psi_n(n_max, e)
    return A * psi

# Complex refractive indices for each wavelength
m_list = np.array([
    5.423 + 1j * 2.9078,         # for 350 nm
    5.623 + 1j * 0.32627,        # for 400 nm
    3.931 + 1j * 0.018521        # for 600 nm
])


wavelengths = np.array([350e-9, 400e-9, 600e-9])

n_points = 100
ka = 2 * pi * r / min(wavelengths)  # max kr for smallest wavelength
kr_values = np.linspace(1e-6, ka + 1, n_points)

results = {}

# Loop over wavelengths
for idx, lam in enumerate(wavelengths):
    m_value = m_list[idx]
    k = 2 * pi / lam
    C_a_prime_list = []

    for kr in kr_values:
        print(f"Computing for kr = {kr:.3e} / {max(kr_values):.3e}")

        n_max = int(kr + 4 * kr**(1/3) + 2)

        if n_max < 1:  # handle very small kr
            C_a_prime_list.append(0)
            continue

        n_array = np.arange(1, n_max + 1)

        eta = m_value * kr

        psi_prime_values = psi_prime_n(n_max, kr)
        psi_star = np.conj(psi_n(n_max, kr))

        cn1_array, cn2_array = mc.cn(n_max, m_value, kr, eta)

        terms = (2 * n_array + 1) * np.real(
            1j * psi_prime_values * psi_star * (
                m_value * np.abs(cn1_array)**2 + np.conj(m_value) * np.abs(cn2_array)**2
            )
        )

        sum_term = np.sum(terms)

        C_a = (2 * pi) / (np.abs(m_value)**2 * k**2) * sum_term
        C_a_prime_list.append(C_a)

    results[lam] = C_a_prime_list
    print(C_a_prime_list)

# Plotting
plt.figure(figsize=(8, 6))
for lam in wavelengths:
    plt.plot(kr_values, results[lam], label=f'λ = {lam * 1e9:.0f} nm')

plt.xlabel('kr')
plt.ylabel("Q'_a")
plt.title("Absorption Cross Section $C'_a$ vs $kr$")
plt.legend()
plt.grid(True)
plt.show()
