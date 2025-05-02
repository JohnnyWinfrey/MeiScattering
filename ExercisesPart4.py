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

pi = np.pi
r = np.linspace(0.1e-6, 5e-6, 500)  # radius in meters

wavelengths = np.array([350e-9, 400e-9, 600e-9])
results = {}

for idx, lam in enumerate(wavelengths):
    print(f"Computing for lambda = {lam}")
    m_value = m_list[idx]
    k = 2 * pi / lam # wavenumber for a given wavelength
    C_a_prime_list = [] # initializing absorption cross-section array
    ka_m = k*r # calculating rho for a given wavenumber
    eta_m = m_value*ka_m
    nn = int(ka_m[-1] + 4 * ka_m[-1]**(1/3) + 2) # n_max at the boundary of the particle
    cn1_array, cn2_array = mc.cn(nn, m_value, ka_m[-1], eta_m[-1]) # c_n at the boundary

    for kr in ka_m:
        print(f"Computing for kr = {kr:.3e} / {max(ka_m):.3e}")

        n_array = np.arange(1, nn + 1)

        eta = m_value * kr

        psi_prime_values = psi_prime_n(nn, eta)
        psi_star = np.conj(psi_n(nn, eta))

        terms = (2 * n_array + 1) * np.real(
            1j * psi_prime_values * psi_star * (
                m_value * np.abs(cn1_array)**2 + np.conj(m_value) * np.abs(cn2_array)**2
            )
        )

        sum_term = np.sum(terms)
        C_a = (2 * pi) / (np.abs(m_value)**2 * k**2) * sum_term
        C_a_prime_list.append(C_a)

    results[lam] = C_a_prime_list

# Plotting
plt.figure(0)
plt.plot(ka_m, results[wavelengths[0]], label=f'λ = {wavelengths[0] * 1e9:.0f} nm')
plt.xlabel('kr')
plt.ylabel("C'_a")
plt.title("Absorption Cross Section $C'_a$ vs $kr$")
plt.legend()
plt.grid(True)

plt.figure(1)
plt.plot(ka_m, results[wavelengths[1]], label=f'λ = {wavelengths[1] * 1e9:.0f} nm')
plt.xlabel('kr')
plt.ylabel("C'_a")
plt.title("Absorption Cross Section $C'_a$ vs $kr$")
plt.legend()
plt.grid(True)

plt.figure(2)
plt.plot(ka_m, results[wavelengths[2]], label=f'λ = {wavelengths[2] * 1e9:.0f} nm')
plt.xlabel('kr')
plt.ylabel("C'_a")
plt.title("Absorption Cross Section $C'_a$ vs $kr$")
plt.legend()
plt.grid(True)
plt.show()
