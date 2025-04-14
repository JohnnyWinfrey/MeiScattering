import numpy as np
import matplotlib.pyplot as plt
import MeiCoefficient as mc

pi = np.pi
r = 5e-6

# Complex refractive indices for each wavelength
m_list = np.array([
    5.423 + 1j * 2.9078,          # for 350 nm
    5.623 + 1j * 3.2627e-01,      # for 400 nm
    3.931 + 1j * 1.8521e-02       # for 600 nm
])

# Wavelengths
wavelengths = np.array([350e-9, 400e-9, 600e-9])

# kr range: from 0.1e-6 to ka + 1
n_points = 10
ka = 2 * pi * r / min(wavelengths)  # max kr for smallest wavelength
kr_values = np.linspace(0.1e-6, ka + 1, n_points)

results = {}

# Loop over wavelengths
for idx, lam in enumerate(wavelengths):
    m_value = m_list[idx]
    k = 2 * pi / lam
    C_a_prime = []

    for kr in kr_values:
        print(f"Printing kr={kr} out of {max(kr_values)}")

        n_max = int(kr + 4 * kr**(1/3) + 2)

        if n_max < 1:  # handle very small kr
            C_a_prime.append(0)
            continue

        n_array = np.arange(1, n_max + 1)

        eta = m_value * kr

        psi_prime = mc.psi_prime(n_max, kr)  # Does not return array up to n_max. Need to figure these out.
        psi_star = np.conj(mc.psi_n(n_max, kr))  # Does not return array up to n_max
        cn1_array, cn2_array = mc.cn(n_max, m_value, kr, eta)  # return arrays up to n_max

        terms = (2 * n_array + 1) * np.real(
            1j * psi_prime * psi_star * (
                m_value * np.abs(cn1_array)**2 + np.conj(m_value) * np.abs(cn2_array)**2
            )
        )

        sum_term = np.sum(terms)

        C_a = (2 * pi) / (np.abs(m_value)**2 * k**2) * sum_term
        C_a_prime.append(C_a)

    results[lam] = C_a_prime

# Plotting
plt.figure(figsize=(8,6))
for lam in wavelengths:
    plt.plot(kr_values, results[lam], label=f'λ = {lam * 1e9:.0f} nm')

plt.xlabel('kr')
plt.ylabel("Q'_a")
plt.title("Absorption Cross Section $C'_a$ vs $kr$")
plt.legend()
plt.grid(True)
plt.show()
