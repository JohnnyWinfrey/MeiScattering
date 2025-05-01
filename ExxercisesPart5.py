import ScatteringMatrix as sm
import helper as eq
import numpy as np
import matplotlib.pyplot as plt

n, e, eta, m, k, a, wavelengths = eq.ExercisePart5()

theta_deg = np.arange(0, 180, 1)
theta_rad = np.cos(np.radians(theta_deg))

for idx, lam in enumerate(wavelengths):
    S11_list = []
    S12_list = []
    S33_list = []
    S34_list = []

    ni = n[idx]
    ei = e[idx]
    etai = eta[idx]
    mi = m[idx]

    for i in range(len(theta_rad)):
        print("Calculating Scattering Matrix for theta =", theta_rad[i])
        S11, S12, S33, S34 = sm.SMatrix(ni, ei, etai, mi, theta_rad[i])

        print(f"n = {ni}")
        print(f"e = {ei}")
        print(f"m = {mi}")
        print(f"eta = {etai}")
        print(f"theta = {theta_rad[i]}")
        print()

        S11_list.append(S11)
        S12_list.append(S12)
        S33_list.append(S33)
        S34_list.append(S34)

    S11_array = np.array(S11_list)
    S12_array = np.array(S12_list)
    S33_array = np.array(S33_list)
    S34_array = np.array(S34_list)

    S11_norm = S11_array / S11_array[0]
    S12_norm = S12_array / S11_array[0]
    S33_norm = S33_array / S11_array[0]
    S34_norm = S34_array / S11_array[0]

    plt.figure(figsize=(10,6))
    plt.plot(theta_deg, S11_norm, label='$S_{11}$')
    plt.plot(theta_deg, S12_norm, label='$S_{12}$')
    plt.plot(theta_deg, S33_norm, label='$S_{33}$')
    plt.plot(theta_deg, S34_norm, label='$S_{34}$')
    plt.xlabel('Theta')
    plt.ylabel('Normalized Scattering Matrix Elements')
    plt.title(f'Normalized Scattering Matrix Elements at λ = {lam * 1e9:.0f} nm')
    plt.legend()
    plt.grid(True)
    plt.show()
