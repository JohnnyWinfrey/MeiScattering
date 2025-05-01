import numpy as np
import pandas as pd

def n_max(k, a):
    """Computes maximum n for recursion based on particle properties."""
    rho = k * a
    return int(round(rho + (4.05 * rho ** (1 / 3) + 2)))

def fuller_values():
    # Example usage
    n_max = 10
    rho = 3.11400000000000  # Size parameter
    eta = 4.253724 + 1j * 1.557e-2
    m = 1.366 + 1j * 5e-3  # Refractive index
    a = 0.272584671033489
    k = rho / a
    return n_max, rho, eta, m, k, a

def conjugate(a):
    return a*np.conjugate(a)

def exercise_values():
    a = 5e-6
    m_real = 1.59
    m_imag = np.array([1e-6, 0.001, 0.1]) * 1j

    m = m_real + m_imag
    e = np.linspace(0.01, 100, 100)
    k = e / a

    eta = m[:, np.newaxis] * e

    n_values = np.array([n_max(k_i, a) for k_i in k])  # Shape: (1000,)

    return n_values, e, eta, m, k, a

def silicon_values():
    a = 5e-6
    with open('silicon_Schinke_2015.csv', 'r') as f:
        lines = f.readlines()

    split_line = None
    for i, line in enumerate(lines):
        if 'wl,k' in line:
            split_line = i
            break

    if split_line is None:
        raise ValueError("Could not find 'wl,k' section in the CSV.")

    lines_n = lines[:split_line]
    lines_k = lines[split_line:]

    data_n = pd.read_csv(pd.io.common.StringIO(''.join(lines_n)))
    wl = data_n['wl'].values
    wl = wl*1e-6
    n = data_n['n'].values

    data_k = pd.read_csv(pd.io.common.StringIO(''.join(lines_k)))
    k = data_k['k'].values

    wn = 2*np.pi / wl

    n_values = []

    for i in wn:
        n_values.append(n_max(i, a))

    return wl, n, k, a, n_values

def ExercisePart4():
    a = np.linspace(0.01e-6, 5e-6, 64)

    m = np.array([5.423 + 1j * 2.9078, 5.623 + 1j * 3.2627e-01, 3.931 + 1j * 1.8521e-02])
    wl = np.array([350e-9, 400e-9, 600e-9])

    k = 2 * np.pi / wl

    e = np.outer(k, a)
    eta = m[:, np.newaxis] * e

    n = np.empty((len(k), len(a)))

    for i in range(len(k)):
        for j in range(len(a)):
            n[i, j] = n_max(k[i], a[j])

    return n, e, eta, m, k, a

def ExercisePart5():
    a = 1e-6

    m = np.array([5.423 + 1j * 2.9078, 5.623 + 1j * 3.2627e-01, 3.931 + 1j * 1.8521e-02])
    wl = np.array([350e-9, 400e-9, 600e-9])

    k = 2 * np.pi / wl

    e = np.outer(k, a)
    eta = m[:, np.newaxis] * e

    n = []

    for i in range(len(k)):
        n.append(n_max(k[i], a))

    return n, e, eta, m, k, a, wl

def ProjectValues():
    e_host = 13.1213853164934
    e_core = 0.358141562509236
    wl = 3
    k = 2*np.pi/wl
    radius_host = 6.26500000000000
    radius_core = 0.171000000000000
    m_host = 1.40900000000000+1j*0.174700000000000
    m_core = 1.59000000000000+1j*0.660000000000000
    n = 32#n_max(k, radius_host)
    return n, e_host, e_core, wl, radius_host, radius_core, m_host, m_core