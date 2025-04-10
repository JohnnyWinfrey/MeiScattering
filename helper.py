import numpy as np
import csv

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
    a = 5e-6  # Particle radius
    m_real = 1.59  # Real part of refractive index
    m_imag = np.array([1e-6, 0.001, 0.1]) * 1j  # Imaginary part as a NumPy array

    m = m_real + m_imag  # Complex refractive indices
    e = np.linspace(0.01, 100, 100)  # Size parameter range
    k = e / a  # Wavenumber

    # Compute eta for each m and each e (shape: (3, 1000))
    eta = m[:, np.newaxis] * e  # Broadcasting to get all combinations

    # Compute n_max for each k (assuming n_max is scalar for each k)
    n_values = np.array([n_max(k_i, a) for k_i in k])  # Shape: (1000,)

    return n_values, e, eta, m, k, a

def silicon_values():
    Col1 = "wl"
    Col2 = "n"
    dictionary={Col1:[], Col2:[]}
    csvFile = csv.reader(open("silicon_Schinke_2015.csv", "rb"))
    for row in csvFile:
        dictionary[Col1].append(row[0])
        dictionary[Col2].append(row[1])

    return dictionary



