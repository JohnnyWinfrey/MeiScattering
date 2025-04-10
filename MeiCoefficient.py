import numpy as np
import scipy.special as sp
import helper as eq

def Pz(n, z):
    """Computes p_n(z) = j_n(z) / j_{n-1}(z) for n in [1, n_max]"""
    p = np.zeros(n + 1, dtype=complex)


    jn = sp.spherical_jn(n, z)
    jn_1 = sp.spherical_jn(n - 1, z)
    p = jn / jn_1

    return p

def Qz(n, e):
    """Computes q_n(z) using chi_n function instead of recursion"""
    q = chi_n(n-1, e)/chi_n(n-2, e)
    return q

def chi_n(n, e):
    """Computes χ_n explicitly using the definition"""
    return e*sp.spherical_yn(n+1, e)

def psi_n(n, e):
    """Computes ψ_n using p recursions and initial conditions"""
    psi_1 = (np.sin(e)/e) - np.cos(e)
    p_values = Pz(n, e)
    return e*sp.spherical_jn(n, e)#psi_1 * np.prod(p_values[1:n])

def A_n(n, z):
    """Computes A_n(z) using Pz"""
    pz = Pz(n, z)
    return -(n / z) + (1 / pz)

def B_n(n, e):
    """Computes B_n(z) using Qz"""
    qz = Qz(n, e)
    return -(n / e) + (1 / qz)

def Y1(n, e, eta, m):
    return A_n(n, eta) - (m*A_n(n,e))

def Y2(n, e, eta, m):
    return (m*A_n(n, eta)) - A_n(n,e)

def mie_coefficients(n_max, e, eta, m):
    """Computes Lorenz-Mie coefficients a_n1 and a_n2"""

    # Initialize array of zeros and establishes complex data type
    a_n1 = np.zeros(n_max, dtype=complex)
    a_n2 = np.zeros(n_max, dtype=complex)

    # Creating top row of variables for the print statement. This is just for debugging really.
    #print(f"\n{'n':>3} {'zp(n) Real':>20} {'zp(n) Imag':>20} {'ZA(n) Real':>20} {'ZA(n) Imag':>20} "
    #      f"{'a1(n) Real':>20} {'a1(n) Imag':>20} {'a2(n) Real':>20} {'a2(n) Imag':>20} {'psi(ka)':>20}{'chi(ka)':>20}")
    #print("=" * 210)



    # For loop to make calculations for each value of n up to n_max.
    for n in range(1, n_max + 1):
        pz = Pz(n, eta)
        psi = psi_n(n, e)
        chi = chi_n(n-1, e)
        A_eta = A_n(n, eta)
        A_e = A_n(n, e)
        B_e = B_n(n, e)

        # I had my chi and psi multiplied into A_eta and not the entire fraction ;_;
        a_n1[n - 1] = 1 / (1 + 1j * ((chi * (A_eta - m * B_e)) / (psi * (A_eta - m * A_e))))
        a_n2[n - 1] = 1 / (1 + 1j * ((chi * (m * A_eta - B_e)) / (psi * (m * A_eta - A_e))))

        # Print calculations
        #print(f"{n:3} {pz.real:20.10E} {pz.imag:20.10E} {A_eta.real:20.10E} {A_eta.imag:20.10E} "
        #      f"{a_n1[n - 1].real:20.10E} {a_n1[n - 1].imag:20.10E} {a_n2[n - 1].real:20.10E} {a_n2[n - 1].imag:20.10E} "
        #      f"{psi.real:20.10E} {chi:20.10E} ")

    return a_n1, a_n2 # Return Mei coefficients!
