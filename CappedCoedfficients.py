import MeiCoefficient as mc
import scipy.special as sp
import helper as eq
import numpy as np

def hankel(n, e):
    return sp.spherical_jn(n, e) + 1j*sp.spherical_yn(n, e)

def q_cap(n, e):
    hn1 = hankel(n, e)
    hn0 = hankel(n-1, e)
    q = hn1/hn0
    return q

def Cn(n, e):
    q = q_cap(n, e)
    C = -(n / e) + (1 / q)
    return C

def zeta_n(n, e):
    return e*hankel(n, e)

def zeta_prime(n, e):
    zeta = zeta_n(n, e)
    zeta_prime = Cn(n, e)*zeta
    return zeta_prime

def an1_1(n_max, e, m):
    eta = e*m
    a_n1 = np.zeros(n_max, dtype=complex)
    for n in range(1, n_max):
        num = zeta_n(n, eta)*(m*Cn(n, e)-Cn(n, eta))
        denom = mc.psi_n(n, eta)*(m*Cn(n, e)-mc.A_n(n, eta))
        a_n1[n] = -1*num/denom
    return a_n1

def an2_1(n, e, eta, m):

    num = zeta_n(n, eta)*m*Cn(n, eta)-Cn(n, e)
    denom = mc.psi_n(n, eta)*m*mc.A_n(n, eta)-Cn(n, e)
    return num/denom

def an1_2(n_max, e, m1, m2):
    a_n1 = np.zeros(n_max, dtype=complex)
    for n in range(1, n_max):
        m = m2/m1
        eta_plus = e*m1
        eta_minus = e*m2
        num = m*mc.psi_prime(n, eta_plus)*mc.psi_n(n, eta_minus) - mc.psi_n(n, eta_plus)*mc.psi_prime(n, eta_minus)
        denom = m*zeta_prime(n, eta_plus)*mc.psi_n(n, eta_minus) - zeta_n(n, eta_plus)*mc.psi_prime(n, eta_minus)
        a_n1[n] = -1*num/denom
    return a_n1

def an2_2(n, e, m1, m2):
    m = m2/m1
    eta_plus = e*m1
    eta_minus = e*m2
    num = m*mc.psi_n(n, eta_plus)*mc.psi_prime(n, eta_minus) - mc.psi_prime(n, eta_plus)*mc.psi_n(n, eta_minus)
    denom = m*zeta_n(n, eta_plus)*mc.psi_prime(n, eta_minus) - zeta_prime(n, eta_plus)*mc.psi_n(n, eta_minus)
    return -num/denom

n, e_host, e_core, wl, radius_host, radius_core, m_host, m_core = eq.ProjectValues()
an11 = an1_1(n, e_host, m_host)

print(an11)