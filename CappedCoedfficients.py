import MeiCoefficient as mc
import scipy.special as sp
import helper as eq
import numpy as np

from MeiCoefficient import psi_prime, psi_n


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

def an2_1(n_max, e, m):
    eta = e * m
    a_n2 = np.zeros(n_max, dtype=complex)
    for n in range(1, n_max):
        num = zeta_n(n, eta)*(m*Cn(n, eta)-Cn(n, e))
        denom = mc.psi_n(n, eta)*(m*mc.A_n(n, eta)-Cn(n, e))
        a_n2[n] = -1*num/denom
    return a_n2

def an1_2(n_max, e_host, e_core, m1, m2):
    m = m2 / m1
    a_n1 = np.zeros(n_max + 1, dtype=complex)
    eta = e_host*m1
    for n in range(1, n_max + 1):

        num = mc.psi_n(n, eta)*zeta_prime(n, eta)-zeta_n(n, eta)*psi_prime(n, eta)
        denom = m1*zeta_prime(n, e_host)*mc.psi_n(n, eta)-zeta_n(n, e_host)*mc.psi_prime(n, eta)

        a_n1[n] = num/denom

    return a_n1


def an2_2(n_max, e_host, e_core, m1, m2):
    eta=e_host*m1
    a_n2 = np.zeros(n_max, dtype=complex)
    for n in range(1, n_max):
        num = mc.psi_n(n, eta)*zeta_prime(n, eta)-zeta_n(n, eta)*mc.psi_prime(n, eta)
        denom = m1*zeta_n(n, e_host)*mc.psi_prime(n, eta)-zeta_prime(n, e_host)*mc.psi_n(n, eta)
        a_n2[n] = -num/denom
    return a_n2


n, e_host, e_core, wl, radius_host, radius_core, m_host, m_core = eq.ProjectValues()
eta = e_host*m_host


a_cap_1n1 = an1_1(n, e_host, m_host)
a_cap_1n2 = an2_1(n, e_host, m_host)
c_cap_1n1 = an1_2(n, e_host, e_core, m_host, m_core)
c_cap_1n2 = an2_2(n, e_host, e_core, m_host, m_core)

a_1n1, a_1n2 = mc.mie_coefficients(n, e_host, eta, m_host)
a_1n1 = -a_1n1
a_1n2 = -a_1n2

c_1n1, c_1n2 = mc.cn(n, m_host, e_host, eta)
c_1n2 = -c_1n2

for i in range(1, n):
    print("n =", i)
    print("a_1n1 =", a_1n1[i-1])
    print("a_1n2 =", a_1n2[i-1])
    print("cn =", c_1n1[i])
    print("dn =", c_1n2[i])
    print("acap_11 =", a_cap_1n1[i])
    print("acap_21 =", a_cap_1n2[i])
    print("acap_12 =", c_cap_1n1[i])
    print("acap_22 =", c_cap_1n2[i])
    print()
