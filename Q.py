import MeiCoefficient as mc
import helper as eq
import numpy as np

n, e, eta, m, k, a= eq.fuller_values()



def Q_direct(n, e, eta, m, k, a):
    a_n1, a_n2 = mc.mie_coefficients(n, e, eta, m)

    Q_abs = (2 / e ** 2) * sum(
        ((2 * n + 1) / abs(mc.psi_n(n, e)) ** 2) * np.real(
            1j * mc.A_n(n, eta) * (
                    m * (a_n1[n - 1] / mc.Y1(n, e, eta, m)) * np.conj(a_n1[n - 1] / mc.Y1(n, e, eta, m)) +
                    np.conj(m) * abs(a_n2[n - 1] / mc.Y2(n, e, eta, m)) ** 2
            )
        )
        for n in range(1, n + 1)
    )
    return Q_abs


def Q_abs(n, e, eta, m, k, a):
    a_n1, a_n2 = mc.mie_coefficients(n, e, eta, m)
    C_sca = (2*np.pi/k**2) * sum((2*n + 1) * (abs(a_n1[n-1])**2 + abs(a_n2[n-1])**2) for n in range(1, n+1))
    C_ext = (2*np.pi/k**2) * sum((2*n + 1) * (a_n1[n-1].real + a_n2[n-1].real) for n in range(1, n+1))
    Q_abs = (C_ext - C_sca) / (np.pi * a**2)
    return Q_abs

def C_ext(n, e, eta, m, k, a):
    a_n1, a_n2 = mc.mie_coefficients(n, e, eta, m)
    C_ext = (2 * np.pi / k ** 2) * sum((2 * n + 1) * (a_n1[n - 1].real + a_n2[n - 1].real) for n in range(1, n + 1))
    return C_ext

def Q_backscatter(n, e, eta, m, k, a):
    a_n1, a_n2 = mc.mie_coefficients(n, e, eta, m)
    Q = abs((1/e)*(sum(((-1)**i)*(2*i+1)*(a_n1[i-1] - a_n2[i-1]) for i in range(1, n+1))))**2
    return Q

def Q_scattering(n, e, eta, m, k, a):
    a_n1, a_n2 = mc.mie_coefficients(n, e, eta, m)
    C_sca = (2*np.pi/k**2) * sum((2*n + 1) * (abs(a_n1[n-1])**2 + abs(a_n2[n-1])**2) for n in range(1, n+1))
    Q_sca = C_sca / (np.pi * a**2)
    return Q_sca

