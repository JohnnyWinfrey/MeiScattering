import numpy as np
import MeiCoefficient as mc
import helper as eq
import tau
import Q as q

n, e, eta, m, k, a= eq.fuller_values()

def SMatrix(n, e, eta, m, theta):
    an1, an2 = mc.mie_coefficients(n, e, eta, m)
    tau1, tau2 = tau.tau_mnp(1, n, theta)

    S1 = 1j/k * sum(((2*n+1)/(n*(n+1)) * ((an1[n-1]*tau1[n])-an2[n-1]*tau2[n]))
                    for n in range(1, n+1))

    S2 = 1j/k * sum(((2*n+1)/(n*(n+1)) * ((-1*an1[n-1]*tau2[n])+an2[n-1]*tau1[n]))
                    for n in range(1, n+1))

    C_ext_S = 4*np.pi / k**2 * S1.real
    C_ext_direct = q.C_ext(n, e, eta, m, k, a) #* (np.pi * a**2)

    #print("C_ext from S1 =", C_ext_S)
    #print("C_ext from direct expansion =", C_ext_direct)
    S11 = (1/2)*(abs(S2)**2 + abs(S1)**2)
    S12 = (1/2)*(abs(S2)**2 - abs(S1)**2)
    S33 = (1/2)*(np.conjugate(S2)*S1 + S2*np.conjugate(S1))
    S34 = (1j/2)*(S1*np.conjugate(S2) - S2*np.conjugate(S1))

    #print()
    #print(f"a_n1 = {an1[0]}")
    #print(f"a_n2 = {an2[0]}")
    #print()
    #print("S1 =", S1)
    #print("S2 =", S2)
    #print()
    print("S_11 = %.2f" % S11)
    print("S_12 = %.2f" % S12)
    print("S_33 = %.2f" % S33.real)
    print("S_34 = %.2f" % S34.real)

    return S11, S12, S33, S34

