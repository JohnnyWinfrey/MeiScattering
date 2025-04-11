import Q
import helper as eq
import matplotlib.pyplot as plt
import numpy as np
import concurrent.futures
import os
from tqdm import tqdm

def Q_absorption(n, e, eta, m, k, a):
    Q_graph = []
    for i in tqdm(range(len(e)), desc=f"Absorption", position=0, leave=False):
        Q_graph.append(Q.Q_abs(n[i], e[i], eta[i], m[i], k[i], a))
    return Q_graph, e

def Q_backscatter(n, e, eta, m, k, a):
    Q_graph = []
    for i in tqdm(range(len(e)), desc=f"Backscatter", position=1, leave=False):
        Q_graph.append(Q.Q_backscatter(n[i], e[i], eta[i], m[i], k[i], a))
    return Q_graph, e

def Q_scattering(n, e, eta, m, k, a):
    Q_graph = []
    for i in tqdm(range(len(e)), desc=f"Scattering", position=2, leave=False):
        Q_graph.append(Q.Q_scattering(n[i], e[i], eta[i], m[i], k[i], a))
    return Q_graph, e

if __name__ == "__main__":
    wl, index, extinction, a, n = eq.silicon_values()

    wn = 2 * np.pi / wl
    e = wn * a
    m = index + (1j * extinction)
    eta = e * m

    print("Size of e =", len(e))
    print("Size of eta =", len(eta))
    print("Size of m =", len(m))
    print("Size of n =", n)

    Q_a, e = Q_absorption(n, e, eta, m, wn, a)
    Q_b, e = Q_backscatter(n, e, eta, m, wn, a)
    Q_s, e = Q_scattering(n, e, eta, m, wn, a)

    plt.figure()
    plt.plot(e, Q_a, label="Q_abs")
    plt.plot(e, Q_b, label="Q_back")
    plt.plot(e, Q_s, label="Q_scat")
    plt.title("Absorption Efficiency")
    plt.xlabel("e")
    plt.ylabel("Q")
    plt.legend()
    plt.show()
