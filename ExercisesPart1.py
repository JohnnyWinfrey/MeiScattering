import Q
import helper as eq
import matplotlib.pyplot as plt
import numpy as np
import concurrent.futures
import os
from tqdm import tqdm

def Q_absorption(n, e, eta, m, k, a):
    Q_graph = []
    for i in tqdm(range(len(e)), desc=f"Absorption m={m}", position=0, leave=False):
        Q_graph.append(Q.Q_abs(n[i], e[i], eta[i], m, k[i], a))
    return Q_graph, e

def Q_backscatter(n, e, eta, m, k, a):
    Q_graph = []
    for i in tqdm(range(len(e)), desc=f"Backscatter m={m}", position=1, leave=False):
        Q_graph.append(Q.Q_backscatter(n[i], e[i], eta[i], m, k[i], a))
    return Q_graph, e

def Q_scattering(n, e, eta, m, k, a):
    Q_graph = []
    for i in tqdm(range(len(e)), desc=f"Scattering m={m}", position=2, leave=False):
        Q_graph.append(Q.Q_scattering(n[i], e[i], eta[i], m, k[i], a))
    return Q_graph, e

def compute_Qs(i, n, e, eta, m_values, k, a):
    m = m_values[i]
    eta_m = eta[i]
    Qa, _ = Q_absorption(n, e, eta_m, m, k, a)
    Qb, _ = Q_backscatter(n, e, eta_m, m, k, a)
    Qs, _ = Q_scattering(n, e, eta_m, m, k, a)
    return Qa, Qb, Qs

def compute_Qs_wrapper(i_and_args):
    i, n, e, eta, m_values, k, a = i_and_args
    return compute_Qs(i, n, e, eta, m_values, k, a)


if __name__ == "__main__":
    n, e, eta, m_values, k, a = eq.exercise_values()

    print("Size of n =", len(n))
    print("Size of e =", len(e))
    print("Size of eta =", len(eta[0]))
    print("Size of m =", len(m_values))
    print("Size of k =", len(k))
    print("Size of a =", a)

    inputs = [(i, n, e, eta, m_values, k, a) for i in range(len(m_values))]

    Q_abs = []
    Q_back = []
    Q_scat = []

    max_workers = os.cpu_count()

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = list(tqdm(executor.map(compute_Qs_wrapper, inputs), total=len(inputs), desc="Processing"))

    for Qa, Qb, Qs in results:
        Q_abs.append(Qa)
        Q_back.append(Qb)
        Q_scat.append(Qs)

    plt.figure(0)
    for i, Qa in enumerate(Q_abs):
        plt.plot(e, Qa, label=f"m = {m_values[i]}")
    plt.title("Absorption Efficiency")
    plt.xlabel("e")
    plt.ylabel("Q_absorption")
    plt.legend()

    plt.figure(1)
    for i, Qs in enumerate(Q_scat):
        plt.plot(e, Qs, label=f"m = {m_values[i]}")
    plt.title("Scattering Efficiency")
    plt.xlabel("e")
    plt.ylabel("Q_scattering")
    plt.legend()

    plt.figure(2)
    for i, Qb in enumerate(Q_back):
        plt.plot(e, Qb, label=f"m = {m_values[i]}")
    plt.title("Backscattering Efficiency")
    plt.xlabel("e")
    plt.ylabel("Q_backscattering")
    plt.legend()

    plt.show()

