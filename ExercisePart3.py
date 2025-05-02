import helper as eq
import matplotlib.pyplot as plt

wl, index, extinction, a, n = eq.silicon_values()
reflectance = []
for n in index:
    r = (1-n)/(1+n)
    reflectance.append(r**2)

plt.figure()
plt.plot(wl, reflectance)
plt.show()