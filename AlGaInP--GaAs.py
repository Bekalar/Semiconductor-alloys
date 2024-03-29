import matplotlib.pyplot as plt
import math as mt

y = 0.53
x = [i / 100.0 for i in range(0, 100)]
T = 300  # temperature in kelvin

# GaInP and AlInP

# gamma
GaP_gamma = 2.886 + 0.1081 * (1 - (mt.cosh(164 / T) * mt.sinh(164 / T)))  # a
AlP_gamma = 3.56  # a
InP_gamma = 1.4236
# L
GaP_L = 2.72  # a
AlP_L = 3.57  # a
InP_L = 2.014

# X
GaP_X = 2.35  # a
AlP_X = 2.52  # a
InP_X = 2.384 - 3.7 * pow(10, -4) * T


def calculate(GaP, AlP, InP):
    eg = []
    for i in x:
        eg.append(i * AlP + y * GaP + (1 - i - y) * InP)
    return eg


# AlGaInP
energy_gamma = calculate(GaP_gamma, AlP_gamma, InP_gamma)
energy_X = calculate(GaP_X, AlP_X, InP_X)
energy_L = calculate(GaP_L, AlP_L, InP_L)

fig, axes = plt.subplots(figsize=(10, 6))
axes.plot(x, energy_gamma, label="\u0393", color="black")
axes.plot(x, energy_X, label="X", color="red")
axes.plot(x, energy_L, label="L", color="blue")
axes.legend(fontsize=14)

plt.xlabel("Composition x")
plt.ylabel("Energy gap [eV]")
plt.title(r'$ Al_{x}Ga_{y}In_{1-x-y}P$', fontsize=20)
plt.grid()
plt.show()
