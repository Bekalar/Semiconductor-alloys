import matplotlib.pyplot as plt
import math as mt

y = [0.53]
x = [i / 100.0 for i in range(0, 100)]
T = 300  # temperature in kelvin
p = y[0]
q = 1 - y[0]
# GaInP and AlInP

# gamma
GaP_gamma = 2.886 + 0.1081 * (1 - (mt.cosh(164 / T) / mt.sinh(164 / T)))  # a
AlP_gamma = 3.63  # a
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

    eg_GaInP = y[0] * GaP + (1 - y[0]) * InP

    for i in x:
        eg.append(i * AlP + (1 - i) * eg_GaInP)

    return eg


# AlGaInP
energy_gamma = calculate(GaP_gamma, AlP_gamma, InP_gamma)
energy_X = calculate(GaP_X, AlP_X, InP_X)
energy_L = calculate(GaP_L, AlP_L, InP_L)

print(energy_gamma)

fig, axes = plt.subplots(figsize=(10, 6))
axes.plot(x, energy_gamma, label="\u0393", color="black")
axes.plot(x, energy_X, label="X", color="red")
axes.plot(x, energy_L, label="L", color="blue")
axes.legend(fontsize=14)

plt.xlabel("Composition x")
plt.ylabel("Energy gap [eV]")
plt.title(r'$ Al_{x}(Ga_{' + str(round(p, 2)) + '}In_{' + str(round(q, 2)) + '})_{1-x}P$', fontsize=20)
# plt.title(r'$ Al_{x}(Ga_{y}In_{1-y})_{1-x}P$', fontsize=20)
plt.grid()
plt.show()
