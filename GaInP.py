import matplotlib.pyplot as plt
import math as mt

x = [i / 100.0 for i in range(0, 100)]
T = 300  # temperature in kelvin

# gamma
C_bowing_gamma = 0.65  # c
GaP_gamma = 2.886 + 0.1081 * (1-(mt.cosh(164/T)*mt.sinh(164/T)))  # a
InP_gamma = 1.4236
# L
C_bowing_L = 1.03  # c bowing
GaP_L = 2.72  # a
InP_L = 2.014

# X
C_bowing_X = 0.20  # c bowing
GaP_X = 2.35  # a
InP_X = 2.384 - 3.7 * pow(10, -4) * T


def calculate(GaP, InP, C_bowing):
    a = GaP
    b = InP - GaP - C_bowing
    c = C_bowing

    eg = []
    for i in x:
        i = 1 - i
        eg.append(a + i * b + i * i * c)
    return eg


energy_gamma = calculate(GaP_gamma, InP_gamma, C_bowing_gamma)
energy_L = calculate(GaP_L, InP_L, C_bowing_L)
energy_X = calculate(GaP_X, InP_X, C_bowing_X)

fig, axes = plt.subplots(figsize=(10, 6))
axes.plot(x, energy_gamma, label="\u0393", color="black")
axes.plot(x, energy_X, label="X", color="red")
axes.plot(x, energy_L, label="L", color="blue")
axes.legend(fontsize=14)

plt.xlabel("Composition x")
plt.ylabel("Energy gap [eV]")
plt.title(r'$ Ga_{x}In_{1-x}P$', fontsize=20)
plt.grid()
plt.show()
