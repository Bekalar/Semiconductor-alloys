import matplotlib.pyplot as plt

x = [i / 100.0 for i in range(0, 100)]
T = 300  # temperature in kelvin

# gamma
C_bowing_gamma = -0.48  # c
AlP_gamma = 3.56  # a
InP_gamma = 1.4236
# L
C_bowing_L = -0.19  # c bowing
AlP_L = 3.57  # a
InP_L = 2.014

# X
C_bowing_X = 0.38  # c bowing
AlP_X = 2.52  # a
InP_X = 2.384 - 3.7 * pow(10, -4) * T


def calculate(AlP, InP, C_bowing):
    a = AlP
    b = InP - AlP - C_bowing
    c = C_bowing

    eg = []
    for i in x:
        i = 1 - i
        eg.append(a + i * b + i * i * c)
    return eg


energy_gamma = calculate(AlP_gamma, InP_gamma, C_bowing_gamma)
energy_L = calculate(AlP_L, InP_L, C_bowing_L)
energy_X = calculate(AlP_X, InP_X, C_bowing_X)

fig, axes = plt.subplots(figsize=(10, 6))
axes.plot(x, energy_gamma, label="\u0393", color="black")
axes.plot(x, energy_X, label="X", color="red")
axes.plot(x, energy_L, label="L", color="blue")
axes.legend(fontsize=14)

plt.xlabel("Composition x")
plt.ylabel("Energy gap [eV]")
plt.title(r'$ Al_{x}In_{1-x}P$', fontsize=20)
plt.grid()
plt.show()
