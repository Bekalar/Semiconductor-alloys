import matplotlib.pyplot as plt

x = [i / 100.0 for i in range(0, 100)]

# gamma
C_bowing_gamma = [(-0.127 + 1.310 * i / 100.0) for i in range(0, 100)]  # c linear bowing parameter
AlAs_gamma = 3.099  # a
GaAs_gamma = 1.519
# L
C_bowing_L = 0  # c bowing
AlAs_L = 2.46  # a
GaAs_L = 1.815

# X
C_bowing_X = 0.055  # c bowing
AlAs_X = 2.24  # a
GaAs_X = 1.981


def calculate(AlAs, GaAs, C_bowing):
    a = AlAs
    b = GaAs - AlAs - C_bowing
    c = C_bowing

    eg = []
    for i in x:
        i = 1 - i
        eg.append(a + i * b + i * i * c)
    return eg


def calculate_list(AlAs, GaAs, C_bowing):
    a = []
    b = []
    c = []

    for i in range(0, 100):
        a.append(AlAs)
        b.append(GaAs - AlAs - C_bowing[i])
        c.append(C_bowing[i])

    eg = []
    for i in range(0, 100):
        j = 1 - i / 100.0
        eg.append(a[i] + j * b[i] + j * j * c[i])
    return eg


energy_gamma = calculate_list(AlAs_gamma, GaAs_gamma, C_bowing_gamma)
energy_L = calculate(AlAs_L, GaAs_L, C_bowing_L)
energy_X = calculate(AlAs_X, GaAs_X, C_bowing_X)

fig, axes = plt.subplots(figsize=(10, 6))
axes.plot(x, energy_gamma, label="\u0393", color="black")
axes.plot(x, energy_X, label="X", color="red")
axes.plot(x, energy_L, label="L", color="blue")
axes.legend(fontsize=14)

plt.xlabel("Composition x")
plt.ylabel("Energy gap [eV]")
plt.title(r'$ Al_{x}Ga_{1-x}As$', fontsize=20)
plt.grid()
plt.show()
