import matplotlib.pyplot as plt

# gamma
C_bowing_gamma = 0.477  # c bowing
InAs_gamma = 0.417  # a
GaAs_gamma = 1.519

# L
C_bowing_L = 0.33  # c bowing
InAs_L = 1.133  # a
GaAs_L = 1.815

# X
C_bowing_X = 1.4  # c bowing
InAs_X = 1.433  # a
GaAs_X = 1.981

x = [i / 100.0 for i in range(0, 100)]


def calculate(InAs, GaAs, C_bowing):
    a = InAs
    b = GaAs - InAs - C_bowing
    c = C_bowing

    eg = []
    for i in x:
        i = 1 - i
        eg.append(a + i * b + i * i * c)
    return eg


energy_gamma = calculate(InAs_gamma, GaAs_gamma, C_bowing_gamma)
energy_L = calculate(InAs_L, GaAs_L, C_bowing_L)
energy_X = calculate(InAs_X, GaAs_X, C_bowing_X)

fig, axes = plt.subplots(figsize=(10, 6))
axes.plot(x, energy_gamma, label="\u0393", color="black")
axes.plot(x, energy_X, label="X", color="red")
axes.plot(x, energy_L, label="L", color="blue")
axes.legend(fontsize=14)

plt.xlabel("Composition x")
plt.ylabel("Energy gap [eV]")
plt.title(r'$ Ga_{x}In_{1-x}As$', fontsize=20)
plt.grid()
plt.show()
