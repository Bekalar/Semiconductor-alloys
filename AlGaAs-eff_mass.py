import matplotlib.pyplot as plt

# Luttinger parameters GaAs
gamma1_GaAs = 6.98
gamma2_GaAs = 2.06
gamma3_GaAs = 2.93
electron_mass_GaAs = 0.067

# Luttinger parameters AlAs
gamma1_AlAs = 3.76
gamma2_AlAs = 0.82
gamma3_AlAs = 1.42
electron_mass_AlAs = 0.15

# composition
x = [i / 100.0 for i in range(0, 100)]


def interpolate(val1, val2):
    temp = []
    for i in x:
        temp.append(val1 * (1 - i) + val2 * i)
    return temp


def hh_zdirection(gamma_1, gamma_2):
    hh = []
    for i in range(len(x)):
        hh.append(1 / (gamma_1[i] - 2 * gamma_2[i]))
    return hh


def lh_zdirection(gamma_1, gamma_2):
    lh = []
    for i in range(len(x)):
        lh.append(1 / (gamma_1[i] + 2 * gamma_2[i]))
    return lh


gamma1 = interpolate(gamma1_GaAs, gamma1_AlAs)
gamma2 = interpolate(gamma2_GaAs, gamma2_AlAs)
gamma3 = interpolate(gamma3_GaAs, gamma3_AlAs)
electron_mass = interpolate(electron_mass_GaAs, electron_mass_AlAs)
hh_mass = hh_zdirection(gamma1, gamma2)
lh_mass = lh_zdirection(gamma1, gamma2)

fig, axes = plt.subplots(figsize=(10, 6))
axes.plot(x, electron_mass, label="Electron mass", color="black")
axes.plot(x, hh_mass, label="Heavy holes mass", color="red")
axes.plot(x, lh_mass, label="Light holes mass", color="blue")
axes.legend(fontsize=14)

plt.xlabel("Composition x")
plt.ylabel("Effective mass")
plt.title(r'$ Al_{x}Ga_{1-x}As$ Particle effective mass ', fontsize=20)
plt.grid()
plt.show()
