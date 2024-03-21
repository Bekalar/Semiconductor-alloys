import matplotlib.pyplot as plt

# gamma
C_bowing_gamma = 0.477  # c bowing
InAs_gamma = 0.417  # a
GaAs_gamma = 1.519

# Constans for tension
ac = -7.17
av = -1.16
ao = 5.65325  # Podłoże
a = 6.0583  # Warstwa
b = -2.0
c11 = 1221
c12 = 566

# valence band
vbo_GaAs = -0.80
vbo_InAs = -0.59

x = [i / 100.0 for i in range(0, 100)]


def calculate(InAs, GaAs, C_bowing):
    a1 = GaAs
    b1 = InAs - GaAs - C_bowing
    c1 = C_bowing

    Eg = []
    for i in x:
        Eg.append(a1 + i * b1 + i * i * c1)
    return Eg


def valence_band(vbo1, vbo2):
    val = []
    for i in x:
        val.append(vbo2 * i + vbo1 * (1 - i))
    return val


def conduction_band(valence, energy_band):
    con = []
    for i in range(0, 100):
        con.append(valence[i] + energy_band[i])
    return con


def calculate_tension_e(a_c, a_a, a_o, energy_band, c_12, c_11):
    epsilon_x = (a_o - a_a) / a_a
    epsilon_z = -2 * (c_12 / c_11) * epsilon_x

    d_Ehc = a_c * (2 * epsilon_x + epsilon_z)

    Ec = [(Eo + d_Ehc) for Eo in energy_band]

    return Ec


def calculate_tension_ehh(a_v, a_a, a_o, energy_band, c_12, c_11):
    epsilon_x = (a_o - a_a) / a_a
    epsilon_z = -2 * (c_12 / c_11) * epsilon_x

    d_Es = - b * (epsilon_z - epsilon_x)
    d_Ehv = a_v * (2 * epsilon_x + epsilon_z)

    Ehh = [(Eo + d_Ehv + d_Es) for Eo in energy_band]

    return Ehh


def calculate_tension_elh(a_v, a_a, a_o, energy_band, c_12, c_11):
    epsilon_x = (a_o - a_a) / a_a
    epsilon_z = -2 * (c_12 / c_11) * epsilon_x

    d_Es = - b * (epsilon_z - epsilon_x)
    d_Ehv = a_v * (2 * epsilon_x + epsilon_z)

    Elh = [(Eo + d_Ehv - d_Es) for Eo in energy_band]

    return Elh


eg = calculate(InAs_gamma, GaAs_gamma, C_bowing_gamma)

vb = valence_band(vbo_GaAs, vbo_InAs)
cb = conduction_band(vb, eg)

cb_E = calculate_tension_e(ac, a, ao, cb, c12, c11)

vb_Ehh = calculate_tension_ehh(av, a, ao, vb, c12, c11)

vb_Elh = calculate_tension_elh(av, a, ao, vb, c12, c11)

fig, axes = plt.subplots(figsize=(10, 6))
axes.plot(x, cb, label=r"$ E_{c}$ without tension", color="black")
axes.plot(x, vb, label=r"$ E_{v}$ without tension", color="black")

axes.plot(x, cb_E, label=r"$ E_{c}$", color="red")
axes.plot(x, vb_Ehh, '--', label=r"$ E_{hh}$", color="red")
axes.plot(x, vb_Elh, '-.', label=r"$ E_{lh}$", color="red")

axes.legend(fontsize=14)

plt.xlabel("Composition x")
plt.ylabel("Energy gap [eV]")
plt.title(r'$ In_{x}Ga_{1-x}As$', fontsize=20)
plt.grid()
plt.show()
