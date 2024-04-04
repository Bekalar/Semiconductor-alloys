import matplotlib.pyplot as plt

# gamma
C_bowing_gamma = 0.477  # c bowing
InAs_gamma = 0.417  # a
GaAs_gamma = 1.519

x = [i / 100.0 for i in range(0, 100)]
a_layer = 6.0583  # InAs lattice constant
ao = 5.65325  # GaAs lattice constant

# parameters
ac_GaAs = -7.17
av_GaAs = -1.16
ac_InAs = -5.08
av_InAs = -1.00
b_GaAs = -2.0
b_InAs = -1.8
vbo_GaAs = -0.80
vbo_InAs = -0.59
c11_GaAs = 1221
c12_GaAs = 566
c11_InAs = 832.9
c12_InAs = 452.6


def interpolate(val1, val2):
    temp = []
    for i in x:
        temp.append(val1 * (1 - i) + val2 * i)
    return temp


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
    epsilon_x = [(a_o - a_a[i]) / a_a[i] for i in range(len(x))]
    epsilon_z = [-2 * (c_12[i] / c_11[i]) * epsilon_x[i] for i in range(len(x))]

    d_Ehc = [a_c[i] * (2 * epsilon_x[i] + epsilon_z[i]) for i in range(len(x))]

    Ec = [(energy_band[i] + d_Ehc[i]) for i in range(len(x))]

    return Ec


def calculate_tension_ehh(a_v, a_a, a_o, energy_band, c_12, c_11):
    epsilon_x = [(a_o - a_a[i]) / a_a[i] for i in range(len(x))]
    epsilon_z = [-2 * (c_12[i] / c_11[i]) * epsilon_x[i] for i in range(len(x))]

    d_Es = [-1 * b[i] * (epsilon_z[i] - epsilon_x[i]) for i in range(len(x))]
    d_Ehv = [a_v[i] * (2 * epsilon_x[i] + epsilon_z[i]) for i in range(len(x))]

    Ehh = [(energy_band[i] + d_Ehv[i] + d_Es[i]) for i in range(len(x))]

    return Ehh


def calculate_tension_elh(a_v, a_a, a_o, energy_band, c_12, c_11):
    epsilon_x = [(a_o - a_a[i]) / a_a[i] for i in range(len(x))]
    epsilon_z = [-2 * (c_12[i] / c_11[i]) * epsilon_x[i] for i in range(len(x))]

    d_Es = [-1 * b[i] * (epsilon_z[i] - epsilon_x[i]) for i in range(len(x))]
    d_Ehv = [a_v[i] * (2 * epsilon_x[i] + epsilon_z[i]) for i in range(len(x))]

    Elh = [(energy_band[i] + d_Ehv[i] - d_Es[i]) for i in range(len(x))]

    return Elh


a = interpolate(ao, a_layer)
b = interpolate(b_GaAs, b_InAs)
ac = interpolate(ac_GaAs, ac_InAs)
av = interpolate(av_GaAs, av_InAs)
c11 = interpolate(c11_GaAs, c11_InAs)
c12 = interpolate(c12_GaAs, c12_InAs)

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
