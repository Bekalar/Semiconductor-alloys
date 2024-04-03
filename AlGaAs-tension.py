import matplotlib.pyplot as plt

# gamma
C_bowing_gamma = [(-0.127 + 1.310 * i / 100.0) for i in range(0, 100)]  # c linear bowing parameter
AlAs_gamma = 3.099  # a
GaAs_gamma = 1.519

x = [i / 100.0 for i in range(0, 100)]
a_layer = 5.66139  # AlAs lattice constant
ao = 5.65330  # GaAs lattice constant

# parameters
ac_GaAs = -7.17
av_GaAs = -1.16
ac_AlAs = -5.64
av_AlAs = -2.47
b_GaAs = -2.0
b_AlAs = -2.3
vbo_GaAs = -0.80
vbo_AlAs = -1.33
c11_GaAs = 1221
c12_GaAs = 566
c11_AlAs = 1250
c12_AlAs = 534


def interpolate(val1, val2):
    temp = []
    for i in x:
        temp.append(val1 * (1 - i) + val2 * i)
    return temp


def calculate(AlAs, GaAs, C_bowing):
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


a = interpolate(a_layer, ao)
b = interpolate(b_GaAs, b_AlAs)
ac = interpolate(ac_GaAs, ac_AlAs)
av = interpolate(av_GaAs, av_AlAs)
c11 = interpolate(c11_GaAs, c11_AlAs)
c12 = interpolate(c12_GaAs, c12_AlAs)

eg = calculate(AlAs_gamma, GaAs_gamma, C_bowing_gamma)

vb = valence_band(vbo_GaAs, vbo_AlAs)
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
plt.title(r'$ Al_{x}Ga_{1-x}As$', fontsize=20)
plt.grid()
plt.show()
