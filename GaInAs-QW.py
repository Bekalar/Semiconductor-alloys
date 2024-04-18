import matplotlib.pyplot as plt

# gamma
C_bowing_gamma = 0.477  # c bowing
InAs_gamma = 0.417  # a
GaAs_gamma = 1.519

# x = [i / 100.0 for i in range(0, 100)]
x = [0.43]
p = x[0]
q = 1 - x[0]
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

# temperature parameters GaAs
alfa_GaAs = 0.0005405
beta_GaAs = 204

# temperature parameters InAs
alfa_InAs = 0.000276
beta_InAs = 93

# temperature
T = 300


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


def calculate_temperature(eg, alf, bet, Tem):
    temp = []
    for i in range(len(x)):
        temp.append(eg[i] - alf[i] * Tem * Tem / (Tem + bet[i]))
    return temp


def valence_band(vbo1, vbo2):
    val = []
    for i in x:
        val.append(vbo2 * i + vbo1 * (1 - i))
    return val


def conduction_band(valence, energy_band):
    con = []
    for i in range(len(x)):
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


# temperature interpolation
alfa = interpolate(alfa_GaAs, alfa_InAs)
beta = interpolate(beta_GaAs, beta_InAs)

# tension interpolation
a = interpolate(ao, a_layer)
b = interpolate(b_GaAs, b_InAs)
ac = interpolate(ac_GaAs, ac_InAs)
av = interpolate(av_GaAs, av_InAs)
c11 = interpolate(c11_GaAs, c11_InAs)
c12 = interpolate(c12_GaAs, c12_InAs)

eg = calculate(InAs_gamma, GaAs_gamma, C_bowing_gamma)
vb = valence_band(vbo_GaAs, vbo_InAs)
cb = conduction_band(vb, eg)

# temperatures
eg_T = calculate_temperature(eg, alfa, beta, T)
cb_T = conduction_band(vb, eg_T)
cb_E = calculate_tension_e(ac, a, ao, cb_T, c12, c11)

vb_Ehh = calculate_tension_ehh(av, a, ao, vb, c12, c11)
vb_Elh = calculate_tension_elh(av, a, ao, vb, c12, c11)

# substrate
vb_substrate = vbo_GaAs
cb_substrate = vbo_GaAs + GaAs_gamma

# medium
vb_medium = cb[0]
# vb_medium = cb_E[0]
cb_medium = vb_Ehh[0]

# width
width = 100  # 100nm
xl_sub = 0
xr_sub = width
xl_med = (100 / 2) - 10
xr_med = xl_med + 20

# plot
fig, axes = plt.subplots(figsize=(10, 6))
# vertical lines
axes.vlines(x=xl_med, ymin=vb_medium, ymax=cb_substrate, color="black")
axes.vlines(x=xl_med, ymin=vb_substrate, ymax=cb_medium, color="black")
axes.vlines(x=xr_med, ymin=vb_medium, ymax=cb_substrate, color="black")
axes.vlines(x=xr_med, ymin=vb_substrate, ymax=cb_medium, color="black")

# horizontal lines substrate
axes.hlines(y=vb_substrate, xmin=xl_sub, xmax=xl_med, color="black")
axes.hlines(y=vb_substrate, xmin=xr_med, xmax=xr_sub, color="black")
axes.hlines(y=cb_substrate, xmin=xl_sub, xmax=xl_med, color="black")
axes.hlines(y=cb_substrate, xmin=xr_med, xmax=xr_sub, color="black")

# horizontal lines medium
axes.hlines(y=cb, xmin=xl_med, xmax=xr_med, label=r"$ E_{c}$ without tension", color="red")
axes.hlines(y=vb, xmin=xl_med, xmax=xr_med, label=r"$ E_{v}$ without tension", color="red")
axes.hlines(y=cb_E, xmin=xl_med, xmax=xr_med, label=r"$ E_{c}$", color="black")
axes.hlines(y=vb_Ehh, xmin=xl_med, xmax=xr_med, label=r"$ E_{hh}$", color="black")
axes.hlines(y=vb_Elh, xmin=xl_med, xmax=xr_med, label=r"$ E_{lh}$", color="gray")

axes.legend(fontsize=14, loc='center right')

plt.xlabel("Width of quantum well [nm]")
plt.ylabel("Energy gap [eV]")
plt.title(r'$ In_{' + str(round(p, 2)) + '}Ga_{' + str(round(q, 2)) + '}As$ / ' + str(T) + 'K', fontsize=20)
plt.grid()
plt.show()
