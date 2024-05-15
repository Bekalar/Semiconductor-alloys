import matplotlib.pyplot as plt

# substrate
# medium
# substrate

# gamma
C_bowing_gamma = [(-0.127 + 1.310 * i / 100.0) for i in range(0, 100)]  # c linear bowing parameter
GaAs_gamma = 1.519
AlAs_gamma = 3.099  # a

# x = [i / 100.0 for i in range(0, 100)]
x = [0.47]
p = x[0]
q = 1 - x[0]
a_layer = 5.66139  # AlAs lattice constant
ao = 5.65325  # GaAs lattice constant

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

# temperature parameters GaAs
alfa_GaAs = 0.0005405
beta_GaAs = 204

# temperature parameters InAs
alfa_AlAs = 0.000885
beta_AlAs = 530

# temperature
T = 300


def interpolate(val1, val2):
    temp = []
    for i in x:
        temp.append(val1 * (1 - i) + val2 * i)
    return temp


def calculate(AlAs, GaAs, C_bowing, x):
    a = []
    b = []
    c = []

    for i in range(len(x)):
        a.append(AlAs)
        b.append(GaAs - AlAs - C_bowing[i])
        c.append(C_bowing[i])

    eg = []
    for i in range(len(x)):
        eg.append(a[i] + x[0] * b[i] + x[0] * x[0] * c[i])
    return eg


def calculate_temperature(eg, alf, bet, Tem):
    temp = []
    for i in range(len(x)):
        temp.append(eg[i] - (alf[i] * Tem * Tem) / (Tem + bet[i]))
    return temp


def calculate_temperature_substrate(egs, alfs, bets, Tems):
    temps = (egs - (alfs * Tems * Tems) / (Tems + bets))
    return temps


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
alfa = interpolate(alfa_GaAs, alfa_AlAs)
beta = interpolate(beta_GaAs, beta_AlAs)

# tension interpolation
a = interpolate(ao, a_layer)
b = interpolate(b_GaAs, b_AlAs)
ac = interpolate(ac_GaAs, ac_AlAs)
av = interpolate(av_GaAs, av_AlAs)
c11 = interpolate(c11_GaAs, c11_AlAs)
c12 = interpolate(c12_GaAs, c12_AlAs)

eg = calculate(AlAs_gamma, GaAs_gamma, C_bowing_gamma, x)
vb = valence_band(vbo_GaAs, vbo_AlAs)
cb = conduction_band(vb, eg)

# temperatures tension substrate
eg_T = calculate_temperature(eg, alfa, beta, T)
cb_T = conduction_band(vb, eg_T)
cb_E = calculate_tension_e(ac, a, ao, cb_T, c12, c11)

vb_Ehh = calculate_tension_ehh(av, a, ao, vb, c12, c11)
vb_Elh = calculate_tension_elh(av, a, ao, vb, c12, c11)

if cb <= cb_E:
    cb_substrate = cb
else:
    cb_substrate = cb_E

if vb_Ehh >= vb_Elh:
    vb_substrate = vb_Ehh
else:
    vb_substrate = vb_Elh

# temperature tension medium
eg_T_GaAs = calculate_temperature_substrate(GaAs_gamma, alfa_GaAs, beta_GaAs, T)
part = eg_T_GaAs / 4  # moving valence band
vb_medium = vbo_GaAs + part
cb_medium = vb_medium + eg_T_GaAs

# width
width = 100  # 100nm
well_width = 10  # quantum well width
xl_sub = 0
xr_sub = width
xl_med = (width / 2) - well_width / 2
xr_med = xl_med + well_width

# plot
fig, axes = plt.subplots(figsize=(10, 6))
# vertical lines
axes.vlines(x=xl_med, ymin=vb_medium, ymax=vb_substrate, color="black")
axes.vlines(x=xl_med, ymin=cb_medium, ymax=cb_substrate, color="black")
axes.vlines(x=xr_med, ymin=vb_medium, ymax=vb_substrate, color="black")
axes.vlines(x=xr_med, ymin=cb_medium, ymax=cb_substrate, color="black")

# horizontal lines substrate
axes.hlines(y=vb_Ehh, xmin=xl_sub, xmax=xl_med, label=r"$ E_{hh}$", color="red")
axes.hlines(y=vb_Elh, xmin=xl_sub, xmax=xl_med, label=r"$ E_{lh}$", color="blue")
axes.hlines(y=vb_Ehh, xmin=xr_med, xmax=xr_sub, color="red")
axes.hlines(y=vb_Elh, xmin=xr_med, xmax=xr_sub, color="blue")
axes.hlines(y=cb_E, xmin=xl_sub, xmax=xl_med, label=r"$ E_{c}$", color="violet")
axes.hlines(y=cb_E, xmin=xr_med, xmax=xr_sub, color="violet")

# without tension substrate
axes.hlines(y=vb, xmin=xl_sub, xmax=xl_med, label=r"$ E_{v}$ no tension", color="grey")
axes.hlines(y=vb, xmin=xr_med, xmax=xr_sub, color="grey")
axes.hlines(y=cb, xmin=xl_sub, xmax=xl_med, label=r"$ E_{c}$ no tension", color="grey")
axes.hlines(y=cb, xmin=xr_med, xmax=xr_sub, color="grey")

# horizontal lines substrate
axes.hlines(y=cb_medium, xmin=xl_med, xmax=xr_med, color="black")
axes.hlines(y=vb_medium, xmin=xl_med, xmax=xr_med, color="black")

axes.legend(fontsize=14, loc='center right')
medium_text = 'GaAs ' + str(round((cb_medium - vb_medium), 2)) + ' eV'
substrate_text = 'AlGaAs ' + str(round((cb_substrate[0] - vb_substrate[0]), 2)) + ' eV'
axes.text(xl_med, 0.2, medium_text, fontsize=12)
axes.text(xl_sub, 0.2, substrate_text, fontsize=12)

plt.xlabel("Width of quantum well [nm]")
plt.ylabel("Energy gap [eV]")
plt.title(r'$ Al_{' + str(round(p, 2)) + '}Ga_{' + str(round(q, 2)) + '}As$ / GaAs ' + str(T) + 'K', fontsize=20)
plt.grid()
plt.show()
