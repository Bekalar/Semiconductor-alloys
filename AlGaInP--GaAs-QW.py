import matplotlib.pyplot as plt
import math as mt

# GaInP and AlInP
# AlGaInP / GaAs

# temperature
T = 300

# gamma
GaP_gamma = 2.886 + 0.1081 * (1 - (mt.cosh(164 / T) / mt.sinh(164 / T)))  # a
AlP_gamma = 3.63  # a
InP_gamma = 1.4236

# x = [i / 100.0 for i in range(0, 100)]
y = [0.51]
x = [0.47]
p_x = x[0]
q_x = 1 - x[0]
p_y = y[0]
q_y = 1 - y[0]
ao = 5.65325  # GaAs lattice constant

# lattice constants
a_GaP = 5.4508
a_AlP = 5.4635
a_InP = 5.8690

# parameters
ac_GaP = -8.2
av_GaP = -1.7
ac_AlP = -5.7
av_AlP = -3.0
ac_InP = -6.0
av_InP = -0.6
b_GaP = -1.6
b_AlP = -1.5
b_InP = -2.0
vbo_GaP = -1.27
vbo_AlP = -1.74
vbo_InP = -0.94
c11_GaP = 1405
c12_GaP = 620.3
c11_AlP = 1330
c12_AlP = 630
c11_InP = 1011
c12_InP = 561

# GaAs parameters
vbo_GaAs = -0.80
GaAs_gamma = 1.519

# temperature parameters GaAs
alfa_GaAs = 0.0005405
beta_GaAs = 204

# temperature parameters GaP
alfa_GaP = 0.0005771
beta_GaP = 372

# temperature parameters AlP
alfa_AlP = 0.0005771
beta_AlP = 372

# temperature parameters InP
alfa_InP = 0.000363
beta_InP = 162


def interpolate(val1, val2, comp):
    temp = []
    for i in comp:
        temp.append(val1 * (1 - i) + val2 * i)
    return temp


def calculate(GaP, AlP, InP):
    eg = []
    eg_GaInP = y[0] * GaP + (1 - y[0]) * InP

    for i in x:
        eg.append(i * AlP + (1 - i) * eg_GaInP)

    return eg


def calculate_temperature(eg, alf, bet, Tem):
    temp = []
    for i in range(len(x)):
        temp.append(eg[i] - alf[i] * Tem * Tem / (Tem + bet[i]))
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
alfa_GaInP = interpolate(alfa_InP, alfa_GaP, y)
beta_GaInP = interpolate(beta_InP, beta_GaP, y)
alfa = interpolate(alfa_GaInP[0], alfa_AlP, x)
beta = interpolate(beta_GaInP[0], beta_AlP, x)

# tension interpolation
a_GaInP = interpolate(a_InP, a_GaP, y)
a_layer = interpolate(a_GaInP[0], a_AlP, x)
a = interpolate(ao, a_layer[0], x)
b_GaInP = interpolate(b_InP, b_GaP, y)
b = interpolate(b_GaInP[0], b_AlP, x)
ac_GaInP = interpolate(ac_InP, ac_GaP, y)
ac = interpolate(ac_GaInP[0], ac_AlP, x)
av_GaInP = interpolate(av_InP, av_GaP, y)
av = interpolate(av_GaInP[0], av_AlP, x)
c11_GaInP = interpolate(c11_InP, c11_GaP, y)
c11 = interpolate(c11_GaInP[0], c11_AlP, x)
c12_GaInP = interpolate(c12_InP, c12_GaP, y)
c12 = interpolate(c12_GaInP[0], c12_AlP, x)

eng = calculate(GaP_gamma, AlP_gamma, InP_gamma)
vb_GaInP = valence_band(vbo_InP, vbo_GaP)
vb = valence_band(vb_GaInP[0], vbo_AlP)
cb = conduction_band(vb, eng)

# temperatures tension substrate
eg_T = calculate_temperature(eng, alfa, beta, T)
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
substrate_text = 'AlGaInP ' + str(round((cb_substrate[0] - vb_substrate[0]), 2)) + ' eV'
axes.text(xl_med, 0.2, medium_text, fontsize=12)
axes.text(xl_sub, 0.2, substrate_text, fontsize=12)

plt.xlabel("Width of quantum well [nm]")
plt.ylabel("Energy gap [eV]")
plt.title(r'$ Al_{' + str(round(p_x, 2)) + '}(Ga_{' + str(round(p_y, 2)) + '}In_{' + str(round(q_y, 2)) + '})_{' + str(
    round(q_x, 2)) + '}P$ / GaAs ' + str(T) + 'K', fontsize=20)
plt.grid()
plt.show()
plt.show()
