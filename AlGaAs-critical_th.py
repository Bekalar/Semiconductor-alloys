import matplotlib.pyplot as plt
import math as mt
import numpy as np


# Matthews - Blakeslee model

def interpolate(val1, val2, x):
    temp = []
    for i in x:
        temp.append(val1 * (1 - i) + val2 * i)
    return temp


def hc_function(hc_val, b_val, v_val, f_val):
    temp = hc_val - ((b_val / (2 * mt.pi * f_val)) * ((1 - 0.25 * v_val) / (1 + v_val)) * (np.log(hc_val / b_val) + 1))
    return temp


def bisection(xl, xr, eps, fun, b_val, v_val, f_val):
    if fun(xr, b_val, v_val, f_val) * fun(xl, b_val, v_val, f_val) > 0:
        return None
    while (xr - xl) / 2.0 > eps:
        xs = (xr + xl) / 2.0
        if fun(xs, b_val, v_val, f_val) * fun(xl, b_val, v_val, f_val) <= 0:
            xr = xs
        else:
            xl = xs
    return (xr + xl) / 2.0


# Number of values
n = 100

# Parameters
a_AlAs = 0.566139  # AlAs lattice constant [nm]
a_GaAs = 0.565325  # GaAs lattice substrate [nm]
c11_GaAs = 1221
c12_GaAs = 566
c11_AlAs = 1250
c12_AlAs = 534

# Composition
x = np.linspace(0, 1, n)
# x = [i / 100.0 for i in range(0, 100)]

# [a,b]
p = 100  # a left value [nm]
q = 90000  # b right value [nm]

# Numerical tolerance parameter
epsilon = 0.00001

# Interpolate parameters
a = interpolate(a_GaAs, a_AlAs, x)
c11 = interpolate(c11_GaAs, c11_AlAs, x)
c12 = interpolate(c12_GaAs, c12_AlAs, x)

b = [a[i] / mt.sqrt(2) for i in range(len(x))]
v = [c12[i] / (c11[i] + c12[i]) for i in range(len(x))]
f = [abs((a_GaAs - a[i]) / a[i]) for i in range(len(x))]
f[0] = epsilon

# Plot figure
fig, (ax1, ax2) = plt.subplots(2, figsize=(12, 6))
fig.suptitle(r'$ Al_{x}Ga_{1-x}As$', fontsize=20)

# Calculate critical thickness for each composition
critical_thickness = np.zeros_like(x)

for i, x_val in enumerate(x):
    hc_values = np.linspace(p, q, n)
    hc_function_values = hc_function(hc_values, b[i], v[i], f[i])
    # print(hc_function_values)
    # Plot hc_fun(hc)
    ax1.plot(hc_values, hc_function_values)
    ax1.set_xlabel('hc [nm]')
    ax1.set_ylabel('hc_function [nm]')
    ax1.set_title('Critical Thickness function of hc [nm]')
    ax1.grid(True)

    zeros = bisection(p, q, epsilon, hc_function, b[i], v[i], f[i])
    # print(zeros)
    if zeros is not None:
        critical_thickness[i] = zeros

# Plot hc(x)
ax2.plot(x, critical_thickness)
ax2.set_ylim([0, 2000])
ax2.set_xlabel('Composition x')
ax2.set_ylabel('Critical Thickness [nm]')
ax2.set_title('Critical Thickness in Composition')
ax2.grid(True)

fig.tight_layout()
plt.show()
