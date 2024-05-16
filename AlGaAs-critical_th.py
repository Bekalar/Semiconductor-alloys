import matplotlib.pyplot as plt
import math as mt

# Matthews - Blakeslee model

# parameters
a_AlAs = 5.66139  # AlAs lattice constant
a_GaAs = 5.65325  # GaAs lattice substrate
c11_GaAs = 1221
c12_GaAs = 566
c11_AlAs = 1250
c12_AlAs = 534

# composition
x = [i / 100.0 for i in range(0, 100)]

# numerical values
granica = 500  # Angstrem
h = 0.00000001
n = int(granica / h)


def thickness_function():
    hc = []
    for i in range(1, n):
        hc[i + 1] = hc[i] * h
    return hc


def interpolate(val1, val2):
    temp = []
    for i in x:
        temp.append(val1 * (1 - i) + val2 * i)
    return temp


def critical_thickness(hc_fun):
    return 0


def bisec(xl, xr, tolx, fun, n_, h_):
    if fun(xr, n_, h_) * fun(xl, n_, h_) > 0:
        return "no zero"
    while (xr - xl) / 2.0 > tolx:
        xs = (xr + xl) / 2.0
        if fun(xr, n_, h_) * fun(xs, n_, h_) < 0:
            xl = xs
        else:
            xr = xs
    return (xr + xl) / 2.0


# calculate parameters
a = interpolate(a_GaAs, a_AlAs)
c11 = interpolate(c11_GaAs, c11_AlAs)
c12 = interpolate(c12_GaAs, c12_AlAs)

b = [a[i] / mt.sqrt(2) for i in range(len(x))]
v = [c12[i] / (c11[i] + c12[i]) for i in range(len(x))]
f = mt.fabs((a_GaAs - a_AlAs) / a_AlAs)