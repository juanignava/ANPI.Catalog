"""
Implementación del método de la secante
"""

import string
from mpmath.libmp.libmpf import to_str
import sympy as sp 
import numpy as np
import matplotlib.pyplot as plt


def secante(f, x0, x1, tol, iterMax):
    x = sp.Symbol('x')
    f1 = sp.sympify(f)
    er = []
    err = tol + 1
    k = 0
    xk_0 = x0
    xk_1 = x1

    while err > tol and k < iterMax:
        k = k+1
        n = (xk_1 - xk_0)*(sp.N(f1.subs(x, xk_1)))
        d = sp.N(f1.subs(x, xk_1)) - sp.N(f1.subs(x, xk_0))
        xk_0 = xk_1
        xk_1 = xk_1 - n/d # realizar validación del denominador
        err = abs(f1.subs(x, xk_1))
        print("valor del x_k: " + str(xk_1))
        print("Valor del error: " + str(err))
        er.append(sp.N(err))

    return [xk_1, k, err]

f = 'exp(-x**2)-x'
x0 = 0
x1 = 1
tol = 10**-9
iterMax = 1000

y = secante(f, x0, x1, tol, iterMax)
print(y)