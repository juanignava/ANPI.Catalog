"""
Implementación del método de la falsa posicion
"""

import sympy as sp 
import numpy as np
import matplotlib.pyplot as plt


def falsa_posicion(f, x0, x1, tol, iterMax):
    x = sp.Symbol('x')
    f1 = sp.sympify(f)
    er = []
    err = tol + 1
    k = 0
    xk_0 = x0
    xk_1 = x1

    if (sp.N(f1.subs(x, xk_0)) * sp.N(f1.subs(x, xk_1)) >= 0):
        xk_1 = 'NA'
        err = 'NA'
        print ("No se cumple el teorema de Bolzano")
        return [xk_1, k, err]

    while err > tol and k < iterMax:
        k = k+1
        n = (xk_1 - xk_0)*(sp.N(f1.subs(x, xk_1)))
        d = sp.N(f1.subs(x, xk_1)) - sp.N(f1.subs(x, xk_0))
        temp = xk_1 - n/d

        if (sp.N(f1.subs(x, xk_0)) * temp < 0):
            xk_1 = temp
        else:
            xk_0 = temp

        err = abs(f1.subs(x, xk_1))
        print("valor del x_k: " + str(xk_1))
        print("Valor del error: " + str(err))
        er.append(sp.N(err))

    return [xk_1, k, err]

f = 'exp(-x)-ln(x)'
x0 = 1
x1 = 2
tol = 10**-9
iterMax = 1000

y = falsa_posicion(f, x0, x1, tol, iterMax)
print("xk: " + str(y[0]))
print("k: " + str(y[1]))
print("error: " + str(y[2]))