"""
implementacion del metodo del punto fijo
"""

from os import error
import sympy as sp 
import numpy as np
import matplotlib.pyplot as plt

def punto_fijo (f, x0, tol, iterMax):

    x = sp.Symbol('x')
    f1 = sp.sympify(f)
    k = 0
    xk = x0
    error = tol+1

    while k < iterMax:

        xk_1 = sp.N(f1.subs(x, xk))

        error = abs((xk_1-xk)/xk_1)
        xk = xk_1
        if error<tol:
            break
        k += 1 
    
    return [xk, k, error]

# Ejemplo numerico
f = "ln(2*x+1)"
x0 = 1.5
tol = 10**-5
iterMax = 1000

y = punto_fijo(f, x0, tol, iterMax)

print("xk: " + str(y[0]))
print("k: " + str(y[1]))
print("error: " + str(y[2]))