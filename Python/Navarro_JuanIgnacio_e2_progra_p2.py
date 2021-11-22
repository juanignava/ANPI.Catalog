"""
EXAMEN PARCIAL 2 - ANPI
Pregunta progrmación 2
Estudiante: Juan Ignacio Navarro Navarro (2019039662)
"""

from os import error
from matplotlib.pyplot import yticks
import numpy as np
from numpy import linalg
from numpy.core.fromnumeric import transpose
import sympy as sp
from sympy.core.evalf import N

def rayleigh(A, x0, iterMax, tol):

    """
    Metodo de Rayleigh para aproximar
    el valor propio mas grande de la 
    matriz A
    """

    # valores iniciales
    N = len(A[0])
    xk = np.matrix(x0)
    xk = np.transpose(xk)
    
    for i in range(iterMax):

        # calculo de sigma
        num = np.dot(transpose(xk), A)
        num = np.dot(num, xk)
        num = linalg.norm(num)
        den = np.dot(transpose(xk), xk)
        den = linalg.norm(den)
        sig_k_1 = num / den

        # calculo de y
        Im = np.eye(N)
        coef = np.subtract(A, sig_k_1*Im)
        yk = linalg.solve(coef, xk)

        # redefinir valor de x
        xk = yk / linalg.norm(yk)

        # calculo del error
        if i>0:
            error = abs(sig_k_1 - sig_k)

            if (error < tol):
                break

        sig_k = sig_k_1

    return sig_k

# Ejemplo numerico del enunciado

A = [[8, 1, 0, 0, 0],
     [1, 8, 1, 0, 0],
     [0, 1, 8, 1, 0],
     [0, 0, 1, 8, 1],
     [0, 0, 0, 1, 8]]

x0 = [1, 1, 1, 1, 1]

iterMax = 50
tol = 10**-12

val_max = rayleigh(A, x0, iterMax, tol)

print("El valor propio de mayor magnitud de A es:")
print(val_max)