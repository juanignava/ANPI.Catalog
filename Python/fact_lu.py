"""
Implementacion del metodo de factorizacion LU
"""

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

def zeros(n, m):
    R = []
    for i in range(0, n):
        R_aux = []
        for j in range(0, m):
            R_aux.append(0)
        R.append(R_aux)
    return R

def eye(n, m):
    R = []
    for i in range(0, n):
        R_aux = []
        for j in range(0, m):
            if i == j:
                R_aux.append(1)
            else:
                R_aux.append(0)
        R.append(R_aux)
    return R

def sust_adelante(A, b):

    m = len(A)
    x = zeros(m, 1)

    for i in range(0,m):
        aux = 0

        for j in range (0, i):
            aux += A[i][j]*x[j]

        x[i] = (1/A[i][i])*(b[i]-aux)

    return x

def sust_atras(A, b):

    m = len(A)
    x = zeros(m, 1)

    for i in range(m-1,-1,-1):
        aux = 0

        for j in range (i+1, m):
            aux += A[i][j]*x[j]

        x[i] = (1/A[i][i])*(b[i]-aux)

    return x

def resize(A, n, m):

    R = []
    for i in range(0, n):
        R_aux = []
        for j in range(0, m):
            R_aux.append(A[i][j])
        R.append(R_aux)
    return R

def sub_prin(A, tol):

    m = len(A)

    for i in range (1, m):
        Ak = A
        Ak.resize(i, i)
        deter = np.linalg.det(Ak)

        if deter < tol:
            return False

    return True

def calculo_lu(A):

    m = len(A)
    L = eye(m,m)
    U = zeros(m,m)

    for j in range(0, m):

        diag_el = A[j][j]
        for fila in range(j+1, m):
            multiplo = A[fila][j]/diag_el

            L[fila][j] = multiplo

            for column in range(0, m):
                A[fila][column] -= multiplo*A[j][column]

    return [L, A]


def fact_lu(A, b):

    l_u = calculo_lu(A)

    L = l_u[0]
    U = l_u[1]

    y = sust_adelante(L, b)

    x = sust_atras(U, y)

    return x

A = [[4, -2, 1], [20, -7, 12], [-8, 13, 17]]
b = [11, 70, 17]


#A = [[1, 1], [0, 3]]
#b = [2, 3]

#x = sust_atras(A, b)

x = fact_lu(A, b)
print(x)