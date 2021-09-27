# Segundo ejericico del examen parcial
# Análisis Numérico para Ingeniería
# Estudiante: Juan Ignacio Navarro Navarro (2019039662)

import numpy as np

### Funciones Auxiliares ###
def eye(n, m):
    """
    Metodo auxiliar para contruir una matriz identidad
    de tamaño n x m
    """
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

def resize(A, n, m):
    """
    Metodo auxiliar para cortar una matriz en su
    n-ésima fila y m-ésima columna. Como salida se tiene
    una matriz de n x m
    """
    R = []
    for i in range(0, n):
        R_aux = []
        for j in range(0, m):
            R_aux.append(A[i][j])
        R.append(R_aux)
    return R

def deter(A):
    """
    Funcion auxiliar para el calculo del determinantes de una
    matriz diagonal superior o inferior.
    Entrada: Una matriz diagonal superior o inferior
    Salida: Determinante de la matriz indicada
    """
    m = len(A)
    det = 1
    for i in range(m):
        det*=A[i][i]
    return det

def sub_prin(A, tol):
    """
    Función auxiliar para la validación del metodo de factorización LU
    por medio del análisis del determinante de las submatrices principales
    """
    m = len(A)

    for i in range (1, m):
        Ak = A
        Ak = resize(A, i, i)
        deter = np.linalg.det(Ak)

        if deter < tol:
            return False

    return True

### Funciones Principales ###

def fact_LU(A):
    """
    Metodo que realiza la respectiva factorización LU de una matriz
    Entarda: Una matriz a la cual sí se le puede realizar factorización LU
    Salida: Un arreglo de matrices que contiene a las reapectivas matrices L y U 
    de la factorización.
    """
    m = len(A)
    L = eye(m,m)

    for j in range(0, m):

        diag_el = A[j][j]
        for fila in range(j+1, m):
            multiplo = A[fila][j]/diag_el

            L[fila][j] = multiplo

            for column in range(0, m):
                A[fila][column] -= multiplo*A[j][column]
    U = A
    return [L, U]

def det_fact_lu(A):
    """
    Funcion que calcula el determinante de una matriz
    por medio de la factorización LU y las propiedades de 
    los determinantes.
    Entarda: Una matriz a la cual sí se le puede realizar factorización LU
    Salida: El determinante de la matriz
    """
    tol = 10**-9
    if sub_prin(A, tol):
        print("La metriz no tiene factorización LU")
        return "NA"
    
    l_u = fact_LU(A)

    L = l_u[0]
    U = l_u[1]

    deterL = deter(L)
    deterU = deter(U)
    deterA = deterL * deterU

    return deterA

### Ejemplo Numérico ###

A = [[2, 3, 0, 1, -1, 0, 1, 0],
     [4, 5, 1, 3, 0, 3, 4, 4], 
     [4, 3, 5, 4, 7, 13, 10, 13], 
     [6, 11, -4, 3, -9, -8, -5, -9], 
     [8, 9, 5, 4, 1, 10, 19, 17], 
     [2, 5, -6, -2, -16, -17, 2, -4], 
     [8, 7, 17, -3, 16, 32, 42, 33], 
     [4, 3, 11, -1, 16, 14, 18, 7]]

det = det_fact_lu(A)
print("El determinante de la matriz es: " + str(det))


