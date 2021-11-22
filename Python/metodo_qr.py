import numpy as np
from numpy.lib.twodim_base import diag
import sympy as sp
import matplotlib.pyplot as plt


def metodo_qr(A0, U0, tol, iterMax):

    Ak = np.matrix(A0)
    Uk = np.matrix(U0)

    print("La matriz con la que se empieza es:\n")
    print(Ak)

    error = 0
    er = []
    k = 0

    while k < iterMax:

        #Valores propios de Ak
        Vpk = np.diag(Ak)

        k = k+1

        Qk, Rk = np.linalg.qr(Ak)
        print("el valor de Q")
        print(Qk)
        print("el valor de R")
        print(Rk)
        Ak = np.matmul(Rk, Qk)
        print(Ak)
        Uk = np.matmul(Uk, Qk)

        #Valores propios de Ak+1
        Vpk1 = np.diag(Ak)

        error = np.linalg.norm(Vpk1-Vpk)/np.linalg.norm(Vpk1)
        er.append(sp.N(error))


        if error < tol:
            break
    
    plt.rcParams.update({'font.size': 14})
    fig, graf = plt.subplots()
    ejex=np.arange(1, k+1, 1)
    graf.plot(ejex, er, 'b', marker='o', markerfacecolor='red', markersize=10)
    graf.set_xlabel('Iteraciones ($k$)')
    graf.set_ylabel('Error')
    graf.set_title('Metodo QR (Error vs Iteraciones)')
    graf.grid(True)
    plt.show()
    
    valoresPropios = np.diag(Ak)
    vectoresPropios = Uk

    return [valoresPropios, vectoresPropios]

# Ejemplo numerico

A0 = [[0, 11, -5],
      [-2, 17, -7],
	  [-4, 26, 10]]

U0 = [[1, 0, 0],
      [0, 1, 0],
	  [0, 0, 1]]

tol = 10**-8
iterMax = 100

r = metodo_qr(A0, U0, tol, iterMax)

print("\n Valores propios y vectores propios obtenidos mediante el Metodo QR \n")

for i in range(0, len(A0)):
    print("Valor propio", i+1, "=", r[0][i])

print("\n Vectores propios: \n")
print(r[1])





