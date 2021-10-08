#from gaussiana import sustitucion_adelante
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import math

def matriz_ceros(n, m):

    """
    Esta funcion verifica crea una matriz de ceros con base en las
    dimensiones recibidas.

    Parametros de entrada : n corresponde a la cantidad de filas de la
                            matriz a crear.

                            m hace referencia a la cantidad de columnas
                            de la matriz a crear.

    Parametros de salida : matriz de ceros creada.
    """

    resultado = []

    for i in range(0, n):

        aux = []

        for j in range(0, m):

            aux.append(0)

        resultado.append(aux)

    return resultado

def matriz_cuadrada(A):
	
	"""
    Esta funcion verifica que una matriz dada sea una matriz cuadrada.

    Parametros de entrada : A corresponde a la matriz mxn.

    Parametros de salida : valor booleano que indica si la matriz es 
							cuadrada.
    """
	m = len(A)

	if (m > 0):
		n = len(A[0])
		if(m == n):
			return True
		else:
			return False
	else:
		print("La matriz ingresada está vacia")

def diagonal_dominante(A):

    """
    Esta funcion verifica que una matriz dada sea diagonal dominante.

    Parametros de entrada : A corresponde a la matriz mxm.

    Parametros de salida : valor booleano que indica si la matriz es 
						   diagonal dominante.
    """

    if(matriz_cuadrada(A)):
        m = len(A)
        for i in range(0, m):
            dom = 0
            aux = 0
            for j in range(0,m):
                if(i == j):
                    dom = abs(A[i][j])
                else:
                    aux = aux + abs(A[i][j])
            if (dom <= aux):
                return False
        return True
    else:
        print("La matriz no es cuadrada")
        return False

def obtener_L_D_U(A):

    """
    Esta funcion descompone la matriz A de tal forma que A = L+D+U.

    Parametros de entrada : A corresponde a la matriz mxm.

    Parametros de salida : L corresponde a la matriz triangular 
                           inferior de A
                           
                           D corresponde a la matriz con la diagonal 
                           principal de A

                           U corresponde a la matriz triangular 
                           superior de A
    """
    
    m = len(A)
    L = matriz_ceros(m,m)
    D = matriz_ceros(m,m)
    U = matriz_ceros(m,m)

    for i in range(0,m):
        for j in range(0,m):
            if (i == j):
                D[i][j] = A[i][j]
            elif(i < j):
                U[i][j] = A[i][j]
            else:
                L[i][j] = A[i][j]
    
    return L, D, U

def sust_adelante(A, b):

    m = len(A)
    x = zeros(m, 1)

    for i in range(0,m):
        aux = 0

        for j in range (0, i):
            aux += A[i][j]*x[j]

        x[i] = (1/A[i][i])*(b[i]-aux)

    return x

def sustitucion_adelante(A, b):

    """
    Esta funcion se encarga de aplicar el metodo de sustitucion
    hacia adelante para resolver un sistema Ax = b.

    Parametros de entrada : A es una matriz de cualquier dimension
                            la cual es la matriz de coeficientes del
                            sistema Ax = b.

                            b es una matriz de Nx1 que corresponde a
                            la matriz de terminos independientes.

    Parametros de salida : matriz x que hace referencia a la solucion
                           del sistema Ax = b.
    """
    
    n = len(b)

    x = matriz_ceros(n, 1)

    for i in range(0, n):

        aux = 0

        for j in range (0, i):

            aux += A[i, j] * x[j]

        x[i] = (b[i] - aux) / A[i, i]

    return x

def gauss_seidel(A, b, x0, tol, iterMax):

    """
    Esta funcion aproxima la solucion de un sistema Ax = b por medio
    del metodo de Gauss_Seidel.

    Parametros de entrada : A es una matriz de cualquier dimension
                            la cual es la matriz de coeficientes del
                            sistema Ax = b.

                            b es una matriz de Nx1 que corresponde a
                            la matriz de terminos independientes.

                            x0 corresponde a una matriz inicial

                            tol es la tolerancia del método

                            iterMax es la cantidad de iteraciones 
                            maximas permitidas por la funcion

    Parametros de salida : xk es la aproximacion a la solución del sistema
                              Ax = b.
                           err es el error obtenido en la ultima iteracion
    """

    if(len(A) != len(b)):
        return "Tamano de b no corresponde con la dimensión de la matriz A"
    if(len(A) != len(x0)):
        return "Tamano de x0 no corresponde con la dimensión de la matriz A"

    if(diagonal_dominante(A)):

        A_matrix = np.array(A)
        b_vector = np.transpose(np.matrix(b))

        r = obtener_L_D_U(A)
        L = np.array(r[0])
        D = np.array(r[1])
        U = np.array(r[2])

        LpD = A_matrix - U

        print("L + D")
        print(LpD)
        y = sustitucion_adelante(LpD, b)
        print(y)
        print("valor de B")
        print(b)

        xk = x0
        err = 0
        er = []
        k = 0

        while k < iterMax:

            k = k+1
            print("valor de U")
            print(U)
            print("valor de xk")
            print(xk)
            print("producto U*xk")
            print(np.matmul(U, xk))
            z = sustitucion_adelante(LpD, np.matmul(U, xk))
            print("valor de z")
            print(z)
            xk = -1*np.array(z) + np.array(y)
            print("valor de xk")
            print(xk)


            err = np.linalg.norm(np.matmul(A_matrix, xk)- np.array(b))
            print("error")
            print(err)
            er.append(sp.N(err))

            if(err < tol):
                break
        
        plt.rcParams.update({'font.size': 14})
        fig, graf = plt.subplots()
        ejex=np.arange(1, k+1, 1)
        graf.plot(ejex, er, 'b', marker='o', markerfacecolor='red', markersize=10)
        graf.set_xlabel('Iteraciones ($k$)')
        graf.set_ylabel('$||Ax-b||$')
        graf.set_title('Metodo de Jacobi (Error vs Iteraciones)')
        graf.grid(True)
        plt.show()

        return [xk, err]

    else:
        return "La matriz ingresada no es diagonal dominante"

A = [[10, 1, -1, 5],
     [0, -7, 1, 1],
     [1, 2, 5, 0],
     [4, 4, 1, 10]]

C = [[10, -3, 1],
     [2, -10, 3],
     [0, -1, 2]]

b = [1, 1, 1, 1]

tol = 10**-5
iterMax = 10

x0 = [1,1,1,1]

R = gauss_seidel(A, b, x0, tol, iterMax)

print("aproximacion: ")
print(R[0])
print("error: ")
print(R[1])
