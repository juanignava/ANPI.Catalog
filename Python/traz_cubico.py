import sympy as sp

def tridiagonal(h):
    """
    Esta funcion permite crear la matriz tridiagonal
    necesaria para el calculo del vector M en la
    formacion de los trazadores cubicos.

    Esta funcion es auxiliar de la llamada traz_cubico
    """
    # Definicion de las constantes iniciales
    coef = []
    n = len(h)

    # Calculo de cada uno de los elementos de la matriz
    for i in range(n-1):
        fila = []
        for j in range(n-1):
            elem = 0
            if i == j:
                elem = 2*(h[j] + h[j+1])
            elif j == i+1:
                elem = h[j]
            elif i == j+1:
                elem = h[i]
            fila.append(elem)
        coef.append(fila)

    return coef

def generar_u(yv, h):
    """
    Esta funcion permite crear el vector U necesitado
    para el calculo del vector M para la formacion
    de los trazadores cubicos.
    Esta funcion es auxiliar de la llamada traz_cubico
    """
    # Definicion de constantes iniciales
    u = []
    n = len(h)

    # Calculo de cada elemento de U
    for i in range(n-1):
        rest1 = (yv[i+2] - yv[i+1]) / h[i+1]
        rest2 = (yv[i+1]-yv[i]) / h[i]
        elem = 6*(rest1-rest2)
        u.append(elem)

    return u

def zeros(n):
    """
    Esta funcion corresponde a una auxiliar
    del metodo de Thomas para generar un vector
    0 de n elementos.
    """
    R = []
    for i in range(0, n):
        R.append(0)

    return R

def thomas(A, d):
    """
    Esta funcion calcula la solucion exacta de
    el sistema de ecuaciones Ax = b para una matriz
    de coeficientes A tridiagonal por medio del
    metodo de Thomas.

    Parametros de entrads:
        A: Matriz de coeficientes tridiagonal
        d: Vector de resultados del sistema.

    Parametro de salida:
        x: Vector solucion del sistema.
    """
    # Definicion de las constantes iniciales
    n = len(d)
    p = zeros(n-1)
    q = zeros(n)
    a = zeros(n)
    b = zeros(n)
    c = zeros(n)
    x= zeros(n)

    # Definicion de los vectores a, b y c
    for i in range(n):
        for j in range(n):
            if i+1==j:
                c[i] = A[i][j]
            if i==j+1:
                temp = A[i][j]
                a[i] = temp
            if i==j:
                b[i] = A[i][j]

    # Definicion de los vectores p y q
    for i in range(n):
        if i == 0:
            p[i] = c[i]/b[i]
            q[i] = d[i]/b[i]
        else:
            if i!=n-1:
                p[i] = c[i]/(b[i]-p[i-1]*a[i])
            q[i] = (d[i]-q[i-1]*a[i]) / (b[i]-p[i-1]*a[i])

    # Calculo de los elementos de x
    for i in range(n-1, -1, -1):
        if i == n-1:
            x[i] = q[i]
        else:
            x[i] = q[i]-p[i]*x[i+1]

    return x


def calc_const(M, h, yv, i):
    """
    Esta funcion es utilizada para calcular las constantes
    utilizadas como coeficientes de cada uno de los trazadores
    especificamente en el trazador i.

    Esta funcion es auxiliar de la llamada traz_cubico
    """
    a = (M[i+1] - M[i]) / (6 * h[i])
    b = M[i] / 2
    c = (yv[i+1] - yv[i]) / h[i] - (h[i] * (M[i+1] + 2 * M[i])) / 6
    d = yv[i]

    const = [a, b, c, d]

    return const

def traz_cubico(xv, yv):
    """
    Esta funcion calcula los trazadores cubicos relacionados
    a los puntos (x0, y0), (x1, y1), ..., (xn, yn) dados
    en las entradas.

    Parametros de entrada:
        xv: corresponde al arreglo contituido por las preimagenes de
        los puntos considerados en el polinomio.

        yv: corresponde al arreglo constituido por las imagenes de los
        puntos considerados en el polinomio.

    Parametros de salida:
        S: Corresponde a un arreglo que contiene un trazador cubico
        por cada una de las secciones de la funcion limitada por el
        conjunto soporte indicado.
    """

    # Definicion de constantes iniciales
    x = sp.Symbol('x')
    n = len(xv)

    # Calculo del arreglo de distancias h
    h = []
    for i in range(n-1):
        distance = xv[i+1] - xv[i]
        h.append(distance)

    # Calculo del arreglo M, este se necesita
    # para el posterior calculo de las constantes
    # a_i, b_i, c_i y d_i
    coef = tridiagonal(h)   # Matriz tridiagonal del metodo
    U = generar_u(yv, h)    # Vector U del metodo
    M = thomas(coef, U)     # M_1 a M_(n-1) es la solucion
                            #   de este sistema
    M = [0] + M + [0]       # M_0 y M_n son 0

    # Calculo de cada uno de los trazadores
    S = []
    for k in range(n-1):

        # Calculo de las constantes necesarias para
        #   la ecuacion k
        constantes = calc_const(M, h, yv, k)
        sum1 = constantes[0] * (x-xv[k])**3
        sum2 = constantes[1] * (x-xv[k])**2
        sum3 = constantes[2] * (x-xv[k])
        sum4 = constantes[3]
        
        S_k = sum1+sum2+sum3+sum4
        S_k = sp.expand(S_k) # distribuir el resultado
        S.append(S_k)

    return S


# Ejemplo numerico
# Con los puntos
#   (1, 2.718282)
#   (1.05, 3.286299)
#   (1.07, 3.527609)
#   (1.1, 3.905416)


xv = [1, 1.05, 1.07, 1.1]
yv = [2.718282, 3.286299, 3.527609, 3.905416]

trazadores = traz_cubico(xv, yv)

print("\nLos trazadores son:\n")

# Mostrat los resultados:
print("\nLos trazadores respectivos son los siguientes:")
k = 0
for trazador in trazadores:
    print("\nTrazador numero " + str(k) + ": ")
    print(trazador)
    k+=1

