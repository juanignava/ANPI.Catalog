# Se importan las librerias a utilizar
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

def variables_simbolicas(variables):

    """
    Esta funcion convierte las variables recibidas de tipo string a varibles de
    tipo simbolico.

    Parametros de entrada: variables = es una cadena de caracteres (string) que
                                       representa cada una de las variables a
                                       utilizar.

    Parametros de salida: v_var = corresponde a un vector de n dimensiones que
                                  que contine cada una de las variables de tipo
                                  simbolico.
    """

    n = len(variables)
    tamano = np.arange(0,n,2)
    v_var = []

    for i in tamano:
        v_var.append(sp.Symbol(variables[i]))

    return v_var

def gradiente(f, v_var):

    """
    Esta funcion calcula el gradiente de la funcion dada con base en las variables
    recibidas.

    Parametros de entrada: f = es una cadena de caracteres (string) que representa
                               la funcion que esta igualada a cero.

                           v_var = corresponde a un vector de n dimensiones que
                                   que contine cada una de las variables de tipo
                                   simbolico.

    Parametros de salida: g = hace referencia a un vector que representa el
                              gradiente de la funcion dada.
    """

    n = len(v_var)
    g = []

    for i in np.arange(0,n,1):

        g.append(sp.diff(f, v_var[i]))

    return g

def coordinado(f, variables, x0, tol, itermMax):

    """
    Esta funcion implementa el metodo del descenso coordinado,
    el cual permite aproximar la optimizacion de una funcion
    por medio de ciclos iterativos

    Parametros de entrada: f = es una cadena de caracteres (string) que representa
                               la funcion que esta igualada a cero.

                           variables = hace referencia a una cadena de caracteres
                                       (string) que contiene cada una de las
                                       variables de la funcion.

                           x0 = corresponde a un vector que contiene los valores
                                iniciales de cada una de las variables requeridas.

                           tol = corresponde a un valor numerico positivo que
                           representa la tolerancia del metodo.

                           iterMax = valor numerico positivo que representa el
                           numero maximo de iteraciones a realizar.

    Parametros de salida: xk = es el vector de n dimensiones donde se ubica el
                               punto optimizado por la funcion.

                          err =  corresponde a un vecotr que contiene la norma
                                 del gradiente de cada una de las iteraciones
                                 realizadas.
    """

    fs = sp.sympify(f)

    # v_var[0] = x | v_var[1] = y
    v_var = variables_simbolicas(variables)

    gs = gradiente(fs, v_var)
    gn = sp.lambdify(v_var,gs)

    x_vector = []
    x_vector.append(x0)

    error = 0
    err = []

    minimos = []

    for k in range(itermMax):

        xk = x_vector[k]

        # Evaluar Y para encontrar X
        expr = fs.subs(v_var[1], xk[1])
        fn = sp.lambdify(v_var[0], expr)
        minimos = optimize.fmin(fn, 1)
        xk[0] = minimos[0]

        # Evaluar X para encontrar Y
        expr = fs.subs(v_var[0], xk[0])
        fn = sp.lambdify(v_var[1], expr)
        minimos = optimize.fmin(fn, 1)
        xk[1] = minimos[0]

        x_vector.append(xk)
        error = np.linalg.norm(gn(xk[0], xk[1]))
        err.append(error)

        # Condicion de parada
        if error < tol:

            break

    # Instrucciones de graficacion
    plt.rcParams.update({'font.size': 14})
    fig, graf = plt.subplots()
    k = len(err)
    ejex = np.arange(0, k, 1)
    graf.plot(ejex, err, 'b', marker='o', markerfacecolor='red', markersize=10)
    graf.set_xlabel('Iteraciones ($k$)')
    graf.set_ylabel('$||∇f(xk)||$')
    graf.set_title('Descenso Coordinado (Iteraciones vs Error)')
    graf.grid(True)
    plt.show()

    return [xk, err]


# Ejemplo numerico
f = '(x-2)^2 + (y+3)^2 + x*y'
variables = 'x y'
x0 = [1,1]
tol = 10**-4
iterMax = 1000

y = coordinado(f,variables, x0, tol, iterMax)
print('Aproximacion = ', y[0])
print('Error = ', y[1])
