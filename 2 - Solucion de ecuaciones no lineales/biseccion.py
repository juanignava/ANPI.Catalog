import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

def biseccion(f, a, b, tol, iterMax):
    """
    Esta funcion aproxima numericamente la solucion de una ecuacion 
     no lineal por medio del metodo de la biseccion para un intervalo dado.
    
    Parametros de entrada: f = una cadena de caracteres (string) que representa
                           la funcion que esta igualada a cero.
                           a, b =  son los extremos del intervalo [a, b] usado 
                           en el metodo de la biseccion.
                           tol = corresponde a un valor numerico positivo que 
                           representa la tolerancia del metodo.
                           iterMax = valor numerico positivo que representa el
                           numero maximo de iteraciones a realizar.

    Parametros de salida: xk = correpsonde a la aproximacion del cero de la
                           ecuacion.
                           err = es |f(xk)|, error absoluto del metodo para xk.
    """
    
    x = sp.Symbol('x')
    f1 = sp.sympify(f)
    xk = a
    err = 0
    er = []
    k = 0

    if sp.N(f1.subs(x, a)) * sp.N(f1.subs(x, b)) < 0: # Teorema de Bolzano

        while k < iterMax: # Iteraciones en el metodo
            k = k+1
            xk = (a+b)/2

            if sp.N(f1.subs(x, a)) * sp.N(f1.subs(x, xk)) < 0:
                b = xk
            else:
                a = xk

            err = abs(f1.subs(x, xk))
            er.append(sp.N(err))

            if err < tol: # Parar si el error es menor a tol
                break
    else:
        xk = 'NA'
        err = 'NA'
        print ("No se cumple el teorema de Bolzano")

    # Instrucciones de graficacion

    plt.rcParams.update({'font.size': 14})
    fig, graf = plt.subplots()
    ejex=np.arange(1, k+1, 1)
    graf.plot(ejex, er, 'b', marker='o', markerfacecolor='red', markersize=10)
    graf.set_xlabel('Iteraciones ($k$)')
    graf.set_ylabel('$|f(x_K)|$')
    graf.set_title('Metodo de la biseccion (Error vs Iteraciones)')
    graf.grid(True)
    plt.show()

    return [xk, err]

# Ejemplo numerico

f = 'exp(x)-6*x-10'
a = 2
b = 6
tol = 10**-5
iterMax = 1000

y = biseccion(f, a, b, tol, iterMax)
print(y)