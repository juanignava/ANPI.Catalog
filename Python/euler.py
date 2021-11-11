import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

#### Para esta funcion se necesita el metodo de lagrange #########
def fun_Lk (xv, k):
    """
    Calculo de los factores de variable de Lagrange. Se recibe
    el vector de preimagenes y la iteracion actual en la cual
    se indefine el polinomio.
    Esta corresponde a una funcion auxiliar de la funcion lagrange()
    """
    x = sp.Symbol('x')
    n = len(xv)
    Lk = 1
    for j in range(n):
        if j != k:
            num = x - xv[j]
            den = xv[k]-xv[j]
            Lk *= num/den
    
    Lk_expand = sp.expand(Lk)

    return Lk_expand

def lagrange(xv, yv):
    """
    Esta funcion aproxima numericamente el polinomio de interpolacion
    que pasa por los puntos (x0, y0), (x1, y1), ..., (xn, yn) dados
    en las entradas por medio del metodo de Lagrange.

    Parametros de entrada:
        xv: corresponde al arreglo contituido por las preimagenes de
        los puntos considerados en el polinomio.

        yv: corresponde al arreglo constituido por las imagenes de los
        puntos considerados en el polinomio.

    Parametros de salida:
        pol: corresponde al polinomio de forma simbolica generado por
        medio de las diferencias divididas de Newton.
    """
    # Definicion de las constantes iniciales
    x = sp.Symbol('x')
    n = len(xv)
    p = 0

    # Formacion del polinomio
    for k in range(n):
        p += yv[k]*fun_Lk(xv, k)
    
    # Distribuir factores del polinomio
    pol = sp.expand(p)
    return pol

#### Fin del metodo de lagrange #########

def euler(f, intervalo, y0, num_puntos):
    
    """
    Esta funcion aproxima la solucion a la ecuacion 
    deferencial de la forma:

    y'(x) = f(x, y(x))
    y(x0) = y0

    por medio del metodo de euler.

    Parametros de entrada:
        f: corresponde a la función de f(x, y(x)) de la
            ecuacion diferencial anterior. Se ingresa en
            formato de string.
        intervalo: corresponde a un arreglo que contiene
            los valores inicial y final del intervalo
            en análisis.
        y0: es la imagen del valor conocido en la ecuacion
            diferencial.
        num_puntos: es el numero de puntos que se utilizara
            en el analisis del meotod de euler.
    """

    # Definicion variables iniciales
    x = sp.Symbol('x')
    y = sp.Symbol('y')
    f1 = sp.sympify(f)

    # calculo de h
    a = intervalo[0]
    b = intervalo[1]
    h = (b-a)/num_puntos

    # Calculo de las preimagenes
    xv = []
    for j in range(num_puntos+1):
        new_value = a + j*h
        xv.append(new_value)

    # calculo de las imagenes
    yv = [y0]
    for i in range(num_puntos):
        new_value = yv[i] + h * sp.N(f1.subs({x:xv[i], y:yv[i]}))
        yv.append(new_value)

    # polinomio de interpolacion
    pol = lagrange(xv, yv)

    # graficacion del polinomio de interpolacion
    plt.rcParams.update({'font.size': 14})
    fig, graf = plt.subplots()
    graf.plot(xv, yv, 'b', marker='o', markerfacecolor='red', markersize=10)
    graf.set_xlabel('Solucion  ($y(x)$)')
    graf.set_ylabel('$x$')
    graf.set_title('Solucion de la ecuacion por medio del metodo de Euler')
    graf.grid(True)
    plt.show()

    return [xv, yv, pol]

# Ejemplo numerico

f = 'y-x^2+1'
intervalo = [0, 5]
num_puntos = 9
y0 = 0.5

R = euler(f, intervalo, y0, num_puntos)

print("\nEl intervalo de preimagenes es:\n")
print(R[0])
print("\nEl intervalo de imagenes es:\n")
print(R[1])
print("\nEl polinomio de interpolacion es:\n")
print(R[2])
