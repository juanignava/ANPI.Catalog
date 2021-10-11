import sympy as sp
import numpy as np
from scipy import optimize

def cota_traz_cubico(f, xv):
    """
    Esta funcion realiza el calculo de la cota de error de un
    trazador cubico para una funcion f con respecto a un
    conjunto soporte xv.

    Parametros de entrada:
        f: corresponde a la funcion de donde se tomaron los
        puntos del calculo de interpolacion.

        xv: corresponde al conjunto sooprte sobre el cual se 
        realiza la interpolacion.

    Parametros de salida:
        cota: el valor numerico de la cota de error calculada.
    """

    # Calculo de h
    n = len(xv)
    dist = []
    for i in range(n-1):
        distance = xv[i+1] - xv[i]
        dist.append(distance)

    # deifinir h como el valor maximo de las distancias
    h = max(dist)

    x = sp.Symbol('x')
    f1  = sp.sympify(f)
    # calculo de la cuarta derivada de la funcion
    df4 = sp.diff(f1, x, 4)
    a = xv[0]
    b = xv[len(xv)-1]

    # calculo del valor maximo de la funcion en 
    # el intervalo [a,b]
    f_aux = -1*abs(df4)
    f_aux_num = sp.lambdify(x, f_aux)
    x_max = optimize.fminbound(f_aux_num, a, b)

    # calculo de la cota
    cota = (5*h**4/384)*sp.N(df4.subs(x, x_max))

    return cota

# Ejemplo numérico

f = 'exp(x/2)'
xv = [1, 1.5, 1.75, 2.15, 2.4, 3]
cota = cota_traz_cubico(f, xv)
print(cota)