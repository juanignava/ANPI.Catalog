import sympy as sp

def dif_div(i, j, xv, yv):
    """
    Calculo del coeficiente de los sumandos del polinomio del metodo
    de diferencias divididas de Newton.
    Esta corresponde a una funcion auxiliar de la funcion dd_newton()
    """
    if i == j:
        return yv[i]

    else:
        n = dif_div(i+1, j, xv, yv) - dif_div(i, j-1, xv, yv)
        d = xv[j]-xv[i]
        return n/d


def dd_newton(xv, yv):
    """
    Esta funcion aproxima numericamente el polinomio de interpolacion
    que pasa por los puntos (x0, y0), (x1, y1), ..., (xn, yn) dados
    en las entradas por medio del metodo de diferencias divididas de 
    Newton.

    Parametros de entrada:
        xv: corresponde al arreglo contituido por las preimagenes de
        los puntos considerados en el polinomio.

        yv: corresponde al arreglo constituido por las imagenes de los
        puntos considerados en el polinomio.

    Parametros de salida:
        pol: corresponde al polinomio de forma simbolica generado por
        medio de las diferencias divididas de Newton.
    """

    # Definicion de constantes iniciales
    x = sp.Symbol('x')
    n = len(xv)
    pol = 0

    # Iteraciones que forman el polinomio
    for j in range(n):

        # Calculo de la diferencia dividida f[x0, ..., xj]
        val = dif_div(0, j, xv, yv)

        # Formacion la parte variable del polinomio
        # (x - x_0)*(x - x_1)*...*(x-x_(n-1))
        pol_var = 1
        for i in range (j):
            pol_var *= (x-xv[i])

        pol = pol + val * pol_var

    # Distribuir factores del polinomio
    pol_expand = sp.expand(pol)
    return pol_expand


# Ejemplo numérico
# usando los puntos:
#   (-2, 0)
#   (0, 1)
#   (1, -1)

xv = [-2, 0, 1]
yv = [0, 1, -1]
polinimio = dd_newton(xv, yv)
print(polinimio)