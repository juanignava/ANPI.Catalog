import sympy as sp


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

    print("\nEl polinomio Lk de la iteracion " + str(k))
    print(Lk_expand)

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


# Ejemplo numerico
# usando los puntos:
#   (-2, 0)
#   (0, 1)
#   (1, -1)
xv = [-2, 0, 1]
yv = [0, 1, -1]
polinimio = lagrange(xv, yv)
print(polinimio)
