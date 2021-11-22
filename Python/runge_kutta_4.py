import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

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

    Esta funcion se utilizara para calcular el polinomio de interpolacion
    una vez que se obtengan los valores de xv y yv mediante el metodo de 
    Runge-Kutta
    """
    # Definicion de las constantes iniciales
    x = sp.Symbol('x')
    n = len(xv)
    p = 0

    # Formacion del polinomio
    for k in range(n):
        p += yv[k]*fun_Lk(xv, k)
    
    # Distribuir factores del polinomio
    polinomio = sp.expand(p)
    return polinomio

def runge_kutta_4(f, intervalo, y0, N):

    """
       Esta  funcion  aproxima  la  solucion a la  ecuacion
       diferencial  de la forma:y’(x) = f(x, y(x))y(x0) = y0
       por  medio  del  metodo  de euler.
       
       Parametros  de  entrada:

            f: corresponde a la  f u n c i n  de f(x, y(x)) de la
               ecuacion  diferencial  anterior. Se  ingresa  en
               formato  de  string.
            
            intervalo: corresponde a un  arreglo  que  contiene
                       los  valores  inicial y final  del  intervalo
                       en analisis.
            
            X0: es un valor inicial conocido
            
            y0: es la  imagen  del  valor  conocido  en la  ecuacion
                diferencial.
            
            N: es el  numero  de  puntos  que se  utilizara
               en el  analisis  del  metodo  de  euler.
    
    """

    # Definicion variables iniciales

    x = sp.Symbol('x')
    y = sp.Symbol('y')
    fun = sp.sympify(f)

    # calculo de h
    a = intervalo[0]
    b = intervalo[1]
    h = (b-a)/(N-1)

    # Calculo de las preimagenes
    xk = []
    for j in range(N):
        xn = a + j*h
        xk.append(xn)

    # calculo de las imagenes
    yk = [y0]
    for i in range(N-1):
        k1 = sp.N(fun.subs({x:xk[i], y:yk[i]}))
        k2 = sp.N(fun.subs({x:(xk[i]+h/2), y:(yk[i]+h*(k1/2))}))
        k3 = sp.N(fun.subs({x:(xk[i]+h/2), y:(yk[i]+h*(k2/2))}))
        k4 = sp.N(fun.subs({x:(xk[i]+h), y:(yk[i]+h*k3)}))
        yn = yk[i] + (h/6) * (k1+2*k2+2*k3+k4)
        yk.append(yn)
    
    # Obtencion del polinomio de interpolacion
    polinomio = lagrange(xk, yk)

    # graficacion del polinomio de interpolacion
    plt.rcParams.update({'font.size': 14})
    fig, graf = plt.subplots()
    graf.plot(xk, yk, 'b', marker='o', markerfacecolor='red', markersize=10)
    graf.set_xlabel('Solucion  ($y(x)$)')
    graf.set_ylabel('$x$')
    graf.set_title('Solucion de la ecuacion utilizando el metodo de Runge-Kutta de orden 4')
    graf.grid(True)
    plt.show()

    return [xk, yk, polinomio]

# Ejemplo numerico

f = '-x*y+(4*x)/y'
intervalo = [0, 1]
N = 11
y0 = 1

R = runge_kutta_4(f, intervalo, y0, N)

print("El intervalo de preimagenes es:")
print(R[0], "\n")
print("El intervalo de preimagenes es:")
print(R[1], "\n")
print("El polinomio de interpolacion es:")
print(R[2], "\n")