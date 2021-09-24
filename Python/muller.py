import  sympy as sp
import  numpy as np
import  matplotlib.pyplot  as plt
import math

def muller(f, x0, x1, x2, tol, iterMax):

    """Esta  funcion  aproxima  numericamente  la  solucion  de una  
       ecuacion no  lineal  por  medio  del metodo  de la  Muller.
       
       Parametros  de  entrada:    f = una  cadena  de  caracteres (string) 
                                   que  representala funcion  que  esta  igualada 
                                   a cero.

                                   x0, x1, x2 =   son  los  valores iniciales que 
                                   utiliza el método de Muller.

                                   tol = corresponde a un valor  numerico  positivo  que
                                   representa  la  tolerancia  del  metodo.

                                   iterMax = valor  numerico  positivo  que  representa  el
                                   numero  maximo  de  iteraciones a realizar.
                                   
       Parametros  de  salida:     aprox = corresponde a la  aproximacion  del  cero 
                                   de la ecuacion.

                                   err = es |f(xk)|, error  absoluto  del  metodo  para xk."""


    x = sp.Symbol('x')
    fun = sp.sympify(f)
    aprox = 0
    err = 0
    er = []
    i = 0

    while (i < iterMax):

        f_x0 = sp.N(fun.subs(x, x0))
        f_x1 = sp.N(fun.subs(x, x1))
        f_x2 = sp.N(fun.subs(x, x2))

        d1 = f_x0 - f_x2
        d2 = f_x1 - f_x2
        h1 = x0 - x2
        h2 = x1 - x2

        a0 = f_x2 # a

        a1 = (((d2 * pow(h1, 2)) -
            (d1 * pow(h2, 2))) /
            ((h1 * h2) * (h1 - h2))) # b
        
        a2 = (((d1 * h2) - (d2 * h1)) /
            ((h1 * h2) * (h1 - h2))) # c
        
        s = ((-2 * a0) / (a1 +
            abs(math.sqrt(a1 * a1 - 4 * a0 * a2)))) # 2c/(b+sgn(b)*(b**2-4ac)**0.5)
        
        p = ((-2 * a0) / (a1 -
            abs(math.sqrt(a1 * a1 - 4 * a0 * a2)))) # 2c/(b+sgn(b)*(b**2-4ac)**0.5)

        if (s >= p):
            aprox = s + x2
        else:
            aprox = p + x2

        x0 = x1
        x1 = x2
        x2 = aprox

        err = abs(fun.subs(x, x2))
        er.append(sp.N(err))
        if err < tol: 
            break
        i += 1

    # Instrucciones de graficacion

    plt.rcParams.update({'font.size': 14})
    fig, graf = plt.subplots()
    ejex=np.arange(1, i+2, 1)
    graf.plot(ejex, er, 'b', marker='o', markerfacecolor='red', markersize=10)
    graf.set_xlabel('Iteraciones ($k$)')
    graf.set_ylabel('$|f(x_2)|$')
    graf.set_title('Metodo de Muller (Error vs Iteraciones)')
    graf.grid(True)
    plt.show()

    return x2


f = 'x**3 + 2*x**2 +10*x-20'
x0 = 0
x1 = 1
x2 = 6
tol = 10**-9
iterMax = 1000
result = muller(f, x0, x1, x2,  tol, iterMax)

print(result)