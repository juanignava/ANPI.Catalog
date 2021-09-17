"""
Implementación del método de Newton - Raphson
"""

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

def newton_raphson(f, x0, tol, iterMax):
    x = sp.Symbol('x')
    f1 = sp.sympify(f)
    df1 = sp.diff(f1,x)
    er = []
    err = tol+1 # Como voy a utilizar un while necesito que se cumpla esa condición
    k = 0
    xk = x0

    while err>tol and k<iterMax:
        k = k+1
        n = sp.N(f1.subs(x, xk)) # esto convierte de simbolico a numerico
        d = sp.N(df1.subs(x, xk))
        """
        Verificar que d sea diferente de 0 o que en valor 
        absoluto sea menor a una tolerancia dada (10**-15)
        """
        xk = xk - n/d
        err = abs(f1.subs(x,xk))
        er.append(sp.N(err))

    plt.rcParams.update({'font.size': 14})
    fig, graf = plt.subplots()
    ejex=np.arange(1, k+1, 1) #Crea un vector que va de 1 a k+1 con saltos de 1
    graf.plot(ejex, er, 'b', marker='o', markerfacecolor='red', markersize=10)
    graf.set_xlabel('Iteraciones ($k$)')
    graf.set_ylabel('$|f(x_K)|$')
    graf.set_title('Metodo de Newton Raphson (Interacciones vs Error)')
    graf.grid(True)
    plt.show() ## Importante usar el .show para que se muestre la ventana emergente

    return [xk, k, err]

f = "exp(x)+x-2" # se debe escribir como texto
x0 = 5
tol = 10**-9
iterMax = 100

y = newton_raphson(f, x0, tol, iterMax)
print(y)