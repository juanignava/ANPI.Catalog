import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

xv = [2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4]
yv = [1, 1.191, 1.5964, 1.8883, 2.1757, 2.4420, 2.6964, 2.9414, 3.1791, 3.4112, 3.6388]

# graficacion del polinomio de interpolacion
plt.rcParams.update({'font.size': 14})
fig, graf = plt.subplots()
graf.plot(xv, yv, 'b', marker='o', markerfacecolor='red', markersize=10)
graf.set_xlabel('Iteraciones ($k$)')
graf.set_ylabel('$|f(x_K)|$')
graf.set_title('Metodo de la biseccion (Error vs Iteraciones)')
graf.grid(True)
plt.show()