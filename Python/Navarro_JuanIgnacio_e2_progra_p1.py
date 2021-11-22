"""
EXAMEN PARCIAL 2 - ANPI
Pregunta progrmación 1
Estudiante: Juan Ignacio Navarro Navarro (2019039662)
"""
import sympy as sp
import matplotlib.pyplot as plt


def runge_kutta_4(a, b, y0, m):

    """
    Metodo de runge kutta de orden 4 para encontrar la solucion al probleme
    y' = (x+y)/x con y(2) = 4
    """

    # definicion de la funcion
    f = '(x+y)/x'
    x = sp.Symbol('x')
    y = sp.Symbol('y')
    fun = sp.sympify(f)

    # calculo de h
    h = (b-a)/ (m-1)

    # calculo de las preimagenes
    xk = []
    for j in range(m):
        xn = a + j*h
        xk.append(xn)

    # calculo de las imagenes
    yk = [y0]
    for i in range(m-1):

        k1 = sp.N(fun.subs({x:xk[i], y:yk[i]}))

        k2 = sp.N(fun.subs({x:(xk[i]+h/2), y:(yk[i]+h*(k1/2))}))

        k3 = sp.N(fun.subs({x:(xk[i]+h/2), y:(yk[i]+h*(k2/2))}))

        k4 = sp.N(fun.subs({x:(xk[i]+h), y:(yk[i]+h*k3)}))

        yn = yk[i] + (h/6) * (k1 + 2*k2 + 2*k3 + k4)

        yk.append(yn)

    # graficacion de los puntos que aproximan la solucion
    plt.rcParams.update({'font.size': 14})
    fig, graf = plt.subplots()
    graf.plot(xk, yk, 'b', marker='o', markerfacecolor='red', markersize=10)
    graf.set_xlabel('Preimagenes  ($xv$)')
    graf.set_ylabel('Imagenes ($yv$)')
    graf.set_title('Puntos que aproximan la solucion del problema en el metodo de Runge-Kutta de orden 4')
    graf.grid(True)
    plt.show()

    return [xk, yk]


# Caso numerico del enunciado
resultado = runge_kutta_4(2, 10, 4, 50)

print("Los el vector de las preimagenes es: \n")
print(resultado[0])
print("Los el vector de las imagenes es: \n")
print(resultado[1])

