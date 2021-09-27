# Primer ejericico del examen parcial
# Análisis Numérico para Ingeniería
# Estudiante: Juan Ignacio Navarro Navarro (2019039662)

import sympy as sp

def metodo_biseccion_mod(a, b, n, tol):
    """
    Esta función permite calcular el cero de la función
    f(x)= 2+cos(exp(x)-2)-exp(x)utilizando el método modificado
    de la bisección

    Parametros de entrada:
            a, b: valores del intervalo en análisis
            tol: tolerancia del criterio de parada
            n: cantidad de subintervalos en los que 
            se dividirá el intervalo
    
    Parámetros de salida:
        xk: aproximación numérica del cero de la función especificada.
    """

    x = sp.Symbol('x')
    f_s = '2 + cos(exp(x)-2)- exp(x)' # Función en análisis
    f = sp.sympify(f_s)
    # contantes iniciales
    error = tol+1
    c = []
    xk = 0

    # Verificación del teorema de Bolzano
    if sp.N(f.subs(x, a)) * sp.N(f.subs(x, b)) >= 0: 
        print("El intervalo inicial no cumple con el teorema de Bolzano")
        return "NA"
    
    while True:
        # Definición del arreglo de particiones c
        h = (b-a)/n
        for j in range(n+1):
            c_j = a + j*h
            c.append(c_j)

        # Busqueda de los términos del arreglo c que contienen
        # el cero de la funcion
        for j in range(n):
            # Condició de Bolzano
            if sp.N(f.subs(x, c[j])) * sp.N(f.subs(x, c[j+1])) < 0:
                # Definición de los valores de la próxima iteración
                a = c[j]
                b = c[j+1]
                xk = (c[j+1] + c[j])/2

        # Calculo del error de la ecuación
        error = abs(f.subs(x, xk))
        if error < tol:
            break
        c = []

    return xk # Retornamos el valor aproximado

# Ejemplo numérico
y = metodo_biseccion_mod(0, 2, 10, 10**-10)
print(y)


