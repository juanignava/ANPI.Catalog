import sympy as sp

def simpson_compuesto(f, m, a, b):

    """
    Esta funcion implementa la Regla Compuesta de Simpson, la cual permite aproximar
    el resultado de una integral definida en un intervalo dado.
    
    Parametros de entrada: f es una cadena caracteres (string) que simboliza la
                           ecuacion que sera evaluada en la integral definida.

                           m hace referencia a la cantidad de puntos necesarios para
                           el calculo de la aproximacion. Debe ser un numero impar
                           
                           a corresponde al valor inicial del intervalo de la
                           integral.
                           
                           b representa el valor inicial del intervalo de la
                           integral.

    Parametros de salida: aprox es el resultado de aproximacion de la integral
                          evaluada en el intervalo dado.

                          error representa la cota de error de la Regla Compuesta
                          de Simpson.
    """

    # Verificacion de m impar
    if(m % 2 == 0):

        return "m debe ser impar"

    x = sp.Symbol('x')

    f1 = sp.sympify(f)

    d1 = sp.diff(f1, x)

    d2 = sp.diff(d1, x)

    d3 = sp.diff(d2, x)

    d4 = sp.diff(d3, x)

    h = (b - a) / (m - 1)

    print("valor de h")
    print(h)

    x0 = a

    xn = b

    xk = []

    # Calculo de los m puntos
    for i in range(1, m - 1):

        xi = a + i * h
        print("valor de preimagen: " + str(i) + " es: " + str(xi))
        xk.append(xi)
    
    sumaPares = 0

    cont = 1

    # Suma de terminos pares
    for i in range(1, len(xk), 2):
        
        sumaPares += sp.N(f1.subs(x, xk[i]))

    sumaImpares = 0
    
    # Suma de terminos impares
    for i in range(0, len(xk), 2):
        
        sumaImpares += sp.N(f1.subs(x, xk[i]))

    print("suma de los pares: " + str(sumaPares))
    print("suma de los impares: " + str(sumaImpares))

    aprox = (h / 3) * (sp.N(f1.subs(x, x0)) + 2 * sumaPares + 4 * sumaImpares +
            sp.N(f1.subs(x, xn)))

    # Verificacion del maximo del intervalo
    if (abs(sp.N(d4.subs(x, a))) > abs(sp.N(d4.subs(x, b)))):

        e = a;     
        
    else:
        
        e = b;

    error = (((b - a) * (h**4)) / 180) * (abs(sp.N(d4.subs(x, e))))

    return [aprox, error]

# Ejemplo numerico

f = 'exp(x)'
m = 7
a = 0
b = 1

y = simpson_compuesto(f, m, a, b)

print("Aproximacion = ", y[0])
print("Cota de Error = ", y[1])
