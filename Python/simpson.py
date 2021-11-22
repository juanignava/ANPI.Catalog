import  sympy as sp

def simpson(fun, a, b):

    """
          Esta  funcion  aproxima  implementa la Regla de Simpson la cual permite
          aproximar el resultado de una integral definida en un intervalo dado por 
          medio de la tercera formula de Newton-Cotes.

          Parametros  de  entrada:    fun = es una cadena de caracteres que 
                                            representa la funcion que será evaluada 
                                            en la integral definida.

                                      a = representa el limite inferior de 
                                          la integral.

                                      b = representa el limite superior de 
                                          la integral.

          Parametros  de  salida:     aprox = resultado obtenido al evaluar la 
                                              integral en el intervalo dado.

                                      error = representa la cota de error de la 
                                              Regla de Simpson.
    
    """

    x = sp.Symbol('x')
    f = sp.sympify(fun)
    d1 = sp.diff(f,x)
    d2 = sp.diff(d1,x)
    d3 = sp.diff(d2,x)
    d4 = sp.diff(d3,x)

    # Obtencion de variables necesarias para obtener el resultado de la aproximacion

    h = (b-a)/2
    x0 = a
    x1 = (a+b)/2
    x2 = b

    fx0 = sp.N(f.subs(x, x0))
    fx1 = sp.N(f.subs(x, x1))
    fx2 = sp.N(f.subs(x, x2))

    aprox = (h/3)*(fx0 + 4*fx1 + fx2)

    # Verificacion del maximo del intervalo

    d4a = sp.N(d4.subs(x, a))
    d4b = sp.N(d4.subs(x, b))

    if abs(d4a) > abs(d4b):
        e = a

    else:
        e = b
    
    d4e = sp.N(d4.subs(x, e))

    error = (h**5/90)*abs(d4e)

    return aprox, error

# Declaracion de las variables a utilizar
fun = 'ln(x)'
a = 2
b = 5

r = simpson(fun, 2, 5)

print("Aproximacion: ", r[0])
print("Cota de error: ", r[1])

