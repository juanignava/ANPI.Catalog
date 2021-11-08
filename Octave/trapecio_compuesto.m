function trapecio_compuesto()
  
  % Ejemplo numérico del mètodo de trapecio compuesto
  
  clc; clear;
  f = 'ln(x)';
  numPuntos = 4;
  intervalo = [2 5];
  
  aprox = metodo_trapecio_compuesto(f, numPuntos, intervalo)
  error = error_trapecio_compuesto(f, numPuntos, intervalo)
  
end

function aprox = metodo_trapecio_compuesto(f, numPuntos, intervalo)

  % Esta funcion implementa la aproximacion de la integral de la
  % función f en el intervalo dado por medio del metodo de el
  % trapecio compuesto para un numero de puntos dado.
  
  % Parametros de entrada:
  %   f ->  corresponde a un conjunto de caracteres que forman la
  %         funcion que se va a integrar.
  %   numPuntos -> corresponde a un entero que representa el numero
  %         de puntos que se utilizan en la aproximacion.
  %   intervalo -> corresponde a un arreglo con dos valores que 
  %         son el intervalo de integracion.
  %
  % Parametros de salida:
  %   aprox -> este es el valor de la aproximacion de la integral
  %         por mecio del metodo del trapecio compuesto.
  
  % Se inicializa el paquete symbolic
  pkg load symbolic;
  syms x;
  f1 = matlabFunction(sym(f));
  
  % Se cargan la variables iniciales
  h = ( intervalo(2) - intervalo(1) ) / (numPuntos-1);
  sum = 0;
  
  % Iteraciones que calculan cada una de las alturas de los trapecios
  for n=0 : numPuntos-2
  
    x_n = intervalo(1) + h*n;  
    sum = sum + f1(x_n) + f1(x_n+h);
    
  endfor
  
  % Aproximacion final considerando la base de cada trapecio
  aprox = (h/2)*sum;
    
end

function error = error_trapecio_compuesto(f, numPuntos, intervalo)
  
  % Esta funcion implementa la aproximacion de la integral de la
  % función f en el intervalo dado por medio del metodo de el
  % trapecio compuesto para un numero de puntos dado.
  
  % Parametros de entrada:
  %   f ->  corresponde a un conjunto de caracteres que forman la
  %         funcion que se va a integrar.
  %   numPuntos -> corresponde a un entero que representa el numero
  %         de puntos que se utilizan en la aproximacion.
  %   intervalo -> corresponde a un arreglo con dos valores que 
  %         son el intervalo de integracion.
  %
  % Parametros de salida:
  %   aprox -> este es el valor de la aproximacion de la integral
  %         por mecio del metodo del trapecio compuesto.
  
  % Se inicializa el paquete symbolic
  pkg load symbolic;
  syms x;
  fs = sym(f);
  
  % Variables iniciales
  a = intervalo(1);
  b = intervalo(2);
  df = diff(fs, 2);
  h = (b - a) / (numPuntos-1);
  
  % definicion del maximo de la funcion derivada
  dfaux=-1*abs(df);
  dfauxnum=matlabFunction(dfaux);
  x_max=fminbnd(dfauxnum,a,b);
  
  % calculo del alpha
  df_num = matlabFunction(df);
  alfa_max = abs(df_num(x_max));
  
  % Calculo del error
  error = (b-a)*h^2*alfa_max/12;
  
end