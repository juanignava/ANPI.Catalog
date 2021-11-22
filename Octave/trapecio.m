function trapecio()
  
  % Ejemplo numerico de la Regla del Trapecio
  
  clc;
  clear;
  
  f = 'ln(x)';
  a = 2;
  b = 5;
  
  [aprox error] = metodo_trapecio(f, a, b)
  
end

function [aprox error] = metodo_trapecio(f, a, b)
  
  % Esta funcion implementa la Regla del Trapecio, la cual permite aproximar
  % el resultado de una integral definida en un intervalo dado por medio de
  % la Segunda Formula de Newton-Cotes.
    
  % Parametros de entrada:    f es una cadena caracteres (string) que simboliza 
  %                           la ecuacion que sera evaluada en la integral
  %                           definida.
  
  %                           a hace referencia al valor inicial del intervalo
  %                           de la integral.
  
  %                           b corresponde al valor final del intervalo de la
  %                           integral.  
  
  % Parametros de salida:     aprox es el resultado de aproximacion de la
  %                           integral evaluada en el intervalo dado.
  %
  %                           error representa la cota de error de la Regla del
  %                           Trapecio.
  
  pkg load symbolic;
  
  syms x;
  
  f1 = matlabFunction(sym(f));
    
  d1 = matlabFunction(diff(f1,x));  
  
  d2 = matlabFunction(diff(d1,x));
  
  h = b - a;
  
  aprox = (h / 2) * (f1(a) + f1(b));
  
  % Verificacion del maximo del intervalo
  if abs(d2(a)) > abs(d2(b))
    
    e = a;     
    
  else
    
    e = b;
    
  end
  
  error = ((h^3)/12) * abs(d2(e));  
  
end
