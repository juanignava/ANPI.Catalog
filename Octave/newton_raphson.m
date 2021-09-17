function newton_raphson()
  
  % Ejemplo numerico del metodo de Newton Raphson
  
  clc; clear;
  
  f = 'exp(x)+x+2';
  x0 = 10;
  tol = 10^-9;
  iterMax = 1000;
  
  [xk error] = metodo_newton_raphson(f, x0, tol, iterMax)
  
end

function [xk error] = metodo_newton_raphson(f, x0, tol, iterMax)
  
  % Esta funcion implementa el metodo de Newton-Raphson, el cual, permite
  % aproximar la solucion de una ecuacion no lineal al aproximar el cero de
  % la funcion sin tener que utilizar el Teorema de Bolzano.
  
  % Parametros de entrada:    f es una cadena caracteres (string) que simboliza 
  %                           la ecuacion no lineal a aproximar sel cero.
  
  %                           x0 es el valor inicial para evaluar la funcion
  
  %                           tol representa la tolerancia del metodo la cual, 
  %                           nos permite tener una condicion de parada
  
  %                           iterMax es la cantidad maxima de iteraciones para 
  %                           evaluar la funcion
  
  
  % Parametros de salida:      xk es el resultado de aproximacion del cero 
  %                           de la funcion
  %
  %                           error es error relativo del cero aproximado
  
  e = [];
  error = tol+1;
  
  pkg load symbolic;
  syms x;
  
  f1 = matlabFunction(sym(f));
  d= matlabFunction(diff(f1,x));
 
  xk = x0;
  for n = 0: iterMax
    if d(xk) < tol
      error = 'NA';
      xk = 'NA';
    endif
    xk = xk - f1(xk)/d(xk);
    error = abs(f1(xk));
    e = [e error];
    
    if error<tol
     break
    endif
     
  end
  
  plot(0:length(e)-1, e)
  title('Error relativo vrs iteraciones')
  xlabel('Iteraciones (k)')
  ylabel('Error relativo')
  
end