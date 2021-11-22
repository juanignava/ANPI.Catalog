function potencia()
  
  % Ejemplo numerico del Metodo de la Potencia
  
  clc;
  clear;
    
  A = [[3 -1 0];[-1 2 -1];[0 -1 3]];
  
  x0 = [1 1 1]';
  
  iterMax = 100;
  
  tol = 10^(-10);
    
  [xk err] = metodo_potencia(A, x0, tol, iterMax)
    
end

function [xk err] = metodo_potencia(A, x0, tol, iterMax)
  
  % Esta funcion implementa el Metodo de la Potencia, el cual permite aproximar
  % el modulo del valor propio de mayor magnitud de una matriz dada.
      
  % Parametros de entrada:    A es una matriz cuadrada cualquiera.
  
  %                           x es un vector inicial de valores aleatorios.
  
  % Parametros de salida:     aprox es la aproximacion de los valores del
  %                           modulo del valor propio.
  
  %                           vectoresPropios son los vectores que contienen
  %                           cada uno de los valores propios del Metodo de la
  %                           Potencia.
  
  xk = x0;
  
  error = [];
  
  vectoresPropios = [xk];
  
  for k = 0:iterMax
    
    yk = A * xk;
    
    ck = norm(yk, inf);
    
    xk_n = (1 / ck) * yk;
    
    err = norm(xk_n - xk);
    
    error = [error err];
    
    xk = xk_n;
    
    vectoresPropios = [vectoresPropios xk];
    
    % Condicion de parada
    if err < tol
      
        break;
        
    end   
    
  end
  
  % Instrucciones de graficacion
  plot(0:length(error)-1, error)
  title('Error absoluto vrs iteraciones')
  xlabel('Iteraciones (k)')
  ylabel('Error absoluto')
  
end