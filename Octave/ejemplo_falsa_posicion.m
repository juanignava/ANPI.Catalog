function ejemplo_falsa_posicion()
  clc;
  clear;
  % Calcular el cero de exp(-x)=ln(x)
  % Usndo el método de la falsa posicion
  
  f = 'exp(-x)-ln(x)';
  a = 1;
  b = 2;
  tol = 10^-9;
  iterMax = 1000;
  
  [xk k error] = falsa_posicion(f, a, b, tol, iterMax)
  
end

function [xk k error] = falsa_posicion(f, a, b, tol, iterMax)
  
  % Parámetros iniciales niciales: ----
  % Parámetros de salida: --- 
  
  pkg load symbolic;
  f1 = matlabFunction(sym(f));
    
  if f1(a)*f1(b) < 0 % Si se cumple el teorema de Bolzano
    
    for k=1:iterMax
      
      xk = b - ((b-a)/(f1(b)-f1(a)))*f1(b);
      
      if f1(a)*f1(xk) < 0   
        b = xk;
      else
        a = xk;
      end
      
      error = abs(f1(xk));
      if error < tol
        break;
      end
      
    end
  
else
  xk= 'NA';
  k='NA';
  error = 'NA'
  display('El intervalo no cumple con el teorema de Bozano');
  
  end
  
end
  