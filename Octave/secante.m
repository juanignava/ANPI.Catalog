
function secante()
  
  % ejemplo del metodo de la secante
  
  clc; clear;
  
  f = 'exp(x)+x+2';
  x0 = 0;
  x1 = 1;
  tol = 10^-9;
  iterMax = 1000;
  
  [xk_1 error] = metodo_secante(f, x0, x1, tol, iterMax)
  
end

function [xk_1 error] = metodo_secante(f, x0, x1, tol, iterMax)
  
  pkg load symbolic;
  syms x;
  
  f1 = matlabFunction(sym(f));
  
  xk_0 = x0;
  xk_1 = x1;
  
  for k=0:iterMax
    
    n = (xk_1 - xk_0) * f1(xk_1);
    d = f1(xk_1) - f1(xk_0);
    
    xk_0 = xk_1;
    
    xk_1 = xk_1 - n/d;
    
    error = abs(f1(xk_1))
    
    if error < tol
      break
    endif
    
  endfor
  
end