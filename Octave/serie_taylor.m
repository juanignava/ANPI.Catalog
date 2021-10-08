function serie_taylor()
  
  % ejemplo numerico serie de taylor
  
  clc; clear;
  
  f = 'exp(x)';
  tol = 10 ^-9;
  iterMax = 1000;
  x_var = 3;
  c = 1;
  
  [sk error] = metodo_serie_taylor(f, c, x_var, tol, iterMax)
  
end

function [sk error] = metodo_serie_taylor(f, c, x_var, tol, iterMax)
  
  error = tol + 1;
  
  pkg load symbolic;
  syms x;
  
  dk = matlabFunction(sym(f))
  sk = 0;
  
  
  for n = 0 : iterMax
    
    deriv = dk(c)
    
    sk0 = sk;
    
    sk = sk + (deriv/factorial(n))*(x_var - c)**n
    
    error = abs(sk - sk0);
    
    if error < tol
      break
    endif
    
    dk = matlabFunction(diff(dk, x));
    
  endfor
  
end