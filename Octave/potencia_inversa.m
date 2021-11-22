function potencia_inversa()
  
  % ejemplo numerico
   
   A = [3 -1 0; -1 2 -1; 0 -1 3];
   
   x0 = [1 1 1]';
   
   iterMax = 1000;
   
   tol = 1e-5;
   
  [Ck xk] = met_potencia_inversa(A, x0, iterMax, tol)
  
end

function [Ck xk] = met_potencia_inversa(A, x0, iterMax, tol)
  
  % valores iniciales
  pkg load symbolic;
  syms x;
  
  xk = x0;
  yk = inverse(A)*xk; 
  Ck = norm(yk, inf);
  
  for i=1:iterMax
    
    xk = (1/Ck) * yk
    yk = A.^-1 *xk;
    Ck_1 = norm(yk, inf)
    
    error = abs(Ck_1 - Ck);
    
    Ck = Ck_1;
    
    if (error < tol)
      break;
    endif
    
  endfor
  
end