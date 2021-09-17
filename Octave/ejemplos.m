% ejemplo con varias funciones en un mismo archivo


function ejemplos()
  [Sk error]=exp_mac(2, 1e-5)
  z=suma(2, 3)
 end
 
 function [Sk error]=exp_mac(a, tol)
  % Explicacion funcion
  % a=
  % tol = 
  
  Sk = 0;
  iterMax = 100000;

  for n = 0: iterMax
    
    Sk_n = Sk + a^n / factorial(n); 
    error = abs(Sk_n - Sk); 
    if error < tol
      
       break      
    end
    
    Sk = Sk_n;
  end
  
 end
 
 function y = suma(a, b)
   y = a+b;
 end
