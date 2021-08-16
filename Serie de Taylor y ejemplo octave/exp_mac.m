% aspectos importantes, no existe un return, màs bien se usan los paréntesis
% luego de la palabra funtion para definir cuáles son mis salidas
% el nombre de la funciño tiene que coincidir con el nombre del archivo

% desde la consola lo que hago para usar la funcion es escribir lo siguiente
% exp_mac(2, 1e-5)
% esto solo me va a dar uno de los valores de salida (el primero Sk_n)
% si quiero que me indique todos entonces debo escribir esto en la consola
% [val1 val2] = exp_mac(2, 1e-5)

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