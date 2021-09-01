function gradiente()
  
  % Ejemplo numérico del método GCNL
  clc; clear;
  f = '(x-2)^4 + (x-2*y)^2';
  x_0 = [0 3];
  tol = 10^-5;
  iterMax = 1000;
  
  [x_k error] = metodo_GCNL(f, x_0, tol, iterMax)
  
 end
 
 
 function [x_k error] = metodo_GCNL(f, x_0, tol, iterMax)
  % implementacion del GCNL
  
  pkg load symbolic;
  syms x y;
  f1 = matlabFunction(sym(f));
  g = gradient(sym(f), [x,y]);
  g_n = matlabFunction(g);
  g_0 = g_n(x_0(1), x_0(2));
  d_0 = -g_0;
  x_k = x_0;
  d_k = d_0;
  g_k = g_0;
  error = 1;
  err = [];
  
  for n = 0:iterMax
    sigma = 0.5;
    a = 1;
    for i = 0:iterMax
      
      izq = f1(x_k(1) + a*d_k(1), x_k(2) + a*d_k(2)) - f1(x_k(1), x_k(2));
      der = sigma*a*[g_k(1)*d_k(1)+g_k(2)*d_k(2)];
      
      if izq <= der
        break;
      endif
      
      a = a/2;
      
    endfor
    
    x_k(1) = x_k(1) + a*d_k(1);
    x_k(2) = x_k(2) + a*d_k(2);
    g_k_prev = g_k;
    g_k = g_n(x_k(1), x_k(2));
    
    error = norm(g_k);
    err = [err error];
    if error < tol
      break;
    endif
    
    B = (norm(g_k)^2)/(norm(g_k_prev)^2);
    
    d_k(1) = -g_k(1) + B*d_k(1);
    d_k(2) = -g_k(2) + B*d_k(2);
    
  endfor
  
  plot(0:length(err)-1, err)
  title('Error absoluto vrs iteraciones')
  xlabel('Iteraciones (k)')
  ylabel('Error absoluto')
  
 end
