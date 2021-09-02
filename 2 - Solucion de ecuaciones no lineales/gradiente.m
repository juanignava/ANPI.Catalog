function gradiente()
  
  % Ejemplo numérico del método GCNL
  
  clc; clear;
  
  f = '(x-2)^4 + (x-2*y)^2';
  x_0 = [0 3];
  tol = 10^-9;
  iterMax = 1000;
  
  [x_k error] = metodo_GCNL(f, x_0, tol, iterMax)
  
 end
 
 
 function [x_k error] = metodo_GCNL(f, x_0, tol, iterMax)
  % Esta funcion implementa el metodo del Gradiente Conjugado
  % No Lineal, el cual permite aproximar la optimizacion de una
  % funcion de forma iterativa.
  
  % Parametros de entrada: f -> es una cadena de caracteres (string)
  %                 que simboliza la funcion a optimiazar.
  
  %                 x_0 -> Representa el vector inicial de n dimen-
  %                 siones del metodo iterativo.
  
  %                 tol -> es la tolerancia del metodo, con la cual
  %                 se puede definir una condicion de parada
  
  %                 iterMax -> es la cantidad maxima de iteraciones
  %                 que se realizaran en el metodo
  
  % Parametros de salida: x_k -> corresponde al vector de n dimensiones
  %                 donde se ubica el punto optimizado de la funcion
  
  %                 error -> corresponde al error absoluto de la 
  %                 optimizacion |gradiente(f(x_k))|
  
  pkg load symbolic;
  syms x y;
  f1 = matlabFunction(sym(f));
  
  g = gradient(sym(f), [x,y]);  % gradiente simbolico
  g_n = matlabFunction(g);      % gradiente numerico
  
  % definiciones iniciales
  g_0 = g_n(x_0(1), x_0(2));
  d_0 = -g_0;
  x_k = x_0;
  d_k = d_0;
  g_k = g_0;
  error = 1;
  err = [];
  
  for n = 0:iterMax % Para encontrar el x_k
    sigma = 0.5;
    a = 1;
    for i = 0:iterMax % Para encontrar el alpha
      
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
    
    if error < tol % condicion de parada
      break;
    endif
    
    B = (norm(g_k)^2)/(norm(g_k_prev)^2);
    
    d_k(1) = -g_k(1) + B*d_k(1);
    d_k(2) = -g_k(2) + B*d_k(2);
    
  endfor
  
  % Instrucciones de graficacion
  eje_x = 0:length(err)-1;
  plot(eje_x,err,'b','LineWidth',3, 'Marker', '*')
  title('Error absoluto vrs iteraciones')
  xlabel('Iteraciones (k)')
  ylabel('Error absoluto')
  
 end
