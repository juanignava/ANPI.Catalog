function p3_bfgs()
  
  % Ejemplo numerico del metodo BFGS
  
  clc; clear;
  
  f = '(x-2)^4 + (x-2*y)^2';
  x0 = [0 0];
  B0 = [5, 1; 1, 4];
  tol = 10^-9;
  iterMax = 1000;

  [k xk error] = bfgs_method(f, x0, B0, tol, iterMax)
  
end

function [k xk error] = bfgs_method(f, x0, B0, tol, iterMax)
  
  pkg load symbolic;
  syms x y;
  f1 = matlabFunction(sym(f));
  
  sig1 = 0.3;
  sig2 = 0.6;
  alpha = 0.2;
  epsilon = 0.21;
  
  Bk = B0;
  g = gradient(sym(f), [x y]);
  g_n = matlabFunction(g);
  g0 = g_n(x0(1), x0(2));
  gk = g0;
  xk = x0';
  
  error = 1;
  err = [];

  for k=0:iterMax
    pk = inv(Bk)*-gk;
    
    n=1;
    lam_k = 1;
    for i=0:iterMax
      
      lam_k = sig2^i;
      i=i+1;
      izq = f1(xk(1)+lam_k*pk(1), xk(2)+lam_k*pk(2));
      sum1 = f1(xk(1), xk(2));
      sum2 = sig1*lam_k*gk'*pk;
      der = sum1 + sum2;
      if izq <= der
        break;
      endif
    endfor
    
    xk_prev = xk;
    xk = xk + lam_k*pk;
    gk_prev = gk;
    gk = g_n(xk(1), xk(2));
     
    sk = xk - xk_prev;
    yk = gk - gk_prev;
    
    izq1 = (yk'*sk)/(norm(sk)^2);
    der1 = epsilon*(norm(gk))^alpha;
    
    if izq1 >= der1
      Bk = Bk - (Bk*sk*sk'*Bk)/(sk'*Bk*sk) + (yk*yk')/(yk'*sk);
    endif
    
    error = norm(gk)
    err = [err error];
    
    if error < tol % condicion de parada
      break;
    endif
    
    k = k + 1;

  endfor
  
  eje_x = 0:length(err)-1;
  plot(eje_x,err,'b','LineWidth',3, 'Marker', '*')
  title('Error absoluto vrs iteraciones')
  xlabel('Iteraciones (k)')
  ylabel('Error absoluto')
  
end