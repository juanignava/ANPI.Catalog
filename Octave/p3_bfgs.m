function p3_bfgs()
  
  % Ejemplo numerico del metodo BFGS
  
  clc; clear;
  
  %f = '(a+1.7*b)*sin(a) - 1.5*c - 0.1*d*cos(d+h-a) + 0.2*h^2-b-1';
  %f = 'a^2+b^2+c^2+d^2+e^2';
  f = 'a^4 + 2*b^4 + 3*c^4 + 4*d^4 + 5*e^4'
  %x0 = [1.12 -0.8 0.5 -1.2 0.1];
  x0 = rand(1,5)*2.56-1.28
  
  
  %B0 = [5, 1, 3, 2, 6; 1, 4, 7, 9, 8; 3, 7, 1, 2, 3; 2, 9, 2, 8, 6; 6, 8, 3, 6, 5];
  B0 = [1, 0, 0, 0, 0; 0, 1, 0, 0, 0; 0, 0, 1, 0, 0; 0, 0, 0, 1, 0; 0, 0, 0, 0, 1];
  tol = 10^-5;
  iterMax = 10000;

  [k xk error] = bfgs_method(f, x0, B0, tol, iterMax)
  
end

function [k xk error] = bfgs_method(f, x0, B0, tol, iterMax)
  
  pkg load symbolic;
  syms a b c d e;
  f1 = matlabFunction(sym(f));
  
  sig1 = 0.3;
  sig2 = 0.6;
  alpha = 0.2;
  epsilon = 0.21;
  
  Bk = B0;
  g = gradient(sym(f), [a, b, c, d, e]);
  g_n = matlabFunction(g);
  g0 = g_n(x0(1), x0(2), x0(3), x0(4), x0(5));
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
      izq = f1(xk(1)+lam_k*pk(1), xk(2)+lam_k*pk(2), xk(3)+lam_k*pk(3), xk(4)+lam_k*pk(4), xk(5)+lam_k*pk(5));
      sum1 = f1(xk(1), xk(2), xk(3), xk(4), xk(5));
      sum2 = sig1*lam_k*gk'*pk;
      der = sum1 + sum2;
      if izq <= der
        break;
      endif
    endfor
    
    xk_prev = xk;
    xk = xk + lam_k*pk;
    gk_prev = gk;
    gk = g_n(xk(1), xk(2), xk(3), xk(4), xk(5));
     
    sk = xk - xk_prev;
    yk = gk - gk_prev;
    
    izq1 = (yk'*sk)/(norm(sk)^2);
    der1 = epsilon*(norm(gk))^alpha;
    
    if izq1 >= der1
      Bk = Bk - (Bk*sk*sk'*Bk)/(sk'*Bk*sk) + (yk*yk')/(yk'*sk)
    endif
    
    error = norm(gk);
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