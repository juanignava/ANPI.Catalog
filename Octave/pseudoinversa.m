function pseudoinversa() 
  
  % ejemplo numerico
  
  A = [1 2 4; 2 -1 1; 1 0 1];
  b = [0 7 0]';
  
  tol = 1e-5;
  iterMax = 1000;
  
  Xk = metodo_pseudoinversa(A, b, tol, iterMax)
  
end


function Xk = metodo_pseudoinversa(A, b, tol, iterMax)
  
  % metodo  
  
  n = length(A);
  Xk = (1/(norm(A)*norm(A))) * A';
  
  xk = Xk * b;
  
  iden = eye(n);
  
  for i=1:iterMax
    
    Xk = Xk * (2 * iden - A * Xk);
    
    xk_n = Xk * b;
    
    error = norm(xk_n - xk)/ norm(xk_n);
    
    xk = xk_n;
    
    if (error < tol) 
      break;
    endif
    
  endfor
  
end