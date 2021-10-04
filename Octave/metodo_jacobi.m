function metodo_jacobi()
  
  clc; clear;
  A=[10 -2 1 -1;1 -20 1 1; 0 4 12 -6; -1 3 3 -8];
  b=[1 1 1 1]';
  iterMax= 3;
  tol=1e-5;
  dom = diag_dom(A);
  if dom ==1
    [xk k error] = jacobi(A, b, tol, iterMax)  
  else
    xk = 'NA'
    k= 'NA'
    error = 'NA'
  endif
  
  
  display('La solcion delmetofdo es')
  inv(A)*b
  
end

function dom = diag_dom(A)
  
  n = length(A);
  dom = 1;
  for i=1:n
    diag_el = A(i,i);
    sum_fila = 0;
    for j=1:n
      if i!=j
        sum_fila = sum_fila + abs(A(i,j));
      endif
    endfor
    
    if abs(diag_el) < sum_fila
      dom =0;
      break
    endif
  endfor
  
end

function [xk k error] = jacobi(A, b, tol, iterMax)
  
  d=diag(A);
  D_inv=diag(1./d);
  display("L + U")
  LpU=A-diag(d)
  display("D inversa por B")
  z=D_inv*b
  

  xk=ones(4,1)

  for k=1:iterMax
    
    display("Calculo de x: ")
    xk=-D_inv*LpU*xk+z
    
    display("operacion de sistema")
    opera= A*xk-b
    
    display("Error con norma")
    error=norm(A*xk-b)
    if error<tol
      break
    end
  end
  
end



