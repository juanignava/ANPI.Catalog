function metodo_jacobi()
  
  clc; clear;
  A=[10 1 -1 5;0 -7 1 1; 1 2 5 0;4 4 1 10];
  b=[1 1 1 1]';
  iterMax=1000;
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
  LpU=A-diag(d);
  z=D_inv*b;

  xk=ones(4,1);

  for k=1:iterMax
    xk=-D_inv*LpU*xk+z;
    error=norm(A*xk-b);
    if error<tol
      break
    end
  end
  
end



