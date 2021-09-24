function metodo_gauss_seidel()
  
  clc; clear;
  A=[10 1 -1 5;0 -7 1 1; 1 2 5 0;4 4 1 10];
  b=[1 1 1 1]';
  iterMax=1000;
  tol=1e-5;
  dom = diag_dom(A);
  if dom ==1
    [xk k error] = gauss_seidel(A, b, tol, iterMax)  
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

function [xk k error] = gauss_seidel(A, b, tol, iterMax)
  n= length(A);
  
  D=diag(diag(A));
  U=triu(A)-D;
  L=tril(A)-D;
  
  y = sust_adelante(L+D, b);

  xk=ones(n,1);

  for k=1:iterMax
    
    z = sust_adelante(L+D, U*xk);
    xk = -z + y;

    error=norm(A*xk-b);
    if error<tol
      break
    end
  end
  
end

function x = sust_atras(A, b)
  
  n = length(b);  
  x = zeros(n,1);
  
  for i = n:-1:1 %sustitución hacia atrás  
    aux = 0;
    for j = i+1:n
      aux = aux + A(i, j)*x(j);
    end
     x(i) = (b(i) - aux)/A(i,i);  
  end  
  
end

function x = sust_adelante(A, b)
  
  n = length(b);  
  x = zeros(n,1);
  
  for i=1:n %sustitución hacia adelante
    aux = 0;
    for j=1:i-1
      aux = aux + A(i,j)*x(j);
    end
     x(i) = (b(i) - aux)/A(i,i); 
  end  
  
end