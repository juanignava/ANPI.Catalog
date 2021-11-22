function ejemplo_cholesky()
  
  A = [[4 2 1];[2 5 2];[1 2 6]];
  b = [1 2 3];
  
  x = fact_cholesky(A, b')
  
end


function pos = positiva_def(A)

  m = length(A);
  pos = 1;
  for k=1:m
    Ak = A(1:k, 1:k);
    deter = det(Ak);
    
    if abs(deter)< 10^-5
      pos = 0;
      break
    endif
  endfor

end

function equal = mat_equal(A, B)
 
 n = length(A);
 
 equal = 1;
 for i=1:n
   for j=1:n
     if A(i,j)!=B(i, j)
       equal = 0;
       break
     endif
   endfor
 endfor
 
end

function x = fact_cholesky(A, b)
  
  equal = mat_equal(A, A');
  
  if equal==1
    
    pos = positiva_def(A);
    
    if pos==1
      L = fact_cholesky_L(A);
      y = sust_adelante(L, b);
      x = sust_atras(L', y);
    else
      x = 'NA'
      display("La matriz no es positiva definida")
    endif
    
  else
    x = 'NA'
    display("La matriz no es simatrica")
  endif
  
end

function L=fact_cholesky_L(A)
  
  m=size(A,1);
  L=zeros(m,m);
  %Los valores de la columna 1 de L
  j=1;
  L(j,j)=sqrt(A(j,j));
  for i=2:m
    L(i,j)=A(i,j)/L(j,j);
  end
  %Los valores de las columna 2,...m, de L
  for j=2:m
    aux1=0;
    for k=1:j-1
      aux1+=(L(j,k))^2;
    end
    L(j,j)=sqrt(A(j,j)-aux1);
    
    for i=j+1:m
      aux2=0;
      for k=1:j-1
        aux2+=L(i,k)*L(j,k);
      end    
      L(i,j)=(A(i,j)-aux2)/(L(j,j));
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
