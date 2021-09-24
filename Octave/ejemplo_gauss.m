function ejemplo_gauss()
  
  A = [[2 -6 12 16];[1 -2 6 6];[-1 3 -3 -7];[0 4 3 -6]];
  b = [70 26 -30 -26];
  
  x = gauss(A, b');
  
  C = [[1 0]; [2 3]];
  D = [2,3];
  
  x = sust_adelante(C, D')
  
end

function x = gauss(A, b)
  
  %Meotodo de eliminaion gaussiana para resolver Ax=b

  
  m = length(b);
  x = zeros(m,1); % genera una matriz de ceros, de tamano nxm
  
  A_b = [A b]; % genera la matriz aumentada
  
  for k=1:m % recorriendo todas las columnas menos la ultima
    
    for i = k+1:m % recorre por cada fila por debajo de la diagonal
      
      fact = A_b(i,k)/A_b(k,k); % factor que genera el ceros
      
      for j=k:m+1 % se le hace el calculo a todas las columnas
        A_b(i, j) = A_b(i,j)- fact*A_b(k,j);
      endfor
      
    endfor
    
  endfor
  
  A_eq = A_b(:,1:m); % tomar de la fila 1 a la n
  B_eq = A_b(:,m+1); % tomar solo la fila n 
  
  % ** OJO en octave la fila 1 es la primera, no es la cero
  
  x = sust_atras(A_eq, B_eq);
  
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