function ejemplo_thomas()
  
  A = [[1 4 0 0];[3 4 1 0];[0 2 3 1];[0 0 1 3]];
  b = [1 2 3 4];
  
  x = thomas(A, b')
  
end

function x = thomas(A, d)

  n = length(d);
  p = zeros(n-1, 1);
  q = zeros(n, 1);
  
  a=zeros(n-1,1);
  b=zeros(n,1);
  c=zeros(n-1,1);
  
  for i=1:n % define a, b and c vetors
  
    for j=1:n
      if i<n
        c(i) = A(i, i+1);
      endif
      if i>1
        a(i) = A(i, i-1);
      endif
      b(i) = A(i,i);
    endfor
    
  endfor
  
  for i=1:n % define p and q vectors
    
    if i==1
      p(i) = c(i)/b(i);
      q(i) = d(i)/b(i);
    else
      if i!=n
        p(i) = c(i)/(b(i)-p(i-1)*a(i));
      endif
      
      q(i) = (d(i)-q(i-1)*a(i)) / (b(i)-p(i-1)*a(i));
    endif
   
  endfor
  
  for i=n:-1:1 % crear el vector x
    if i==n
      x(i) = q(i);
    else
      x(i) = q(i)-p(i)*x(i+1);
    endif
  endfor
  
  
end