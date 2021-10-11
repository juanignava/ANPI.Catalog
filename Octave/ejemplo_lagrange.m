function ejemplo_lagrange()
  
  clc; clear;
  
  pkg load symbolic 
  xv=[-2 0 1];
  yv=[0 1 -1];

  p=metodo_lagrange(xv,yv)
  
end

function p=metodo_lagrange(xv,yv)
  syms x
  n=length(xv)-1;
  p=0;
  for k=0:n
    p=p+yv(k+1)*fun_Lk(xv,k);
  end
  p=expand(p);
end

function Lk=fun_Lk(xv,k)
  syms x
  %k=0,1,....,n
  n=length(xv)-1;
  Lk=1;
  for j=0:n
    if j~=k
      Lk=Lk*(x-xv(j+1))/(xv(k+1)-xv(j+1));
    end    
  end
  Lk=expand(Lk);
end