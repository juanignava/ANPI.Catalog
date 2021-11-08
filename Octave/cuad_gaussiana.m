function cuad_gaussiana()
  
  % Ejemplo numerico para el metodo de Cuadratura Gaussiana y Cota de Error
  pkg load symbolic;
  clc;
  clear;
  f = 'ln(z)'
  a = 2
  b = 5
  n = 3
  
  [aprox error] = cuadratura_gaussiana(f,n,a,b)

end

function[aprox error] = cuadratura_gaussiana(f,n,a,b)
  
  % Esta funciion implementa el metodo de Cuadratura Gaussiana, el cuar permite
  % aproximar el resultado de una integral definida en un intervalo dado.
  
  % Parametros de entrada: f es una cadena de caracteres (string) que simboliza
  %                        la ecuacion que sera evaluada en la integral.
  
  %                        n corresponde al numero de orden
  
  %                        a corresponde al limite inferior de la integral
  
  %                        b corresponde al limite superior de la integral
  
  % Parametros de salida: aprox es el resultado de la aproximacion de la
  %                       integral en un intervalo dado.
  
  %                       error representa la cota de error del metodo de 
  %                       Cuadratura Gaussiana
  
  
  syms z
  
  fs=sym(f);
  
  % Obtencion de la funcion g(x) para poder utilizan el metodo de cuadratura 
  % gaussiana para intervalos distinto a [-1, 1].
  
  y=((b-a)*z+(b+a))/2;
  
  gs=(b-a)/2*subs(fs,z,y)
  
  gn=matlabFunction(gs);
  
  [x, w] = ceros_pesos_cuad_gauss(n)
  
  aprox = 0;
  
  % Obtencion de la aproximacion 
  
  for i=1:n
    aprox = aprox + w(i)*gn(x(i))
  endfor
  

  df = diff(gs,n);

  
  dfaux=-1*abs(df);
  dfauxnum=matlabFunction(dfaux);
  x_max=fminbnd(dfauxnum,a,b);
  
  df_num = matlabFunction(df);
  alfa_max = abs(df_num(x_max));
  
  
  error = alfa_max/135;
  
end

function [x, w] = ceros_pesos_cuad_gauss(n)
  
  if n>10 && n<2
    x=0; w=0;
    display('El valor de n debe ser menor o igual a 10 y mayor o igual a 2')
    
  else 
    
    x=[]; w=[];
    
    if n == 2
      x(1)=-0.5773502692; x(2)=0.5773502692;
      w(1)=1; w(2)=1;
      
    elseif n == 3
      x(1) = -0.7745966692; 
      x(2) = 0; 
      x(3) = 0.7745966692;
      w(1) = 0.5555555555; 
      w(2) = 0.888888888; 
      w(3) = w(1);
      
    elseif n == 4
      x(1) = -0.86113631159405;
      x(2) = -0.339981043584856;      
      x(3) = 0.339981043584856;      
      x(4) = -0.86113631159405;
      w(1) = 0.652145154862546;
      w(2) = 0.347854845137454;
      w(3) = w(1);
      w(4) = w(2);  
      
    elseif n == 5
      x(1) = -0.9061798459
      x(2) = -0.5384693101
      x(3) = 0
      x(4) = 0.9061798459
      x(5) = 0.5384693101
      w(1) = 0.2369268851
      w(2) = 0.4786286705
      w(3) = 0.5688888889
      w(4) = w(1)
      w(5) = w(2)
      
    elseif n == 6
      x(1) = -0.9324695142
      x(2) = -0.6612093865
      x(3) = -0.2386191861
      x(4) = 0.9324695142
      x(5) = 0.6612093865
      x(6) = 0.2386191861
      w(1) = 0.1713244924
      w(2) = 0.3607615730
      w(3) = 0.4679139346
      w(4) = w(1)
      w(5) = w(2)
      w(6) = w(3)
      
    elseif n == 7
      x(1) = -0.9491079123
      x(2) = -0.7415311856
      x(3) = -0.4058451514
      x(4) = 0
      x(5) = 0.9491079123
      x(6) = 0.7415311856
      x(7) = 0.4058451514
      w(1) = 0.1294849662
      w(2) = 0.2797053915
      w(3) = 0.3818300505
      w(4) = 0.4179591837
      w(5) = w(1)
      w(6) = w(2)
      w(7) = w(3)
      
    elseif n == 8
      x(1) = -0.9602898565
      x(2) = -0.7966664774
      x(3) = -0.5255324099
      x(4) = -0.1834346425
      x(5) = 0.9602898565
      x(6) = 0.7966664774
      x(7) = 0.5255324099
      x(8) = 0.1834346425
      w(1) = 0.1012285363
      w(2) = 0.2223810345
      w(3) = 0.3137066459
      w(4) = 0.3626837834
      w(5) = w(1)
      w(6) = w(2)
      w(7) = w(3)
      w(8) = w(4)
      
    elseif n == 9
      x(1) = -0.9681602395
      x(2) = -0.8360311073
      x(3) = -0.6133714327
      x(4) = -0.3242534234
      x(5) = 0
      x(6) = 0.9681602395
      x(7) = 0.8360311073
      x(8) = 0.6133714327
      x(9) = 0.3242534234
      w(1) = 0.0812743883
      w(2) = 0.1806481607
      w(3) = 0.2606106964
      w(4) = 0.3123470770
      w(5) = 0.3302393550
      w(6) = w(1)
      w(7) = w(2)
      w(8) = w(3)
      w(9) = w(4)

    elseif n == 10
      x(1) = -0.9739065285
      x(2) = -0.8650633667
      x(3) = -0.6794095683
      x(4) = -0.4333953941
      x(5) = -0.1488743390
      x(6) = 0.9739065285
      x(7) = 0.8650633667
      x(8) = 0.6794095683
      x(9) = 0.4333953941
      x(10) = 0.1488743390
      w(1) = 0.0666713443
      w(2) = 0.1494513492
      w(3) = 0.2190863625
      w(4) = 0.2692667193
      w(5) = 0.2955242247
      w(6) = w(1)
      w(7) = w(2)
      w(8) = w(3)
      w(9) = w(4)
      w(10) = w(5)
      
    end
  end
end

