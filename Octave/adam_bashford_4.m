function adam_bashford_4 ()
  
  % Ejemplo numerico
  clc; clear;
  f = '1 + (x-y)^2';
  inter = [2 4];
  num_puntos = 6;
  val_inicial = [1 1.191 1.5964 1.8883]; % puntos calculados inicialmente con el
                                         % metodo de Euler
  [xv, yv, pol] = metodo_adam_bashford_4(f, inter, num_puntos, val_inicial)
  
  % Grafica de la solucion por medio del metodo (rojo)
  hold on
  x_g=inter(1):0.0001:inter(2);
  pol1 = matlabFunction(sym(pol));
  y_p=pol1(x_g);
  plot(x_g,y_p,'r')
  
  % Grafica de la solucion analitica del problema (azul)
  y_s=@(x) -1*(x-1).^(-1)+x
  x_g=inter(1):0.0001:inter(2);
  y_g=y_s(x_g);
  plot(x_g,y_g,'b')
  title('Graficas de la solucion analitica (azul) y la del polinomio de interpolacion por medio del metodo (rojo)')
  xlabel('x')
  ylabel('y(x)')
  
end

function [xv, yv, pol] = metodo_adam_bashford_4(f, inter, num_puntos, val_inicial)
  
  % Esta funcion calcula la solucion al problema de Cauchy de la forma
  % y'(x) = f(x, y(x))
  % y(x0) = y0
  % por medio del metodo de adam bashford
  % Parametros de entrads:
  %   f: la funcion f definida anteriormente.
  %   inter: es el intervalo que contiene el valor inicial y final
  %     de analisis de la ecuacion.
  %   num_puntos: corresponde al numero de puntos con el que se
  %     realizara el analisis.
  %   val_inicial: corresponde al conjunto de valores iniciales
  %     que necesita este metodo para ejecutarse
  
  % Parametros de salida:
  %   xv, yv: son lo vectores que componen los pares ordenados de los
  %     valores aproximados de la solucion del problema.
  %   pol: es el polinomio de interpolacion calculado por medio del 
  %     metodo de lagrange para los puntos indicados.
  
  % valores iniciales
  pkg load symbolic;
  syms x y;
  f1 = matlabFunction(sym(f));
  
  % calculo de h
  a = inter(1);
  b = inter(2);
  
  h=(b-a)/(num_puntos-1);
  xv=a:h:b;
  
  % calculo de yv
  y0 = val_inicial(1);
  y1 = val_inicial(2);
  y2 = val_inicial(3);
  y3 = val_inicial(4);
  yv = [y0 y1 y2 y3];

  for n=4:num_puntos-1  
    fk = f1(xv(n), yv(n));
    fk_1 = f1(xv(n-1), yv(n-1));
    fk_2 = f1(xv(n-2), yv(n-2));
    fk_3 = f1(xv(n-3), yv(n-3));
    yv(n+1)= yv(n)+ (h/24) * (55*fk - 59*fk_1 + 37*fk_2 - 9*fk_3);
  end
  
  % uso de lagrange para calcular polinomio de interpolacion
  pol = metodo_lagrange(xv, yv);
  
end

function p=metodo_lagrange(xv,yv)
  
  % Metodo de Lagrange para el caclculo del polinomio
  % de interpolacion
  
  % Parametros de entrada:
  %   xv, yv: vetores de los pares ordenador a partir de
  %     los cuales se calculara el polinomio.
  
  % p =  polinomio de interpolacion de forma simbolica.
  
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
  n=length(xv)-1;
  Lk=1;
  for j=0:n
    if j~=k
      Lk=Lk*(x-xv(j+1))/(xv(k+1)-xv(j+1));
    end    
  end
  Lk=expand(Lk);
end