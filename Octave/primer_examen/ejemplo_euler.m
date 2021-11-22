%% Ejemplo del metodo de Euler
clc; clear; close all
a=0; b=5;
num_pt=11;
f=@(x,y) y-x.^2+1;
h=(b-a)/(num_pt-1);
xv=a:h:b;
yv=[0.5];
for n=1:num_pt-1  
  yv(n+1)=yv(n)+h*f(xv(n),yv(n));
end


%grafica de la aproximacion
hold on
plot(xv,yv,'r')

%solucion analitica
y_s=@(x) (x+1).^2-0.5*exp(x);
%Graficar la solucion
x_g=a:0.0001:b;
y_g=y_s(x_g);
plot(x_g,y_g,'b')

title('Graficas de la solucion analitica (azul) y por medio del metodo (rojo)')
xlabel('x')
ylabel('y(x)')