%Ejemplo para calcular la cota de error
%del trazafor cúbico
clc; clear;
f='exp(x/2)';
xv=[1 1.5 1.75 2.15 2.4 3]; %Intervalo [1, 3]

%Cálculo de h
n=length(xv);
dist=zeros(n-1,1);
for i=1:n-1
  dist(i)=xv(i+1)-xv(i);
end
h=max(dist)

%Calculo del max de la 4ta derivada
pkg load symbolic
syms x

fs=sym(f);
fs4=diff(fs,4);
a=xv(1); b=xv(end); %Valores del intervalo
faux=-1*abs(fs4);
fauxnum=matlabFunction(faux);
x_max=fminbnd(fauxnum,a,b);

fn4=matlabFunction(fs4);

cota=(5*h^4/384)*fn4(x_max)