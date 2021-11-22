function ejemplo_biseccion()
  clc; clear;
  % Calcular el cero de exp(x)-2*x-10=0
  % Paso 1: Conocer el intervalo donde se encuentra el cero
  %         Sugerencia: graficar la función (podemos usar el comando plot que
  % implementamos abajo)
  %         f(x) =  exp(x)-2*x-10
  
  %xv = -5:0.1:5;
  %yv = exp(xv)-2*xv-10;
  %plot(xv, yv)
  %grid on
  
  % Concluimos que la funci{on tiene un cero en el intervalo [2,4]
  
  f = 'exp(x)-2*x-10';
  a=2;
  b=4;
  tol=10^-5;
  iterMax =  1000;
  [xk k error]=biseccion(f, a, b, tol, iterMax)
 end
 
 
 function [xk k error] = biseccion(f, a, b, tol, iterMax)
   % Esta funcion aproxima la solucion de la ecuacion f(x)=0,
   % utilizando el método de la bisección
   %
   % Sintaxis: [xk k error] = biseccion(f, a, b, tol, iterMax)
   %
   % Parámetros iniciales:
   %    f = una cadena de caracteres (string) que representa a la funcion
   %    a,b = son los extremos del intervalo [a,b]
   %    tol = un numero positivo que representa la tolerancia.
   %    iterMax = cantidad de iteraciones máximas
   %
   % Parámetros de salida:
   %    x= aproximación del cero de la función f
   %    k = la cantidad de iteraciones
   %    error = |f(xk)|
   
   %%% Se debe instalar el paquete Symbolic
   %  Paso 1: Descargar el archivo symbolic-win-py-bundle-2.9.0.tar.gz
   %   de la pàgina https://github.com/cbm755/octsympy/releases
   % Paso2: Escribir en la Ventana de Comandos la instruccion
   %   pkg install symbolic-win-py-bundle-2.9.0.tar.gz
   % Paso3 :cargar el paquete cargado 
   % pkg load symbolic

    pkg load symbolic; % IMPORTANTE PONER PUNTO Y COMA AQUI
    f1 = matlabFunction(sym(f)); % importante usar esto
    % la funcion f1 a la cual nosotros vamos a comparar el metodo de la biseccion
    
    if f1(a)*f1(b)< 0
      % Se cumple la condicion de bolzano
      for k=1:iterMax
        xk = (a+b)/2; % tenemos dos intervalos
        
        if f1(a)*f1(xk) < 0 % se cumple la condicion en el intervalo 1
          b = xk;
          
        else % se cumple la condicion en el intervalo 2
          a = xk;
        end
        
        error = abs(f1(xk));
        
        if error < tol
          break;
        endif
      end
      
    else
      
      xk = 'NA';
      k = 'NA';
      error = 'NA';
      display ('El intervalo seleccionado no cumple las condiciones del teorema de Bolzano');
      
    end
    
    
end