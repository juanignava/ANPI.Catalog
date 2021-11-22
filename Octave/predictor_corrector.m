function predictor_corrector()
  
  % Ejemplo numerico del Metodo Predictor-Corrector
  
  clc;
  clear;
  warning ("off");
  
  f = 'y-x^2+1';
  a = 0;
  b = 2;
  n = 11;
  y0 = 0.5;
  
  [xn yn polinomio] = metodo_predictor_corrector(f, a, b, n, y0)
    
end

function [xn yn polinomio] = metodo_predictor_corrector(f, a, b, n, y0)
  
  % Esta funcion implementa el Metodo de Predictor-Corrector, el cual permite
  % mejorar la estimacion de una pendiente por medio de la determinacion y
  % promediado de dos derivadas para un intervalo dado.
    
  % Parametros de entrada:    f es una cadena caracteres (string) que simboliza 
  %                           la pendiente a estimar.
  
  %                           a hace referencia al valor inicial del intervalo.
  
  %                           b corresponde al valor final del intervalo.
  
  %                           n es la cantidad de puntos del intervalo.
  
  %                           y0 valor inicial del eje y de coordenadas.
  
  % Parametros de salida:     xn corresponde a los n punto del intervalo en el
  %                           eje x de coordenadas.
  
  %                           yn hace referencia a la aproximacion de la
  %                           derivada de los n puntos en eje x de coordenadas.
  %
  %                           polinomio representa el polinomio de interpolacion
  %                           del Metodo de Predictor-Corrector.
  
  pkg load symbolic;
  
  syms x y;
  
  f1 = matlabFunction(sym(f));
  
  xn = [];
  
  yn = [y0];
      
  h = (b - a)/(n - 1);
  
  % calulo de los n puntos en x
  for i = 0: n - 1
    
    x = a + (i * h);
           
    xn = [xn x];    
    
  end 
    
  % calulo de los n puntos en y
  for i = 1: n - 1
    
    % valor predictor
    yn_vp = yn(i) + (h * f1(xn(i), yn(i)));
    
    num = f1(xn(i), yn(i)) + f1(xn(i + 1), yn_vp);
    
    % valor corrector
    yn_vc = yn(i) + h * (num / 2);
    
    yn = [yn yn_vc];
        
  end
  
  polinomio = lagrange(xn, yn);
  
  % Instrucciones de graficacion
  plot(xn, yn, 'r')
  title('Metodo de Predictor-Corrector')
  xlabel('xn')
  ylabel('yn')
    
end

function Lk_expand = fun_Lk(xv, k)
    
    % Calculo de los factores de variable de Lagrange. Se recibe el vector de
    % preimagenes y la iteracion actual en la cual se indefine el polinomio.
    % Esta corresponde a una funcion auxiliar de la funcion lagrange().
    
    syms x;
        
    n = length(xv);
    
    Lk = 1;
    
    for j = 1: n
      
        if j != k
          
            num = x - xv(j);
            
            den = xv(k) - xv(j);
            
            Lk = Lk * (num / den);
            
        end
            
    end
        
    Lk_expand = expand(Lk);
    
end

function pol = lagrange(xv, yv)
    
    % Esta funcion aproxima numericamente el polinomio de interpolacion que
    % pasa por los puntos (x0, y0), (x1, y1), ..., (xn, yn) dados en las
    % entradas por medio del metodo de Lagrange.

    % Parametros de entrada:    xv corresponde al arreglo contituido por las
    %                           preimagenes de los puntos considerados en el
    %                           polinomio.

    %                           yv hace referencia al arreglo constituido por
    %                           las imagenes de los puntos considerados en el
    %                           polinomio.

    % Parametros de salida:     pol representa el polinomio de forma simbolica
    %                           generado por medio de las diferencias divididas
    %                           de Newton.
    
    syms x;
     
    n = length(xv);
    
    p = 0;

    % Formacion del polinomio
    for k = 1: n
      
        p += yv(k) * fun_Lk(xv, k);
        
    end
        
    pol = expand(p);
    
end
