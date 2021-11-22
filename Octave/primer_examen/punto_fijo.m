function ejemplo_numerico_punto_fijo()
  
  % Ejemplo numerico 

  clc; clear;
  f = 'ln(2*x+1)';
  x0 = 1.5;
  tol = 10^-5;
  iterMax = 1000;
  
  [aprox error] = punto_fijo(f, x0, tol, iterMax)

end

function [aprox error] = punto_fijo(f, x0, tol, iterMax)
  
  % Esta funcion aproxima numericamente el valor de convergencia del punto fijo
  % de una ecuacion no lineal que cumpla el criterio de existencia (Teorema de 
  % Brouwer)y unicidad.
  %  
  % Parametros de entrada: f es una cadena caracteres (string) que simboliza la
  %                        ecuacion no lineal a aproximar su respectivo punto
  %                        fijo.
  %
  %                        x0 corresponde al valor inicial de la aproximacion.
  %
  %                        tol representa el valor de la tolerancia del error
  %                        relativo.
  %
  %                        iterMax hace referencia a la cantidad maxima de
  %                        iteraciones que realizara el metodo para aproximar
  %                        el punto fijo de la funcion dada
  %  
  % Parametros de salida: aprox es la aproximacion del punto fijo de la ecuacion
  %                       no lineal.
  %
  %                       error corresponde al error relativo del punto fijo
  %                       aproximado.
                                      
  pkg load symbolic % Importar el paquete simbolico
    
  f1 = matlabFunction(sym(f));
  
  x_k = x0;
  
  e = []; % Arreglo de errores relativos
  
  for k = 0 : iterMax % Iteraciones del metodo
            
    x_k1 = f1(x_k);
    
    aprox = x_k1;
    
    error = abs(f1(x_k1) - f1(x_k)) / abs(f1(x_k1));
    
    e = [e error];
    
    if error < tol % Condicion de parada del ciclo iterativo
          
      break;
      
    end
   
    x_k = x_k1; 
    
  end 
  
  % Configuracion del grafico de Iteraciones vs Error Relativo  
  x = 0 : length(e) - 1;
  plot(x,e,'b','LineWidth',3, 'Marker', '*')
  title('Gráfico de Iteraciones vs Error Relativo del Método del Punto Fijo')
  xlabel('Iteraciones')
  ylabel('Error Relativo')
  
end
