
function cota_poly_inter()
  
  % Ejemplo numerico
  clc; clear;
  f = 'sin(pi*x/2)';
  xv = [-1 0 1 2];
  num_eval = 0.54;
  error = cota_poly(f, xv, num_eval)
end

function error = cota_poly(f, xv, num_eval)
  % Esta funcion es utilizada paraa el calculo de la 
  % cota de error de un polinomio de interpolacion
  % definido para una funcion f y un conjunto soporte xv
  
  % Parametros de entrada:
  %   f: corresponde a la funcion de donde se tomaron los
  %         puntos del calculo de la interpolacion.
  %   xv: corresponde al conjunto sooprte sobre el cual se
  %         realiza la interpolacion.
  %   num_eval: corresponde al numero sobre el cual se esta
  %         realizando el calculo de la cota de interpolacion
  %         en la funcion.
  
  % Parametros de salida:
  %   error: es el valor numerico de la cota de error calculada.
  
  % definicion de constantes iniciales
  pkg load symbolic;
  syms x;
  n = length(xv);
  fs = sym(f);
  
  % calculo de la n-esima derivada de la funcion
  % donde n es la cantidad de elementos en el conjunto soporte
  df = diff(fs, n);
   
  a=xv(1); b=xv(end); %Valores del intervalo
  
  % definicion del maximo de la funcion derivada
  faux=-1*abs(fs);
  fauxnum=matlabFunction(faux);
  x_max=fminbnd(fauxnum,a,b);
  
  % Calculo de la constante del error alpha/n!
  fdn = matlabFunction(df);
  alpha = fdn(x_max);
  const = (alpha/factorial(n));
  
  % definicion del polinomio por el que se multiplica
  % la constante para encontrar el error
  var = 1;
  for i=1:n
    var = var * (x - xv(i));
  endfor
  var_num = matlabFunction(var);
  
  % calculo de la cota de error para el valor num_eval
  error = const * var_num(num_eval);
  
end