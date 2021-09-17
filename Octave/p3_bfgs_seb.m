function ejemplo_newton_raphson()
  clc; clear; 
  %Calcular el cero de x^5 -2=0 mediante el metodo de Newton-Raphson
  %f="(x1 - 1.7*x2)*sin(x1) - 1.5*x3 - 0.1*x4*cos(x4 + x5 - x1) + 0.2*(x5)^2 -x2 -1";
  f= " x1^2 + x2^2 + x3^2 + x4^2 + x5^2"; %Funcion que se desea aproximar
  %f= " (x1^2 + x2^2 + x3^2 + x4^2 + x5^2)^2";
  variables = [ 'x', 'y'];
  %f = '(x1)**2 +(x2)**3 + (x3)**4 +(x4)**4 +(x5)**6'
  %f = '(x- 2)**2 + (y+3)**2 + x*y';
  x0 = 1; %Valor inicial
  tol=10^-5; %Tolerancia
  iterMax=1000; %Numero de iteraciones maximas
  [xk error]= newton_raphson(f, tol, iterMax)
  
end

function [xk error]= newton_raphson(f, tol, iterMax)
%Esta funcion aproxima la solucion de la ecuacion f(x)=0, 
%utilizando el metodo de la biseccion
    %
    %Sintaxis:  [xk k error]= newton_raphson(f, x0, tol, iterMax)
    % 
    %Parametros Iniciales: 
    %            f = una  cadena de caracteres (string) que representa a la funcion f
    %            x0 = valor inicial 
    %            tol = un numero positivo que representa a la tolerancia para el criterio |f(xk)|<tol
    %            iterMax = cantidad de iteraciones maximas
    %            
    %Parametros de Salida:                           
    %            xk = aproximacion del cero de la funcion f
    %            k = numero de iteraciones realizados
    %            error =  |f(xk)|
    %%%% Se debe instalar el paquete Symbolic
    %%%% Paso 1: Descargar el archivo symbolic-win-py-bundle-2.9.0.tar.gz
    %%%%         de la pagina https://github.com/cbm755/octsympy/releases
    %%%% Paso 2: Escribir en la Venta de Comandos de Octave la instruccion
    %%%%         pkg install symbolic-win-py-bundle-2.9.0.tar.gz  
    
    %%% Cargar el paquete symbolic
    warning('off', 'all');
    pkg load symbolic
    syms x1 x2 x3 x4 x5;
    variables = [x1 ; x2 ; x3 ; x4 ; x5];
    f=sym(f); 
    error = tol + 1;
    k = 0;
    err = [];
    
    

    
    %------------------------------------------------------------------------------------
    
    xk = [0.1; 0.1; 0.1; 0.1; 0.1];
    n = length(xk);
    lambda = 1;
    sigma1 = 0.5;
    Bk = eye(n,n);
    pk = 0;
   
 
    g = gradient(f,variables);
 
   
    while (k < iterMax) %Se verifica el punto de parada
      
      gk = double(subs(g, variables, xk));
      pk = double(inv(Bk)*-gk);
     
      f_izq = double(subs(f, variables, xk + lambda*pk));
      f_der = double(subs(f, variables, xk) + sigma1*lambda*transpose(gk)*pk);
      
      rho = 2;
      i = 0;
      while(f_izq <= f_der)
      
        lambda = rho^i;
        f_izq = double(subs(f, variables, xk + lambda*pk));
        f_der = double(subs(f, variables, xk) + sigma1*lambda*transpose(gk)*pk);
        i++; 
        
      endwhile
      
      xk1 = double(xk + lambda*pk);
      sk = double(xk1 - xk); 
      yk = double(subs(g, variables, xk1)- subs(g, variables, xk));
      
      sk_t = double(transpose(sk));
      yk_t = double(transpose(yk));
  
      cond = double(yk_t * sk);
      if (cond > 0)
        Bk = double(Bk - (Bk*sk*sk_t*Bk)/(sk_t*Bk*sk) + (yk*yk_t)/(yk_t*sk));
      endif
      
      xk = xk1;
      k += 1;
      error = double(norm(subs(g, variables, xk)));
      display(error);
      
      
      
      if (error < tol)
        break;
      endif
 

      
      
    
    
    endwhile

    xk = double(xk);
    %Graficacion
    plot(0:length(err)-1,err,'g','LineWidth',2)
    title('Grado del Polinomio vrs Error Relativo')
    xlabel('Grado del Polinomio (k)')
    ylabel('Error Relativo')
    
    
    
end
