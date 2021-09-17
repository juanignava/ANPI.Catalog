% Ejemplo de evaluar un número a en la función exponencial, utilizando el
% polinomio de Taylor.
clc; clear;

tol = 1e-15;
a = 2;
y = 0;
Sk = 0;
e = [];

iterMax = 100000;

for n = 0: iterMax
  
  Sk_n = Sk + a^n / factorial(n); 
  error = abs(Sk_n - Sk);
  
  e = [e error];
  if error < tol
    
     break      
  end
  
  Sk = Sk_n;
end

plot(0:length(e)-1, e, "g", 'LineWidth', 2) 
title("Grado del polinomio vs error relativo")
xlabel("Grado del polinomio")
ylabel ("Error relativo")