// Se importan las librearías a utilizar
#include <iostream>
#include <math.h>
#include <vector>
//#include "matplotlibcpp.h"

using namespace std;
//namespace plt = matplotlibcpp;

// Definicion de la funcion del problema
double function(double x){

    return exp(x)+ pow(x,2) -6;

}

double falsa_posicion(double a, double b, double tol, int iterMax) {
    /*
    Esta  funcion  aproxima  numericamente  la  solucion  de una  ecuacion no
    lineal  por  medio  del  metodo  de la falsa posicion  para un  intervalo  dado

    Parametros de entrada: a, b = valores que limitan el intervalo de analisis 
                            [a, b]
                           tol = valor determinado como la tolerancia mínima del
                            error aplicado en el metodo.
                           iterMax = iteraciones maximas que se le aplicaran al 
                           metodo iterativo

    Parametros de salida: y = valor de la aproximacion del cero de la funcion dentro 
                          del intervalo (solo si se cumple el teorema de Bolzano)
                          iteracion = cantidad de iteraciones necesarias para obetener
                          el resultado.
                          error = error absoluto obtenido con dicha iteracion.
    */

   double error;
   double xk;
   vector < double > errors;
   vector < double > iter;
   double iteracion;
   if (function(a)*function(b) < 0) { // Se verifica teorema de Bolzano
   
        for (int k = 0; k < iterMax; k++) {
       
            iteracion = 1.0*k;
            xk = b - ((b-a)/(function(b)-function(a)))*function(b);

            if (function(a)*function(xk) < 0) {
                b = xk;
            }
            else {
                a = xk;
            }

            error = abs(function(xk));
            errors.push_back(error);
            iter.push_back(iteracion);

            if (error < tol) {
                break;
            }
            
        } 
        cout << "xk=" << xk  << "\n";
        cout << "k=" << iteracion << "\n";
        cout << "error=" << error << "\n";
   }
   else {
       cout << "No se cumle el teorema de Bolzano" << "\n";
   }

    /*
    plt::named_plot("Error |function(x)|", iter, errors);
    plt::title("Error Falsa Posicion");
    plt::legend();
    plt::show();
    */
    return xk;
}


int main(){

    // Ejemplo numerico
    double a = 1;
    double b = 3;
    double tol = 0.00001;
    int iterMax = 1000;
    double y = falsa_posicion(a, b, tol, iterMax);
}