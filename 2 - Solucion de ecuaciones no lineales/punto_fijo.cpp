//Se importan las librerias a utilizar.
#include <iostream>
#include <math.h>
#include <vector>
#include "matplotlibcpp.h"
//Fin de la importacion de librerias.
using namespace std;
namespace plt = matplotlibcpp;


double funcion(double x){

  return log(2*x+1);

}


double punto_fijo(double x0, double tol, int iterMax){
  // Esta funcion aproxima numericamente el valor de 1/a, donde a es diferente de 0
  // Parametros de entrada: a=valor para aproximar su inverso multiplicativo, tol=tolerancia
  // Parametros de salida: y=aproximacion numerica de 1/a
    double error;
    double xk=x0;
    vector< double > errors;
    vector< double > iter;
    double iteracion;
    for(int k = 1; k<iterMax; k+=1){
      iteracion=1.0*k;
      double xk_n=funcion(xk);
      error=abs(xk_n-xk);
      xk=xk_n;
      errors.push_back(error);
      iter.push_back(iteracion);
      if (error<tol){
        break;
      }
    }

    cout << "xk=" << xk  << "\n";
    cout << "k=" << iteracion << "\n";
    cout << "error=" << error << "\n";

    plt::named_plot("Error |x_{k+1}-x_k|", iter,errors);
    plt::title("Error Punto Fijo");
    plt::legend();
    plt::show();

    return 0;
}

int main(){
  // Ejemplo Numerico
    double x0=7;
    double tol=0.0000001;
    int iterMax=100;
    double y = punto_fijo(x0, tol, iterMax);
    return 0;
}
