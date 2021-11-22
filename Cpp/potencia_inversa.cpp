#include <iostream>
#include <armadillo>
#include <vector>
#include <tuple>
#include "matplotlibcpp.h"

using namespace std;
using namespace arma;
using namespace std;
namespace plt = matplotlibcpp;

tuple<double, mat> potencia_inversa(mat matriz, mat v_inicial, int iterMax, float tol) {

    /*
    Esta funcion aproxima el valor propio de menor magnitud y el vector propio
    asociado para una matriz especificada por medio del metodo de la potencia inversa.

    Parametros de entrada:
        matriz: es la matriz a la cual se le va a calcular el valor y vector propio.
        v_inicial: corresponde al vector inicial que utiliza el metodo de la potencia inversa.
        iterMax: es la cantidad maxima de iteraciones que se van a definir para el metodo.
        tol: es la tolerancia definida para el metodo.

    Parametros de salida:
        Ck: corresponde a la aproximacion del valor propio de menor magnitud de la matriz.
        xk: corresponde al vector propio correspondiente a la aproximacion.
    */

   // valores iniciales
   mat xk = v_inicial;
   mat yk = arma::solve(matriz, xk);
   double Ck = norm(yk, "inf");
   double Ck_1;
   double error = tol+1;
   vector<double> err;
   vector< double > iter;
   double iteracion;

    // Metodo iterativo
    for (int i = 0; i < iterMax; i++)
    {
        iteracion = i;

        // calculos del metodo de la potencia inversa
        xk = (1/Ck) * yk;
        yk = arma::solve(matriz, xk);
        Ck_1 = norm(yk, "inf");
            
        // calculo dele error
        error = abs(Ck_1-Ck);
        err.push_back(error);
        iter.push_back(iteracion);

        Ck = Ck_1;

        // condicion de parada
        if (error < tol) break;
    }

    // intrucciones para la graficacion
    plt::named_plot("Error |Ck_1 - Ck|", iter,err);
    plt::title("Error Potencia Inversa");
    plt::legend();
    plt::show();

    return {Ck, xk};

}

int main() {

    // Ejemplo numerico
    
    mat A = { {3, -1, 0},
              {-1, 2, -1},
              {0, -1, 3}};

    mat x0 = {1, 1, 1};
    x0 = x0.t();

    int iterMax = 100;
    float tol = 0.0000001;

    tuple<double, mat> result = potencia_inversa(A, x0, iterMax, tol);

    cout << "Valor propio de menor magnitud: " << get<0>(result) << endl;

    cout << "Vector propio asociaco" << get<1>(result) << endl;

}