// Se importan las librerias a utilizar
#include <iostream>
#include <math.h>
#include <vector>
#include "matplotlibcpp.h"

using namespace std;
namespace plt = matplotlibcpp;

// Definicion de la funcion del problema
double funcion(double x) {

    return exp(-pow(x,2)) - x; // e^(-(x^2)) - x = 0

}

double secante(double x0, double x1, double tol, int iterMax) {

    /*
        Esta funcion aproxima numericamente la solucion de una ecuacion no lineal
        por medio del metodo de la secante partiendo de dos valores iniciales hasta
        encontrar el cero de la funcion en cuestion.

        Parametros de entrada: x0 corresponde al primer valor inicial necesario para
                               comenzar el metodo iterativo.

                               x1 representa el segundo valor inicial necesario para
                               comenzar el metodo iterativo.

                               tol hace referencia al valor de tolerancio del error
                               absoluto.

                               iterMax es la cantidad maxima de iteraciones que realizara
                               el metodo para aproximar el cero de la funcion dada.

        Parametros de salida: aprox corresponde a la aproximacion del cero de la funcion.

                              error representa el error absoluto del cero aproximado.
                                
    */

    vector < double > iterVector; // Vector de iteraciones
    vector < double > errorVector; // Vector de errores

    double aprox;
    double error;
    double x_k_anterior = x0;
    double x_k_actual = x1;
    double x_k_siguiente;

    for (int k = 0; k < iterMax; k++) {

        x_k_siguiente = x_k_actual - ((x_k_actual - x_k_anterior)*funcion(x_k_actual)) /
                        (funcion(x_k_actual) - funcion(x_k_anterior));

        aprox = x_k_siguiente;

        error = abs(funcion(x_k_siguiente));

        iterVector.push_back(k);

        errorVector.push_back(error);

        if (error < tol) { // Condicion de parada del ciclo iterativo

            break;

        }

        x_k_actual = x_k_siguiente;

    }

    // Impresion de resultados de Iteraciones vs Error Absoluto del Metodo de la Secante
    cout << "Iteraciones = " << iterVector.back() << endl;
    cout << "Aproximacion = " << x_k_siguiente << endl;
    cout << "Error absoluto = " << error << endl;
    
    // Configuracion del grafico 
    plt::named_plot("Error Absoluto", iterVector, errorVector);
    plt::title("Grafico de Iteraciones vs Error Absoluto del Metodo de la Secante");
    plt::legend();
    plt::show();

    return x_k_siguiente;

}

int main() {

    // Ejemplo numerico
    double x0 = 0;
    double x1 = 1;
    double tol = pow(10,-3);
    int iterMax = 1000;

    double aprox = secante(x0, x1, tol, iterMax);

    return 0;
}