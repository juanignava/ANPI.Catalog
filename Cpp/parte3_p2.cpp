#include <iostream>
#include <armadillo>
#include <vector>

using namespace std;
using namespace arma;

mat pseudoinversa(mat A, mat b, double alpha_1, double alpha_2, double tol, int iterMax) {

    /*
    Esta funcion realiza una aproximacion iterativa de la solucion del sistema de ecuaciones
    A*x = b donde A es una que no necesariamente posee inversa por lo que se utiliza el 
    calculo de la pseudoinversa de la matriz y el resultado corresponde al valor mas cercano
    a una solucion.

    Parametros de entrada:
            A: es una matriz de cualquier dimension a la cual se le calculara la pseudoinversa
            alpha_1, alpha_2: valores numericos de las constantes a utilizar en el metodo iterativo
                deben ser positivos.
            tol: valor numerico de la tolerancia utilizada para la condicion de parada.
            iterMax: valor entero de la cantidad de iteraciones maximas del metodo.

    Parametros de salida:
            xk: Corresponde a la matriz que produce el error del sistema mas pequeno posible (aproximado
            a una solucion del sistema).
    */

    // validacion de constantes alpha
    if (alpha_1 <= 0 || alpha_2 <= 0)
    {
        cout << "Los valores de alpha deben ser positivos" << endl;
        return A;
    }
    
    // definicion de valores iniciales de la iteracion
    mat X0 = alpha_1 * A.t();
    mat X1 = alpha_2 * A.t();
    mat Xk = X1;
    mat Xk_prev = X0;
    double error = tol+1;

    mat xk = Xk * b;
    mat xk_n;

    // iteraciones para calcular la pseudoinversa
    for (int k = 0; k < iterMax; k++)
    {

        mat X_temp = Xk;
        Xk = Xk_prev + Xk - Xk_prev * A * Xk;
        Xk_prev = X_temp;

        // actualizacion xk
        xk_n = Xk * b;

        // calculo del error
        double error = norm(xk_n - xk) / norm(xk_n);

        xk = xk_n;

        // condicion de parada
        if (error < tol) {
            cout << "Error asociado: " << error << endl;
            break;
        }
        
    }
    
    return xk;
    
}

int main(){

    /*
    Ejemplo numerico, punto mas cercano a los puntos:
    (0, 0)
    (1, 3)
    (1.5, 1)
    */
   mat A = {{0.667,-1},
            {4,1},
            {3,-1}};

    mat b = {0,7,0};
    b = b.t();

    double alpha_1 = 5 * pow(10, -10);
    double alpha_2 = 2* pow(10, -11);
    double tol = pow(10, -5);
    int iterMax = 200;

    mat x = pseudoinversa(A, b, alpha_1, alpha_2, tol, iterMax);

    cout << "solucion al sistema" << endl;
    cout << x << endl;

    return 0;
}