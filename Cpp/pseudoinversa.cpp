#include <iostream>
#include <armadillo>
#include <vector>

using namespace std;
using namespace arma;

mat pseudoinversa(mat A, mat b, double tol, int iterMax) {

    /*
    Esta  funcion  aproxima  numericamente, en forma iterativa, la solucion
    de un sistema de la forma Ax = b, por medio del metodo de la pseudoinversa
    basado en el metodo de Newton-Schultz.

    Parametros de entrada: A es un matriz pseudoinversa de cualquier dimension.

                           b corresponde a un vector columna que contiene los
                           resultados del sistema Ax = b.
    
                           tol hace refencia a la tolerancia del metodo iterativo.
                        
                           iterMax representa iteraciones maximas que se le
                           aplicaran al metodo iterativo.

    Parametros de salida: y es el vector que contiene la aproximacion de la solucion
                          del sistema.
    */

    int n = A.n_rows;

    mat Xk = (1 / (norm(A) * norm(A))) * A.t();

    mat xk = Xk * b;

    mat xk_n;

    mat identidad = eye(n, n);

    int iteracion;

    vector < int > iteraciones;

    double error;
        
    vector < double > errores;

    for(int k = 0; k < iterMax; k++) {

        iteracion = k;

        iteraciones.push_back(iteracion);

        Xk = Xk * (2 * identidad - A * Xk);

        xk_n = Xk * b;

        error = norm(xk_n - xk) / norm(xk_n);

        errores.push_back(error);   

        xk = xk_n;     

        if (error < tol) { // condicion de parada

            break;

        }       

    }

    cout << "Xk = " << Xk  << endl;
    cout << "k = " << iteracion << endl;
    cout << "error = " << error << endl;

    // instrucciones de graficacion
    /*
    plt::named_plot("Grafico de iteraciones vs error", iteraciones, errores);
    plt::title("Metodo de la Pseudoinversa");
    plt::legend();
    plt::show();
    */
    
    return xk;

}

int main(){

    // Ejemplo numerico
    /*
    mat A = {{0, 1, 4},
             {1, 2, 5},
             {4, 5, 8},
             {9, 10, 13},
             {16, 17, 20}};
    */

   mat A = {{0.667,-1},
            {4,1},
            {3,-1}};

    mat b = {0,7,0};
    b = b.t();

    double tol = 0.0000001;
    int iterMax = 1000;

    mat x = pseudoinversa(A, b, tol, iterMax);

    cout << "solucion al sistema" << x << endl;

    return 0;
}