#include <iostream>
#include <armadillo>
#include <vector>

using namespace std;
using namespace arma;

mat generar_matriz(int rows, int columns) {

    mat Result = zeros(rows, columns);

    
    for (int i = 0; i < rows; i++)
    {
        
        for (int j = 0; j < columns; j++)
        {
            Result(i, j) = pow(i, 2) + pow(j, 2);

        }
        
    }
     

    return Result;

}

double norma_frobenius(mat A) {

    int n = A.n_rows;
    int m = A.n_cols;

    double sum = 0;

    for (int i = 0; i < n; i++)
    {
        
        for (int j = 0; j < m; j++)
        {
            sum = sum + pow(A(i, j), 2);

        }
        
    }

    double root = sqrt(sum);

    return root;

}

mat pseudoinversa(mat A, double alpha_1, double alpha_2, double tol, int iterMax) {

    /*
    Esta funcion realiza una aproximacion iterativa de la pseudoinversa de una
    matriz A.

    Parametros de entrada:
            A: es una matriz de cualquier dimension a la cual se le calculara la pseudoinversa
            alpha_1, alpha_2: valores numericos de las constantes a utilizar en el metodo iterativo
                deben ser positivos.
            tol: valor numerico de la tolerancia utilizada para la condicion de parada.
            iterMax: valor entero de la cantidad de iteraciones maximas del metodo.

    Parametros de salida:
            Xk: es la matriz aproximada de la pseudoinversa de la matriz A.
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

    // iteraciones para calcular la pseudoinversa
    for (int k = 0; k < iterMax; k++)
    {

        mat X_temp = Xk;
        Xk = Xk_prev + Xk - Xk_prev * A * Xk;
        Xk_prev = X_temp;

        // calculo del error con norma de frobenius
        mat error_mat = A * Xk * A - A;
        error = norma_frobenius(error_mat);

        // condicion de parada
        if (error < tol) {
            cout << "Error asociado: " << error << endl;
            break;
        }
        
    }
    
    return Xk;
    
}

int main() {

    // Ejemplo numerico
    mat A = generar_matriz(45, 30);

    double alpha_1 = 5 * pow(10, -10);
    double alpha_2 = 2* pow(10, -11);
    double tol = pow(10, -5);
    int iterMax = 200;

    mat A_pseudo = pseudoinversa(A, alpha_1, alpha_2, tol, iterMax);

    return 0;

}