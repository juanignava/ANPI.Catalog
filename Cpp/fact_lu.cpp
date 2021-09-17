#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

vector <mat> calculo_lu(mat A){
    /*
    Funcion auxiliar de la funcion fact_lu que calcula las matrices L y U
    de la respectiva factorizacion.

    Parametros de entrada: A: La amtriz de coeficientes del sistema de
                            ecuaciones.

    Parametros de salida: l_u: un vector que contiene a:
                            L -> La matriz triangular inferior.
                            U -> La matriz triangular superior.
                          Tales que LU = A.
    */
    int m = A.n_rows;
    mat L = eye<mat> (m, m); // matriz identidad
    mat U = zeros<mat> (m, m); // matriz cero

    // por cada elemento en la diagonal
    for (int j = 0; j < m-1; j++) 
    {
        double diag_el = A(j, j);

        // por cada fila en la matriz
        for (int fila = j+1; fila < m; fila++) 
        {
            // multiplo que se necesita en las operaciones
            double multiplo = A(fila, j)/diag_el;

            // este multiplo corresponde a un termino de J
            L(fila, j) = multiplo;

            // por cada columna en la matriz
            for (int column = 0; column < m; column++)
            {
                // operacion de la fila para igualar a cero
                // los terminos abajo de la diagonal.
                A(fila, column) -= multiplo*A(j, column);
            }            
        }
        
    }

    vector <mat> l_u = {L, A};
    return l_u;
    
}

vec sust_adelante(mat A, mat b){
    /*
    Funcion auxiliar de la funcion fact_lu que encuantra la solucion de 
    un sistema triengular inferior por medio de la sustitucion hacia 
    adelante.
    Parametros de entrada: A: una matriz triangular inferior.
                           b: un vetor columna con los terminos 
                            independientes.

    Parametros de salida: x: un vector columna con la solucion del 
                            sistema descrito.
    */
    int m=A.n_rows;
    mat x=zeros<mat>(m,1);

    for (int i = 0; i<m;i++){
        double aux=0;

        for (int j=0; j<=i-1;j++){
            aux+=A(i,j)*x(j);
        }

        x(i)=(1/A(i,i))*(b(i)-aux);
    }

    return x;
}

vec sust_atras(mat A, mat b){
    /*
    Funcion auxiliar de la funcion fact_lu que encuantra la solucion de 
    un sistema triengular superior por medio de la sustitucion hacia 
    atras.
    Parametros de entrada: A: una matriz triangular superior.
                           b: un vetor columna con los terminos 
                            independientes.

    Parametros de salida: x: un vector columna con la solucion del 
                            sistema descrito.
    */
    int m=A.n_rows;
    mat x=zeros<mat>(m,1); 

    for (int i = m-1; i>=0;i--){
        double aux=0;

        for (int j=i+1; j<=m-1;j++){
            aux+=A(i,j)*x(j);
        }

        x(i)=(1/A(i,i))*(b(i)-aux);
    }

    return x;
}

int sub_princ(mat A, double tol){
    /*
    Funcion auxiliar de la funcion fact_lu que verifica que las submatrices
    principales de A sean invertibles.
    Parametros de entrada: A: Una matriz cuadrada.
                           tol: la tolerancia del determinante
    Parametros de salida: un entero 1 o 0 que indica si la matriz cumple la
                        condicion (1) o si no la cumple (0).
    */
    int m = A.n_rows;

    for (int i = 1; i <= m; i++)
    {
        mat Ak = A;
        Ak.resize(i, i);
        double deter = abs(arma::det(Ak));

        // Verifica que el determinante no sea cero.
        if (deter < tol){
            return 0;
        }
    }

    return 1;
}

mat fact_lu(mat A, mat b) {
    /*
    Esta funcion encuantra la solucion exacta el sistema de ecuaciones cuya
    matriz de coeficientes es A y su vector de resultados es b. El calculo
    se realiza mediante el metodo de factorizacion LU

    Parametros de entrada: A: corresponde a una matriz cuadrada invertible 
                            y con cada una de sus submatrices principales
                            invertibles.

                           b: corresponde a un vector columna con los 
                           terminos independientes.

    Parametros de salida: x: corresponde a la solucion de la matriz de 
                            incognitas del sistema.
    */

    vector <mat> l_u = calculo_lu(A); // retorna la matriz L y U

    // Se guardan las matrices L y U
    mat L = l_u[0]; 
    mat U = l_u[1];

    // Se resuelve el sistema Ly=b
    mat y = sust_adelante(L, b);

    // Se resulve el sistema Ux=y
    mat x = sust_atras(U, y);

    return x;

}

int main(){

    // Ejemplo numerico
    mat A = {{4, -2, 1}, {20, -7, 12}, {-8, 13, 17}};
    mat b = {11, 70, 17};
    b = b.t();

    // Verificacion de condicion de la matriz
    // La verificacion es que todas sus submatrices principales
    // sean invertibles.
    double det_tol = pow(10, -5);
    int ver = sub_princ(A, det_tol);

    if(ver == 0){
        cout << "La matriz no tiene una unica factorizacion LU" << endl;
        return 0;
    }

    // Llamado del metodo de factorizacion LU
    mat x = fact_lu(A, b);
    cout << "el vector x" << x << endl;

    return 0;
}