#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

/// Función Sustitución Hacia Atrás (Matriz Triangular Superior)
vec sust_atras(mat A, mat b){
  int m=A.n_rows;
  mat x=zeros<mat>(m,1); // Crea matriz de m filas y 1 columna
  for (int i = m-1; i>=0;i--){
    double aux=0;
    for (int j=i+1; j<=m-1;j++){
      aux+=A(i,j)*x(j);
    }
    x(i)=(1/A(i,i))*(b(i)-aux);
  }
  return x;
}

int main()
  {
  mat A={{1,-3,6,8},{0, 1, 0,-2},{0,0,3,1},{0,0,0,1}}; // Matriz de 4x4
  mat b1={35,-9,5,5}; // Vector de 1 fila y 4 columnas
  mat b=b1.t(); // Transpuesta para obtener vector columna

  mat x = sust_atras(A,b);
  ///mat B=join_rows(A,b);


  cout << x << endl;
  //cout << B << endl;

  return 0;
  }