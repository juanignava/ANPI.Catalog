#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main(){

  mat A(6, 5, fill::randu);
  mat B(6, 5, fill::randu);
  
  cout << A*B.t() << endl;
  
  return 0;

}